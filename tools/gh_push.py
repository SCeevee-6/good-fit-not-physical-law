"""Push a folder to a GitHub repo as one commit through the REST API (no local git needed).

usage: python tools/gh_push.py <repo-name> <folder> "<commit message>" [--private] [--description "text"]
token: env GITHUB_TOKEN. Skips whatever .gitignore in <folder> lists (simple patterns only).
Creates the repo if it does not exist. Prints one line: url, commit, file count.
"""
import os, sys, json, base64, fnmatch, hashlib, urllib.request, urllib.error

API = 'https://api.github.com'


def call(method, path, token, body=None):
    data = None if body is None else json.dumps(body).encode()
    req = urllib.request.Request(API + path, method=method, data=data, headers={
        'Authorization': 'Bearer ' + token, 'Accept': 'application/vnd.github+json',
        'X-GitHub-Api-Version': '2022-11-28', 'User-Agent': 'push-script'})
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            t = r.read()
            return r.status, (json.loads(t) if t else {})
    except urllib.error.HTTPError as e:
        t = e.read()
        try:
            return e.code, json.loads(t)
        except Exception:
            return e.code, {'raw': t[:200].decode('utf-8', 'replace')}


def ignored(rel, pats):
    parts = rel.split('/')
    for p in pats:
        p = p.strip().rstrip('/')
        if not p or p.startswith('#'):
            continue
        if fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(parts[-1], p) or p in parts or rel.startswith(p + '/'):
            return True
    return False


def collect(folder):
    pats = ['.git']
    gi = os.path.join(folder, '.gitignore')
    if os.path.exists(gi):
        pats += open(gi, encoding='utf-8').read().splitlines()
    out = []
    for root, dirs, files in os.walk(folder):
        rel_root = os.path.relpath(root, folder).replace('\\', '/')
        base = '' if rel_root == '.' else rel_root + '/'
        dirs[:] = [d for d in dirs if not ignored(base + d, pats)]
        for f in files:
            rel = base + f
            if not ignored(rel, pats):
                out.append(rel)
    return sorted(out)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    name, folder, msg = args[0], args[1], args[2]
    private = '--private' in sys.argv
    desc = ''
    if '--description' in sys.argv:
        desc = sys.argv[sys.argv.index('--description') + 1]
        if desc in args:
            args.remove(desc)
    token = os.environ['GITHUB_TOKEN']
    code, me = call('GET', '/user', token)
    assert code == 200, ('auth failed', code)
    login = me['login']
    code, repo = call('GET', f'/repos/{login}/{name}', token)
    if code == 404:
        code, repo = call('POST', '/user/repos', token, {'name': name, 'description': desc, 'private': private, 'auto_init': True})
        assert code in (200, 201), ('create failed', code, repo.get('message'))
    branch = repo['default_branch']
    code, ref = call('GET', f'/repos/{login}/{name}/git/ref/heads/{branch}', token)
    assert code == 200, ('no branch', code)
    head = ref['object']['sha']
    code, cur = call('GET', f'/repos/{login}/{name}/git/trees/{head}?recursive=1', token)
    have = {e['path']: e['sha'] for e in cur.get('tree', []) if e.get('type') == 'blob'} if code == 200 else {}
    entries, sent = [], 0
    for rel in collect(folder):
        with open(os.path.join(folder, rel), 'rb') as fh:
            data = fh.read()
        sha = hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()   # git blob id; unchanged files are not re-sent
        if have.get(rel) != sha:
            code, blob = call('POST', f'/repos/{login}/{name}/git/blobs', token, {'content': base64.b64encode(data).decode(), 'encoding': 'base64'})
            assert code == 201, ('blob failed', rel, code)
            sha = blob['sha']; sent += 1
        entries.append({'path': rel, 'mode': '100644', 'type': 'blob', 'sha': sha})
    tree = None                                   # full snapshot, built in chunks (one large tree call can 502)
    for i in range(0, len(entries), 300):
        body = {'tree': entries[i:i + 300]}
        if tree:
            body['base_tree'] = tree['sha']
        code, tree = call('POST', f'/repos/{login}/{name}/git/trees', token, body)
        assert code == 201, ('tree failed', code, tree.get('message'))
    code, commit = call('POST', f'/repos/{login}/{name}/git/commits', token, {'message': msg, 'tree': tree['sha'], 'parents': [head]})
    assert code == 201, ('commit failed', code)
    code, _ = call('PATCH', f'/repos/{login}/{name}/git/refs/heads/{branch}', token, {'sha': commit['sha']})
    assert code == 200, ('ref update failed', code)
    print(repo['html_url'], commit['sha'][:8], len(entries), 'files', sent, 'sent')


if __name__ == '__main__':
    main()
