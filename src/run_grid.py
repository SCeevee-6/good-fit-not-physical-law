"""Run a list of grid cells in parallel (multiprocessing spawn). Usage: python run_grid.py <workers> <key:noise_idx:seed> ..."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import multiprocessing as mp


def main():
    nw = int(sys.argv[1]); specs = sys.argv[2:]
    outdir = os.path.join(HERE, '..', 'results', 'cells'); os.makedirs(outdir, exist_ok=True)
    logpath = os.path.join(HERE, '..', 'logs', 'grid.log')
    jobs = []
    for s in specs:
        key, n, sd = s.split(':'); jobs.append((key, int(n), int(sd), outdir, logpath))
    import run_cell
    with mp.get_context('spawn').Pool(nw) as pl:
        for cid, st in pl.imap_unordered(run_cell.run_cell, jobs):
            print(cid, st, flush=True)


if __name__ == '__main__':
    main()
