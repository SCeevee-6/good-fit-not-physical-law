"""Numpy copy of scipy.signal.savgol_filter (mode='interp', axis=0) for machines where scipy.signal will not load.
Interior points use the least-squares convolution coefficients; the first and last window//2 points use a
polynomial fit to the edge window, as scipy does. Only the interior is used downstream (data.smooth trims the edges).
"""
from math import factorial
import numpy as np


def savgol_coeffs(window_length, polyorder, deriv=0, delta=1.0):
    half = window_length // 2
    x = np.arange(-half, window_length - half, dtype=float)
    A = x[None, :] ** np.arange(polyorder + 1)[:, None]
    y = np.zeros(polyorder + 1)
    y[deriv] = factorial(deriv) / delta ** deriv
    return np.linalg.lstsq(A, y, rcond=None)[0]          # 'dot' order: out[n] = sum_i c[i] x[n - half + i]


def _edge(x, w0, w1, i0, i1, polyorder, deriv, delta, out):
    t = np.arange(w1 - w0, dtype=float)
    P = np.polyfit(t, x[w0:w1], polyorder)
    if deriv:
        P = np.array([np.polyder(P[:, k], deriv) for k in range(P.shape[1])]).T if P.ndim == 2 else np.polyder(P, deriv)
    tt = np.arange(i0 - w0, i1 - w0, dtype=float)
    V = np.vander(tt, P.shape[0]) @ P
    out[i0:i1] = V / delta ** deriv


def savgol_filter(x, window_length, polyorder, deriv=0, delta=1.0, axis=0, mode='interp'):
    assert axis == 0 and mode == 'interp'
    x = np.asarray(x, dtype=float)
    squeeze = x.ndim == 1
    if squeeze:
        x = x[:, None]
    n, w, half = x.shape[0], window_length, window_length // 2
    c = savgol_coeffs(w, polyorder, deriv, delta)
    out = np.zeros_like(x)
    for i in range(w):
        out[half:n - half] += c[i] * x[i:n - w + 1 + i]
    _edge(x, 0, w, 0, half, polyorder, deriv, delta, out)
    _edge(x, n - w, n, n - half, n, polyorder, deriv, delta, out)
    return out[:, 0] if squeeze else out
