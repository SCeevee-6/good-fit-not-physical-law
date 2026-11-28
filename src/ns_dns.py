"""Protocol S1: 2D incompressible Navier-Stokes in vorticity form on [0, 2pi]^2, pseudo-spectral (docs/protocol_S1.md).

    w_t + u . grad w = nu lap w + F,  F = sin(4 y),  u = (psi_y, -psi_x),  lap psi = -w.
Grid N x N, 2/3 dealiasing, RK4 with fixed dt. Arrays are indexed [y, x]; rfft2 gives (ky, kx >= 0).
F = sin(4y) is the vorticity forcing, i.e. the velocity forcing f = (cos(4y)/4, 0) with curl f = F, so the energy
input is P = <u_x cos(4y)> / 4 and the dissipation is D = nu <w^2>; for the dealiased system dE/dt = P - D exactly.
"""
import numpy as np


class Grid:
    def __init__(self, N=64, nu=0.02):
        self.N, self.nu = N, nu
        k = np.fft.fftfreq(N, 1.0 / N); kr = np.fft.rfftfreq(N, 1.0 / N)
        self.KY, self.KX = np.meshgrid(k, kr, indexing='ij')
        self.K2 = self.KX ** 2 + self.KY ** 2; self.K2i = np.where(self.K2 > 0, 1.0 / np.where(self.K2 > 0, self.K2, 1), 0.0)
        cut = N / 3.0
        self.mask = (np.abs(self.KX) < cut) & (np.abs(self.KY) < cut)
        x = np.arange(N) * 2 * np.pi / N
        self.Y, self.X = np.meshgrid(x, x, indexing='ij')
        self.Fh = np.fft.rfft2(np.sin(4 * self.Y)) * self.mask
        self.fx = np.cos(4 * self.Y) / 4.0

    def fields(self, wh):
        """u, v, w, w_x, w_y in physical space from vorticity coefficients."""
        ph = wh * self.K2i
        u = np.fft.irfft2(1j * self.KY * ph, s=(self.N, self.N)); v = np.fft.irfft2(-1j * self.KX * ph, s=(self.N, self.N))
        wx = np.fft.irfft2(1j * self.KX * wh, s=(self.N, self.N)); wy = np.fft.irfft2(1j * self.KY * wh, s=(self.N, self.N))
        return u, v, wx, wy

    def nonlinear(self, wh):
        """Dealiased transform of u . grad w."""
        u, v, wx, wy = self.fields(wh)
        return np.fft.rfft2(u * wx + v * wy) * self.mask

    def rhs(self, wh):
        return -self.nonlinear(wh) - self.nu * self.K2 * wh + self.Fh

    def diagnostics(self, wh):
        u, v, wx, wy = self.fields(wh)
        w = np.fft.irfft2(wh, s=(self.N, self.N))
        E = 0.5 * float(np.mean(u ** 2 + v ** 2)); P = float(np.mean(u * self.fx)); D = self.nu * float(np.mean(w ** 2))
        return E, P, D


def initial(grid, seed, amp=0.1, kmax=8):
    rng = np.random.default_rng([seed, 4242])
    wh = (rng.normal(size=grid.K2.shape) + 1j * rng.normal(size=grid.K2.shape)) * ((grid.K2 > 0) & (grid.K2 <= kmax ** 2))
    w = np.fft.irfft2(wh, s=(grid.N, grid.N)); w *= amp / max(np.std(w), 1e-300)
    return np.fft.rfft2(w) * grid.mask


def rk4(grid, wh, dt):
    k1 = grid.rhs(wh); k2 = grid.rhs(wh + 0.5 * dt * k1); k3 = grid.rhs(wh + 0.5 * dt * k2); k4 = grid.rhs(wh + dt * k3)
    return wh + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def run(grid, wh, dt, nsteps, every, on_sample=None):
    """Advance nsteps. Budget terms at every step (trapezoid of P - D against the change in E); on_sample(t_index, wh)
    every `every` steps. Returns arrays E, P, D at every step (length nsteps + 1) and the final state."""
    E = np.empty(nsteps + 1); P = np.empty(nsteps + 1); D = np.empty(nsteps + 1)
    E[0], P[0], D[0] = grid.diagnostics(wh)
    for i in range(1, nsteps + 1):
        wh = rk4(grid, wh, dt)
        E[i], P[i], D[i] = grid.diagnostics(wh)
        if on_sample is not None and i % every == 0:
            on_sample(i, wh)
    return E, P, D, wh
