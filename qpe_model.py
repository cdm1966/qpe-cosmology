"""
qpe_model.py --- Reference implementation of the QPE distance-redshift relation
used in Chen (2026), "Cosmological Tests of Scale-Dependent Quantum Potential
Energy".

Everything in this module is a direct transcription of the equations of the
paper; the equation numbers quoted in the comments refer to the published
version.  No data are hard-coded except the thirteen DESI DR2 measurements,
which are reproduced in `DESI_DR2` below.

The two ingredients are:

  1. Geometry (Section 3.4).  The metric is static,

         ds^2 = g1^2 dt^2 - g1^{-2} dr^2 - C(r)^2 dOmega^2 ,
         g1(x) = sqrt(1 - 2 sinh x),   x = H0 r / c ,

     with 1 + z = 1/g1.  The scale dependence is carried by the areal radius

         C(r) = [1+z]^{gamma(z)-2} D_prop(z) ,                        Eq. (86)

     so that  D_A = C ,  D_L = (1+z)^gamma D_prop ,  D_M = (1+z)^{gamma-1} D_prop .
     The exponent is three constant segments, Eq. (90),

         gamma(z) = g1 + (g2-g1)/2 [1+tanh((z-z1)/dz1)]
                       + (g3-g2)/2 [1+tanh((z-z2)/dz2)] .

  2. CMB projection (Section 6).  The emission kernel is

         W(x) ~ u^{1-p/2} exp[-tau(x)] ,   u = g1^2 ,                 Eq. (148)

     with  tau = tau2 * A(x)/A(x_*) ,  A = int u^{-(q+1)/2} dx  and x_* at z=2.
     The angular power spectrum is Eq. (156),

         C_l = (2/pi) int k^2 dk P(k) [ int_0^{Cmax} G(u) j_l(ku) du ]^2 ,

     where the change of variable u = C(r) is two-valued (C is not monotonic)
     and G(u) sums the two branches, Eq. (158).

A note on units: distances are in Mpc, H0 in km/s/Mpc, k in 1/Mpc.  Internally
the projection uses dimensionless x = H0 r / c and k' = k c / H0.
"""

import numpy as np
from scipy import integrate
from scipy.special import spherical_jn

# --------------------------------------------------------------------------
# Constants and published parameter values
# --------------------------------------------------------------------------

C_LIGHT = 299792.458                  # speed of light, km/s
MAX_X = np.arcsinh(0.5)               # horizon, x_h = ln(phi) = 0.4812

#: Three-segment exponent, Eq. (91) of the paper.
GAMMA_1, GAMMA_2, GAMMA_3 = 1.1886, 1.0912, 0.9409

#: Transition redshifts and widths of Eq. (90).
Z1, DZ1, Z2, DZ2 = 1.5, 0.10, 100.0, 15.0

#: SN Ia fit (Section 4.2), chi2 = 745.6 for N = 1701.
SN_PARAMS = dict(H0=72.664, gamma1=1.1888, gamma2=1.0935, gamma3=GAMMA_3)

#: BAO fit (Section 7.3), chi2 = 8.0 for 13 observables.
BAO_PARAMS = dict(H0=72.66, rd=132.3, gamma_M=(2.246, 2.115),
                  gamma_H=(1.415, 0.303), z_split=1.0, width=0.15)

#: CMB fit (Section 6.7), chi2 = 1464.5 for N = 836, thirteen parameters
#: including the cutoff exponent beta of Eq. (108).
CMB_PARAMS = dict(A_pk=5.8245e6, H0=73.048, A0_low=5911.0, A0_high=1124.0,
                  n_low=2.6091, n_mid=4.0714, n_high=3.8059,
                  k_omega=0.059745, k_break=1.4397, p=3.0874, q=17.173,
                  tau2=164.05, beta=5.81)

#: DESI DR2 BAO measurements (Section 7.2).
#: columns: z, observable (1 = D_M/r_d, 2 = D_H/r_d, 0 = D_V/r_d), value, sigma
DESI_DR2 = np.array([
    [0.295, 0, 7.944, 0.075],
    [0.510, 1, 13.587, 0.169], [0.510, 2, 21.863, 0.427],
    [0.706, 1, 17.347, 0.180], [0.706, 2, 19.458, 0.332],
    [0.934, 1, 21.574, 0.153], [0.934, 2, 17.641, 0.193],
    [1.321, 1, 27.605, 0.320], [1.321, 2, 14.178, 0.217],
    [1.484, 1, 30.519, 0.758], [1.484, 2, 12.816, 0.513],
    [2.330, 1, 38.988, 0.531], [2.330, 2, 8.632, 0.101]])


# --------------------------------------------------------------------------
# Geometry: gamma(z), D_prop, C(r) and D_L  (Section 3)
# --------------------------------------------------------------------------

def gamma_of_z(z, g1=GAMMA_1, g2=GAMMA_2, g3=GAMMA_3):
    """Three-segment exponent gamma(z), Eq. (90)."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return (g1 + (g2 - g1) * 0.5 * (1.0 + np.tanh((z - Z1) / DZ1))
            + (g3 - g2) * 0.5 * (1.0 + np.tanh((z - Z2) / DZ2)))


def x_of_z(z):
    """Dimensionless radial coordinate, from 1+z = 1/sqrt(1-2 sinh x)."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    arg = np.clip(0.5 * (1.0 - (1.0 + z) ** -2), 0.0, np.sinh(MAX_X) - 1e-12)
    return np.arcsinh(arg)


def z_of_x(x):
    """Inverse of x_of_z."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    u = np.clip(1.0 - 2.0 * np.sinh(np.clip(x, 0.0, MAX_X)), 1e-30, None)
    return 1.0 / np.sqrt(u) - 1.0


def d_prop_of_z(z, H0, nx=20000):
    """Light-travel distance D_prop = (c/H0) int dx/u, Eq. (79).  [Mpc]"""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    xg = np.linspace(0.0, MAX_X * 0.99999, nx)
    ug = np.clip(1.0 - 2.0 * np.sinh(xg), 1e-30, None)
    dp = integrate.cumulative_trapezoid(1.0 / ug, xg, initial=0.0)
    return np.interp(x_of_z(z), xg, dp) * C_LIGHT / H0


def C_of_z(z, H0, g1=GAMMA_1, g2=GAMMA_2, g3=GAMMA_3):
    """Areal radius C(r), Eq. (86).  [Mpc]"""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return (1.0 + z) ** (gamma_of_z(z, g1, g2, g3) - 2.0) * d_prop_of_z(z, H0)


def D_L_of_z(z, H0, g1=GAMMA_1, g2=GAMMA_2, g3=GAMMA_3):
    """Luminosity distance D_L = (1+z)^gamma D_prop, Eq. (89).  [Mpc]"""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return (1.0 + z) ** gamma_of_z(z, g1, g2, g3) * d_prop_of_z(z, H0)


def mu_of_z(z, H0, g1=GAMMA_1, g2=GAMMA_2, g3=GAMMA_3):
    """Distance modulus, Eq. (84)."""
    return 5.0 * np.log10(D_L_of_z(z, H0, g1, g2, g3)) + 25.0


def K_gamma(x_query, gamma_lo, gamma_hi, z_split=1.0, width=0.15, nx=4000):
    """Cumulative path integral K_gamma(x) of Eq. (162).  [dimensionless]

    The exponent is gamma_lo below z_split and gamma_hi above it, joined by a
    tanh transition of the given width in ln z.
    """
    xg = np.linspace(1e-8, MAX_X * 0.99999, nx)
    u = np.clip(1.0 - 2.0 * np.sinh(xg), 1e-30, None)
    zg = 1.0 / np.sqrt(u) - 1.0
    t = 0.5 * (1.0 + np.tanh((np.log(np.maximum(zg, 1e-8)) - np.log(z_split)) / width))
    gam = gamma_lo + (gamma_hi - gamma_lo) * t
    integrand = u ** (-0.5 * gam)
    K = integrate.cumulative_trapezoid(integrand, xg, initial=0.0)
    return np.interp(np.atleast_1d(x_query), xg, K)


# --------------------------------------------------------------------------
# Type Ia supernovae (Section 4)
# --------------------------------------------------------------------------

def sn_chi2(data_file, params=None, z_col="zHD", mu_col="MU_SH0ES",
            err_col="MU_SH0ES_ERR_DIAG"):
    """chi^2 of the SN Ia Hubble diagram, Eq. (94).

    `data_file` is the Pantheon+SH0ES catalogue (not redistributed here; see
    README.md for the download instructions).
    """
    p = params or SN_PARAMS
    d = np.genfromtxt(data_file, names=True)
    z = d[z_col]
    mu = d[mu_col]
    err = d[err_col]
    m = np.isfinite(z) & np.isfinite(mu) & np.isfinite(err) & (err > 0) & (z > 0)
    z, mu, err = z[m], mu[m], err[m]
    model = mu_of_z(z, p["H0"], p["gamma1"], p["gamma2"], p["gamma3"])
    return float(np.sum(((mu - model) / err) ** 2)), len(z)


# --------------------------------------------------------------------------
# Baryon acoustic oscillations (Section 7)
# --------------------------------------------------------------------------

def boltzmann_lcdm(z, H0=67.4, om=0.315):
    """Flat-LCDM distances for the comparison of Fig. 5.  [Mpc]"""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    zz = np.linspace(0.0, max(z.max() * 1.2, 3.0), 5000)
    E = np.sqrt(om * (1.0 + zz) ** 3 + (1.0 - om))
    chi = integrate.cumulative_trapezoid(1.0 / E, zz, initial=0.0) * (C_LIGHT / H0)
    DM = np.interp(z, zz, chi)
    DH = C_LIGHT / (H0 * np.sqrt(om * (1.0 + z) ** 3 + (1.0 - om)))
    return DM, DH


def bao_model(z, params=None):
    """(D_M/r_d, D_H/r_d) of Eqs. (152)-(154)."""
    p = params or BAO_PARAMS
    x = x_of_z(z)
    KM = K_gamma(x, p["gamma_M"][0], p["gamma_M"][1], p["z_split"], p["width"])
    KH = K_gamma(x, p["gamma_H"][0], p["gamma_H"][1], p["z_split"], p["width"])
    DM = (C_LIGHT / p["H0"]) * KM / p["rd"]
    DH = (C_LIGHT / p["H0"]) * KH / (np.atleast_1d(z) * p["rd"])
    return DM, DH


def bao_chi2(params=None, data=None):
    """chi^2 of the BAO measurements, Eq. (166), diagonal uncertainties."""
    p = params or BAO_PARAMS
    data = DESI_DR2 if data is None else data
    z, kind, val, sig = data[:, 0], data[:, 1].astype(int), data[:, 2], data[:, 3]
    DM, DH = bao_model(z, p)
    pred = np.where(kind == 1, DM, np.where(kind == 2, DH,
                    np.cbrt(z * DM ** 2 * DH)))
    return float(np.sum(((val - pred) / sig) ** 2)), len(z)



# --------------------------------------------------------------------------
# Metric-perturbation kernel (Appendix B)
# --------------------------------------------------------------------------

def metric_kernel(xg, ug, a0, beta=5.81, use_cut=True):
    """W_met(x), Eq. (B28) of the paper, in the same units as W_em.

        W_met(x) = -x^2 int_x^{x_h} K(x') dx' / int_0^{x_h} g1^2 dx'
        K(x)     = 2 g1^2 G/(g1 x) [4 pi rho_c + 8 pi A0 rho_c sinh x]
                 = 3 g1 (1 + 2 A0 f sinh x)/x        (units H0 = 1)
    """
    g1 = np.sqrt(np.clip(ug, 1e-30, None))
    fac = (np.clip(1.0 - (np.clip(xg, 0.0, MAX_X) / MAX_X) ** beta, 0.0, None)
           if use_cut else 1.0)
    K = 3.0 * g1 * (1.0 + 2.0 * a0 * fac * np.sinh(xg)) / np.maximum(xg, 1e-12)
    IK = integrate.cumulative_trapezoid(K, xg, initial=0.0)
    D = integrate.trapezoid(ug, xg)
    return -xg ** 2 * (IK[-1] - IK) / D

# --------------------------------------------------------------------------
# CMB (Section 6)
# --------------------------------------------------------------------------

def build_geometry(nx=4000, gamma_fn=gamma_of_z):
    """Tabulate x, z, C(x) and dC/dx on a grid refined near the horizon.

    The grid is uniform at small x and logarithmic in (x_h - x) near the
    horizon, where both D_prop and C vary rapidly (Section 6.5).
    """
    x1 = np.linspace(1e-6, 0.85 * MAX_X, nx // 2)
    x2 = MAX_X - np.logspace(np.log10(0.15 * MAX_X), -14.0, nx // 2)
    xg = np.unique(np.concatenate([x1, x2]))
    ug = np.clip(1.0 - 2.0 * np.sinh(xg), 1e-30, None)
    zg = 1.0 / np.sqrt(ug) - 1.0
    dp = integrate.cumulative_trapezoid(1.0 / ug, xg, initial=0.0)
    Cx = (1.0 + zg) ** (gamma_fn(zg) - 2.0) * dp
    i_turn = int(np.argmax(Cx))
    dCdx = np.gradient(Cx, xg)
    return xg, ug, zg, Cx, dCdx, i_turn


def emission_kernel(xg, ug, x_star, p, q, tau2):
    """Emission kernel W(x), Eq. (148).

    tau is normalised at x_*, the radius at which z = 2, so that
    tau(x_*) = tau2 and the optical depth grows towards the horizon.
    """
    integ = ug ** (-0.5 * (q + 1.0))
    A = integrate.cumulative_trapezoid(integ, xg, initial=0.0)
    A_star = float(np.interp(x_star, xg, A))
    tau = tau2 * A / max(A_star, 1e-30)
    W = ug ** (1.0 - 0.5 * p) * np.exp(-np.clip(tau, 0.0, 80.0))
    return W


def n_of_k(k, n_low, n_mid, n_high, k_omega, sigma=0.30):
    """Three-segment spectral index, Section 6.6 (transitions at k_w and 3k_w)."""
    k = np.maximum(np.asarray(k, dtype=float), 1e-30)
    k1, k2 = max(k_omega, 1e-6), max(3.0 * k_omega, 1e-6)
    t1 = 0.5 * (1.0 + np.tanh((np.log(k) - np.log(k1)) / sigma))
    t2 = 0.5 * (1.0 + np.tanh((np.log(k) - np.log(k2)) / sigma))
    return n_low + (n_mid - n_low) * t1 + (n_high - n_mid) * t2


def cutoff_volume(x):
    """f(x) = 1 - (x/x_h)^3: the local quantum potential is switched off as the
    volume of the perturbation approaches that of the causal horizon, so that
    the metric of a local system rejoins the cosmological metric at infinity."""
    return np.clip(1.0 - (np.asarray(x, dtype=float) / MAX_X) ** 3, 0.0, None)


def cutoff_none(x):
    return np.ones_like(np.asarray(x, dtype=float))


def cutoff_horizon(x, beta=None):
    """f(x) = 1 - (x/x_h)^beta, Eq. (108) of the paper."""
    b = 5.81 if beta is None else beta
    return np.clip(1.0 - (np.clip(np.asarray(x, dtype=float), 0.0, MAX_X) / MAX_X) ** b,
                   0.0, None)


def power_spectrum(k, p, cutoff=None):
    """Matter power spectrum P(k), Eq. (131) of the paper.  [Mpc^3]

    `cutoff` is an optional callable f(x) multiplying A_0.  The default (None)
    reproduces the fit with a constant A_0, which is what Section 6.7 of the
    paper reports as the no-cutoff comparison; pass `cutoff_horizon` to
    include the horizon cutoff with the fitted beta.
    """
    k = np.maximum(np.asarray(k, dtype=float), 1e-30)
    nk = n_of_k(k, p["n_low"], p["n_mid"], p["n_high"], p["k_omega"])
    cos2 = np.cos(k / p["k_omega"]) ** 2
    x_of_k = p["H0"] / (np.maximum(k, p["H0"] / (0.481 * C_LIGHT)) * C_LIGHT)
    a_t = 0.5 * (1.0 + np.tanh((np.log(k) - np.log(p["k_omega"])) / 0.30))
    A0k = p["A0_low"] + (p["A0_high"] - p["A0_low"]) * a_t
    if cutoff is not None:
        A0k = A0k * cutoff(x_of_k)
    om2 = np.maximum(2.0 * A0k * np.sinh(np.clip(x_of_k, 0.0, 20.0)) - 1.0, 0.05)
    P_low = p["A_pk"] * k ** (-nk) / np.sqrt(om2) * cos2
    return np.where(k < p["k_break"], P_low,
                    p["A_pk"] * k ** (-nk) * np.exp(-k / p["k_break"]))


class CMBProjector:
    """Two-branch projection of Eq. (156)-(150) on the Planck multipoles."""

    def __init__(self, data_file, nu=900, nk=300, ell_subsample=3,
                 ell_max=2508.0, nx=4000, nproc=1):
        self.nu, self.nk = int(nu), int(nk)
        raw = np.loadtxt(data_file)
        ell = raw[:, 0].astype(int)
        Dl = raw[:, 1]
        err = 0.5 * (raw[:, 2] + raw[:, 3])
        m = ell <= ell_max
        self.ells = ell[m][::ell_subsample]
        self.Dl, self.err = Dl[m][::ell_subsample], err[m][::ell_subsample]

        self.XG, self.UG, self.ZG, self.CX, self.DCDX, self.I = build_geometry(nx)
        self.CMAX = self.CX[self.I]
        self.x_star = float(np.interp(2.0, self.ZG, self.XG))

        # uniform grid in the areal radius u, Eq. (158)
        self.U = np.linspace(0.0, self.CMAX * 0.9997, self.nu)
        # Branch 1: C increases from 0 to Cmax.  Branch 2: C decreases from
        # Cmax back to 0, and is tabulated in reverse so that the abscissa is
        # ascending, as np.interp requires.  The explicit left/right values fix
        # the behaviour for u below the value of C at the horizon, where the
        # second branch has no resolved solution; they reproduce the choice
        # made in the published analysis.
        self.XA = np.interp(self.U, self.CX[:self.I + 1], self.XG[:self.I + 1],
                            left=self.XG[0], right=self.XG[self.I])
        self.XB = np.interp(self.U, self.CX[self.I:][::-1], self.XG[self.I:][::-1],
                            left=self.XG[self.I], right=self.XG[-1])
        floor = 1e-3 * np.median(np.abs(self.DCDX[np.abs(self.DCDX) > 1e-6]))
        self.DA = np.maximum(np.abs(np.interp(self.XA, self.XG, self.DCDX)), floor)
        self.DB = np.maximum(np.abs(np.interp(self.XB, self.XG, self.DCDX)), floor)
        self.DU = np.gradient(self.U)
        self.KPG = np.logspace(0.0, np.log10(4e5), self.nk)
        self._jl = None

    def build_bessel(self, nproc=1):
        """Cache j_l(k' u) on the (ell, k', u) grid; ~1 GB and ~10 min."""
        def rows(chunk):
            return np.stack([spherical_jn(int(e), self.KPG[:, None] * self.U[None, :])
                             for e in chunk])
        if nproc > 1:
            import multiprocessing as mp
            chunks = [c for c in np.array_split(np.asarray(self.ells), nproc) if len(c)]
            with mp.get_context("fork").Pool(nproc) as pool:
                self._jl = np.concatenate(pool.map(rows, chunks), axis=0).astype(np.float32)
        else:
            self._jl = rows(np.asarray(self.ells)).astype(np.float32)

    def G_of_u(self, p, with_metric=False, a0_met=None):
        """Two-branch kernel G(u), Eq. (158).  Optionally adds W_met."""
        W = emission_kernel(self.XG, self.UG, self.x_star, p["p"], p["q"], p["tau2"])
        if with_metric:
            a0m = p["A0_low"] if a0_met is None else a0_met
            W = W + metric_kernel(self.XG, self.UG, a0m, p.get("beta", 5.81))
        G = np.interp(self.XA, self.XG, W) / self.DA + np.interp(self.XB, self.XG, W) / self.DB
        s = np.sum(G * self.DU)
        return None if (not np.isfinite(s) or s <= 0) else G / s

    def Dl_model(self, p, cutoff=None, with_metric=False, a0_met=None):
        """D_ell = l(l+1)C_l/2pi of the model.  [muK^2]"""
        if self._jl is None:
            self.build_bessel()
        G = self.G_of_u(p, with_metric=with_metric, a0_met=a0_met)
        k_phys = self.KPG * (p["H0"] / C_LIGHT)
        Pk = power_spectrum(k_phys, p, cutoff=cutoff)
        pref = (p["H0"] / C_LIGHT) ** 3 * self.KPG ** 2 * Pk
        win = np.tensordot(self._jl, (G * self.DU).astype(np.float32), axes=([2], [0]))
        integ = pref[None, :] * win.astype(np.float64) ** 2
        Cl = (2.0 / np.pi) * np.sum(integ * np.gradient(self.KPG)[None, :], axis=1)
        return self.ells * (self.ells + 1) / (2.0 * np.pi) * Cl

    def chi2(self, p, cutoff=None, with_metric=False, a0_met=None):
        model = self.Dl_model(p, cutoff=cutoff, with_metric=with_metric, a0_met=a0_met)
        return float(np.sum(((self.Dl - model) / self.err) ** 2)), len(self.ells)
