"""
tests for tb_gauge = true: conversion of real-space operators from the wannier90 convention to the tb gauge
(util/get_oper_R.py:to_tb_gauge) and the k-space transforms with orbital positions (util/utility.py).

Reference model: an orthonormal frame V(k) = sum_R e^{ik.R} V_R (D x W, V^† V = 1; the first W columns of
U_0 E_1(k) U_1 E_2(k) ... with unitary U_j and E_j(k) = 1 - P_j + P_j e^{ik.n_j}, P_j a projector) and a Hamiltonian
H(k) = V h(k) V^† + Q H_out Q (Q = 1 - V V^†) for which span V(k) is invariant, as for isolated bands.
The operators in the wannier90 convention are, e.g., A_a(k) = i V^† d_a V, BB_a(k) = i V^† H d_a V,
CC_ab(k) = (d_a V)^† H d_b V; all are trigonometric polynomials, so their real-space matrices are obtained exactly by a
discrete Fourier transform on a large enough mesh. In the tb gauge the frame is V^I(k) = V(k) e^{ik.τ}; its operators are
evaluated directly, with derivatives by central finite differences.
"""

import numpy as np
import pytest

from symclosestwannier.analyzer.get_response import (
    berry_get_imfgh_klist,
    berry_get_js_k,
    wham_get_D_h,
    wham_get_deleig,
)
from symclosestwannier.util.constants import hbar_SI
from symclosestwannier.util.get_oper_R import get_HH_R, get_v_R, to_tb_gauge
from symclosestwannier.util.utility import (
    fourier_transform_r_to_k,
    fourier_transform_r_to_k_new,
    fourier_transform_r_to_k_vec,
    tb_gauge_positions,
    wigner_seitz,
)

A = np.array([[2.5, 0.0, 0.0], [-1.25, 2.165, 0.0], [0.3, 0.4, 7.0]])
B = 2 * np.pi * np.linalg.inv(A).T
# orbital positions (fractional), one outside [0, 1).
TAU = np.array([[0.1, 0.2, 0.05], [0.6, 0.3, -0.2], [1.3, 0.7, 0.4]])
NS = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
DIM, NUM_WANN = 5, 3
# harmonics of the operators lie in [-3, 3]^3, so a 7^3 mesh with R in [-3, 3]^3 represents them exactly.
NMESH = 7
H_STEP = 1e-5


# ==================================================
class Model:
    def __init__(self, seed=0):
        rng = np.random.default_rng(seed)

        def herm(n):
            M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            return M + M.conj().T

        def unitary(n):
            return np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]

        def projector(n, rank):
            Q = unitary(n)[:, :rank]
            return Q @ Q.conj().T

        def mul(X, Y):
            Z = {}
            for r1, x in X.items():
                for r2, y in Y.items():
                    r = tuple(np.add(r1, r2))
                    Z[r] = Z.get(r, 0) + x @ y
            return Z

        # U(k) as {R: U_R}
        U = {(0, 0, 0): unitary(DIM)}
        for n in NS:
            P = projector(DIM, 2)
            U = mul(U, {(0, 0, 0): np.eye(DIM) - P, n: P})
            U = mul(U, {(0, 0, 0): unitary(DIM)})

        self.NV = np.array(list(U.keys()))
        self.V_R = np.array([U[n][:, :NUM_WANN] for n in U])
        self.nv_cart = self.NV @ A

        # h(k) = h_0 + sum_j (h_j e^{ik.n_j} + h.c.)
        h1 = rng.normal(size=(NUM_WANN, NUM_WANN)) + 1j * rng.normal(size=(NUM_WANN, NUM_WANN))
        self.h_R = {(0, 0, 0): herm(NUM_WANN), (1, 0, 0): h1, (-1, 0, 0): h1.conj().T}
        self.H_out = herm(DIM)
        self.sigma = np.array([herm(DIM) for _ in range(3)])

        irvec = [(i, j, l) for i in range(-3, 4) for j in range(-3, 4) for l in range(-3, 4)]
        self.irvec = np.array(irvec)
        self.ndegen = np.ones(len(irvec), dtype=int)

    # k: fractional coordinates.
    def V(self, k):
        return np.einsum("r,rdw->dw", np.exp(2j * np.pi * self.NV @ k), self.V_R)

    def dV(self, k):
        """analytic cartesian derivative, (3, D, W)."""
        return np.einsum("r,ra,rdw->adw", np.exp(2j * np.pi * self.NV @ k), 1j * self.nv_cart, self.V_R)

    def VI(self, k):
        return self.V(k) * np.exp(2j * np.pi * TAU @ k)[None, :]

    def dVI(self, k):
        """finite-difference cartesian derivative of the tb-gauge frame, (3, D, W)."""
        dk = H_STEP * np.linalg.inv(B)  # rows: fractional steps along cartesian x, y, z
        return np.array([(self.VI(k + dk[a]) - self.VI(k - dk[a])) / (2 * H_STEP) for a in range(3)])

    def H(self, k):
        """Hamiltonian in the full space, span V(k) is invariant."""
        V = self.V(k)
        h = sum(np.exp(2j * np.pi * np.dot(k, R)) * hR for R, hR in self.h_R.items())
        Q = np.eye(DIM) - V @ V.conj().T
        return V @ h @ V.conj().T + Q @ self.H_out @ Q

    # ==================================================
    def operators_k(self, k, V, dV):
        H, S = self.H(k), self.sigma
        Vd = V.conj().T
        return {
            "HH": Vd @ H @ V,
            "AA": np.array([1j * Vd @ dV[a] for a in range(3)]),
            "BB": np.array([1j * Vd @ H @ dV[a] for a in range(3)]),
            "CC": np.array([[dV[a].conj().T @ H @ dV[b] for b in range(3)] for a in range(3)]),
            "SS": np.array([Vd @ S[a] @ V for a in range(3)]),
            "SR": np.array([[1j * Vd @ S[a] @ dV[b] for b in range(3)] for a in range(3)]),
            "SH": np.array([Vd @ S[a] @ H @ V for a in range(3)]),
            "SHR": np.array([[1j * Vd @ S[a] @ H @ dV[b] for b in range(3)] for a in range(3)]),
        }

    def operators_R(self):
        """
        wannier90-convention operators, O(R) = 1/N_k sum_k e^{-ik.R} O(k), not divided by ndegen (= 1).
        """
        if not hasattr(self, "_operators_R"):
            kpoints = np.array([[i, j, l] for i in range(NMESH) for j in range(NMESH) for l in range(NMESH)]) / NMESH
            ops = [self.operators_k(k, self.V(k), self.dV(k)) for k in kpoints]
            phase = np.exp(-2j * np.pi * kpoints @ self.irvec.T) / len(kpoints)
            d = {}
            for key in ops[0]:
                Ok = np.array([o[key] for o in ops])  # (k, ..., m, n)
                Ok = np.moveaxis(Ok, 0, -3)  # (..., k, m, n)
                d[f"{key}_R"] = np.einsum("kR,...kmn->...Rmn", phase, Ok)
            self._operators_R = d

        return {k: v.copy() for k, v in self._operators_R.items()}


KPOINTS = np.array([[0.13, 0.27, 0.41], [-0.31, 0.05, 0.22], [0.5, -0.17, 0.08]])


def to_k(O_R, model, k, atoms_frac):
    """Fourier transform of an operator with any number of leading cartesian indices."""
    shape = O_R.shape[:-3]
    O = O_R.reshape((-1,) + O_R.shape[-3:])
    Ok = np.array([fourier_transform_r_to_k(o, np.array([k]), model.irvec, model.ndegen, atoms_frac)[0] for o in O])
    return Ok.reshape(shape + Ok.shape[-2:])


# ==================================================
def test_reference_model_is_consistent():
    """the analytic real-space operators reproduce the frame expressions (wannier90 convention, no orbital phase)."""
    model = Model()
    ops_R = model.operators_R()
    for k in KPOINTS:
        np.testing.assert_allclose(model.V(k).conj().T @ model.V(k), np.eye(NUM_WANN), rtol=0, atol=1e-12)
        ref = model.operators_k(k, model.V(k), model.dV(k))
        for key, Ok in ref.items():
            np.testing.assert_allclose(to_k(ops_R[f"{key}_R"], model, k, None), Ok, rtol=0, atol=1e-10)


# ==================================================
def test_to_tb_gauge_matches_tb_gauge_frame():
    """
    to_tb_gauge followed by the transform with orbital positions gives the operators of the frame V(k) e^{ik.τ},
    at k points that are not on any mesh.
    """
    model = Model()
    ops_I = to_tb_gauge(model.operators_R(), model.irvec, A, TAU)

    for k in KPOINTS:
        ref = model.operators_k(k, model.VI(k), model.dVI(k))
        for key, Ok in ref.items():
            np.testing.assert_allclose(to_k(ops_I[f"{key}_R"], model, k, TAU), Ok, rtol=0, atol=1e-6, err_msg=key)

    # BB, CC, SR and SHR do change.
    ops_II = model.operators_R()
    for key in ("AA_R", "BB_R", "CC_R", "SR_R", "SHR_R"):
        assert np.max(np.abs(ops_I[key] - ops_II[key])) > 0.1, key
    for key in ("HH_R", "SS_R", "SH_R"):
        np.testing.assert_array_equal(ops_I[key], ops_II[key])


# ==================================================
def test_to_tb_gauge_does_not_modify_input_and_checks_requirements():
    model = Model()
    ops = model.operators_R()
    copy = {k: v.copy() for k, v in ops.items()}
    to_tb_gauge(ops, model.irvec, A, TAU)
    for key in ops:
        np.testing.assert_array_equal(ops[key], copy[key])

    with pytest.raises(ValueError, match="HH_R"):
        to_tb_gauge({"BB_R": ops["BB_R"]}, model.irvec, A, TAU)
    with pytest.raises(ValueError, match="BB_R"):
        to_tb_gauge({"HH_R": ops["HH_R"], "CC_R": ops["CC_R"]}, model.irvec, A, TAU)
    with pytest.raises(ValueError, match="SS_R"):
        to_tb_gauge({"SR_R": ops["SR_R"]}, model.irvec, A, TAU)
    with pytest.raises(ValueError, match="SH_R"):
        to_tb_gauge({"SHR_R": ops["SHR_R"]}, model.irvec, A, TAU)


# ==================================================
def test_derivatives_with_orbital_positions():
    """
    k derivatives taken with orbital positions in the phase: dH/dk (fourier_transform_r_to_k_new), the velocity of
    get_v_R, and the curl of A (fourier_transform_r_to_k_vec, pseudo=True).
    """
    model = Model()
    ops_I = to_tb_gauge(model.operators_R(), model.irvec, A, TAU)
    cwi = {"unit_cell_cart": A, "irvec": model.irvec}
    v_R = get_v_R(cwi, ops_I["HH_R"], TAU)
    dk = H_STEP * np.linalg.inv(B)

    def H_I(k):
        return model.operators_k(k, model.VI(k), model.dVI(k))["HH"]

    def A_I(k):
        return model.operators_k(k, model.VI(k), model.dVI(k))["AA"]

    for k in KPOINTS:
        _, delHH = fourier_transform_r_to_k_new(ops_I["HH_R"], np.array([k]), A, model.irvec, model.ndegen, TAU)
        dH = np.array([(H_I(k + dk[a]) - H_I(k - dk[a])) / (2 * H_STEP) for a in range(3)])
        np.testing.assert_allclose(delHH[:, 0], dH, rtol=0, atol=1e-5)
        np.testing.assert_allclose(to_k(v_R, model, k, TAU) * hbar_SI, dH, rtol=0, atol=1e-5)

        _, Omega = fourier_transform_r_to_k_vec(ops_I["AA_R"], np.array([k]), model.irvec, model.ndegen, TAU, A, pseudo=True)
        dA = np.array([(A_I(k + dk[a]) - A_I(k - dk[a])) / (2 * H_STEP) for a in range(3)])
        curl = np.array([dA[1][2] - dA[2][1], dA[2][0] - dA[0][2], dA[0][1] - dA[1][0]])
        np.testing.assert_allclose(Omega[:, 0], curl, rtol=0, atol=1e-4)


# ==================================================
def make_cwi(model, tb_gauge):
    return {
        "tb_gauge": tb_gauge,
        "atoms_frac": {("X", i + 1): list(t) for i, t in enumerate(TAU)},
        "nw2n": list(range(NUM_WANN)),
        "num_wann": NUM_WANN,
        "num_fermi": 1,
        "fermi_energy_list": [0.0],
        "unit_cell_cart": A,
        "irvec": model.irvec,
        "ndegen": model.ndegen,
        "zeeman_interaction": False,
        "shc_alpha": 1,
        "shc_gamma": 3,
    }


def operators_for(model, tb_gauge):
    ops = model.operators_R()
    if tb_gauge:
        ops = to_tb_gauge(ops, model.irvec, A, TAU)
    return ops


# ==================================================
def test_orbital_magnetization_terms_do_not_depend_on_tb_gauge():
    """
    -2Im f, -2Im g, -2Im h (Berry curvature and orbital magnetization, with the LLambda term from CC) at each k
    are the same in both gauges. Only the sums of the J0, J1 and J2 terms are gauge invariant.
    """
    model = Model()
    occ = np.array([[1.0, 0.0, 0.0]] * len(KPOINTS))

    res = {}
    for tb_gauge in (False, True):
        cwi = make_cwi(model, tb_gauge)
        res[tb_gauge] = berry_get_imfgh_klist(cwi, operators_for(model, tb_gauge), KPOINTS, imf=True, img=True, imh=True, occ=occ)

    for f_II, f_I in zip(res[False], res[True]):
        f_II, f_I = f_II.sum(axis=2), f_I.sum(axis=2)
        assert np.max(np.abs(f_II)) > 0.1
        np.testing.assert_allclose(f_I, f_II, rtol=0, atol=1e-6)


# ==================================================
def test_spin_current_does_not_depend_on_tb_gauge():
    """
    matrix elements of the spin current (QZYZ18 Eq. (23)) between eigenstates agree up to the phases of the eigenvectors.
    """
    model = Model()
    alpha = 0

    res = {}
    for tb_gauge in (False, True):
        cwi = make_cwi(model, tb_gauge)
        ops = operators_for(model, tb_gauge)
        atoms_frac = tb_gauge_positions(cwi)
        HH, delHH = fourier_transform_r_to_k_new(ops["HH_R"], KPOINTS, A, model.irvec, model.ndegen, atoms_frac)
        E, U = np.linalg.eigh(HH)
        delE = wham_get_deleig(delHH, E, U)
        D_h = wham_get_D_h(delHH, E, U)
        res[tb_gauge] = berry_get_js_k(cwi, ops, KPOINTS, E, delE[alpha], D_h[alpha], U)

    assert np.max(np.abs(res[False])) > 0.1
    np.testing.assert_allclose(np.abs(res[True]), np.abs(res[False]), rtol=1e-7, atol=1e-7)


# ==================================================
def test_get_HH_R_does_not_depend_on_tb_gauge():
    """
    get_oper_R gives operators in the wannier90 convention: with tb_gauge = true the Hamiltonian from the DFT mesh is
    the same and interpolates correctly between mesh points.
    """
    model = Model()
    N = 6
    kpoints = np.array([[i / N, j / N, l / N] for i in range(N) for j in range(N) for l in range(N)])
    Hk = np.array([model.operators_k(k, model.V(k), model.dV(k))["HH"] for k in kpoints])
    Ek, Uk = np.linalg.eigh(Hk)

    irvec, ndegen = wigner_seitz(A, [N, N, N])
    for tb_gauge in (False, True):
        cwi = {**make_cwi(model, tb_gauge), "Ek": Ek, "Uk": Uk.transpose(0, 2, 1).conj(), "kpoints": kpoints, "irvec": irvec}
        HH_R = get_HH_R(cwi)
        for k in KPOINTS:
            Hk_ref = model.operators_k(k, model.V(k), model.dV(k))["HH"]
            Hk_int = fourier_transform_r_to_k(HH_R, np.array([k]), irvec, ndegen)[0]
            np.testing.assert_allclose(Hk_int, Hk_ref, rtol=0, atol=1e-8)


# ==================================================
@pytest.fixture
def restore_state(monkeypatch):
    """
    CWManager changes the current directory and appends it to sys.path.
    """
    import sys

    monkeypatch.chdir(".")
    monkeypatch.setattr(sys, "path", list(sys.path))


@pytest.mark.parametrize("tb_gauge", [False, True])
@pytest.mark.parametrize("berry_task", ["morb", "shc"])
def test_response_converts_operators_once(tmp_path, monkeypatch, restore_state, tb_gauge, berry_task):
    """
    Response.set_operators takes the operators of get_oper_R (wannier90 convention) and converts them once for tb_gauge.
    """
    import symclosestwannier.analyzer.response as response
    from symclosestwannier.cw.cw_manager import CWManager

    model = Model()
    ops = model.operators_R()

    def fake_get_oper_R(name, cwi, *args):
        if name == "SHC_R":
            return ops["SR_R"].copy(), ops["SHR_R"].copy(), ops["SH_R"].copy()
        return ops[name].copy()

    monkeypatch.setattr(response, "get_oper_R", fake_get_oper_R)

    cwi = {
        **make_cwi(model, tb_gauge),
        "seedname": "model",
        "use_tb_approximation": False,
        "berry": True,
        "berry_task": berry_task,
        "gyrotropic": False,
        "spin_moment": False,
        "spin_decomp": False,
    }
    cwm = CWManager(topdir=str(tmp_path), verbose=False, parallel=False, formatter=False)
    res = response.Response(cwi, cwm)
    res.set_operators()

    expected = to_tb_gauge(ops, model.irvec, A, TAU) if tb_gauge else ops
    keys = ("HH_R", "AA_R", "BB_R", "CC_R") if berry_task == "morb" else ("HH_R", "AA_R", "SS_R", "SR_R", "SHR_R", "SH_R")
    for key in keys:
        np.testing.assert_allclose(res[key], expected[key], rtol=0, atol=1e-12, err_msg=key)
