#!/usr/bin/env python3
"""
Semi-analytical uniaxial solution of the 4-fiber-family kinematic-growth model
(MBE-KG), matching the custom FEBio material FEKinematicGrowthZ.cpp on the
single-element uniaxial coupon (aneurysm-knockout/uniax-fractions/uniax-*.feb).

Model (per FEKinematicGrowthZ.cpp, hyperels_mod=3), extra 2nd-PK stress:
    Se = phi_e * mu * Ge^2                              (elastin, neo-Hookean)
       + phi_m * Sh_m * (N1 x N1)                       (SMC, circumferential)
       + phi_c * [ Sh_ct*beta_t*(N1 x N1)               (collagen, 4 families)
                 + Sh_cz*beta_z*(N2 x N2)
                 + Sh_cp*beta_d*(Np x Np)
                 + Sh_cn*beta_d*(Nn x Nn) ]
    Sh_i = c_i (l_i^2 - 1) exp(d_i (l_i^2-1)^2) G_i^2 , l_i^2 from Fe (elastic).
Directions: N1=circumferential, N2=axial (z), Np/Nn=diagonal at +/-alpha from
axial in the circ-axial plane.

Coupon geometry (AXIAL coupon, 2026-07-14):
the .feb builds the fibre CS from the Gauss point's position about the global z axis
(Nt = (-y,x,0)/r; Nz = z from the identity local CS). The cube is placed away from that
axis, so the basis is UNIFORM over the element: N_r = x, N_theta = y, N_z = z. The
sample is pulled along z, i.e. along the AXIAL direction, which is the canonical coupon
(the fibre angle alpha_0 is itself measured from the axial direction).

Axes here are ordered (pull, free, free) = (axial, circumferential, radial), so the
referential structural tensors are exactly diagonal and need no azimuthal averaging:
    <N2 x N2> = diag(1, 0, 0)                 (axial, the pull direction)
    <N1 x N1> = diag(0, 1, 0)                 (circumferential)
    <Np x Np> = <Nn x Nn> = diag(cos^2a, sin^2a, 0)   (diagonal, in the axial-circ plane)
    no fibre family lies along the radial direction, which is carried by elastin alone.
This SUPERSEDES the earlier azimuthally-averaged tensors, which were an artefact of the
cube straddling the z axis (the fibre directions then spread 180 deg within one element).

Loading: axial force along x (global), lateral faces y,z traction-free
(sigma_yy = sigma_zz = 0).  Passive: theta=1, near-incompressible penalty.
MBE: isotropic growth Fg = J^(1/3) I, Lagrange multiplier p = sigma_v^x - sigma_vo
holds the homeostatic volumetric Cauchy stress sigma_vo.

Validation vs FEBio (uniax-MBEdense/PASSdense, 8 F_o run):
  passive sigma_xx agrees to < 0.1% and lam_y/lam_z to 4 decimals over the whole
  branch; MBE matches once sigma_vo is taken as FEBio stores it. Residual sub-1%
  differences trace to the near-incompressibility PENALTY used in the FE simulation
  (k_J=1e6, Jdep=0.9999), which admits a small compressibility (J~1.0001) that a
  strictly self-consistent analytical sigma_vo does not.

Validation on the AXIAL coupon (2026-07-14): passive 0.08%, MBE 0.00% over the whole
F_o -> 8 F_o branch, with NO calibration factor. See verify_coupon.py.
"""
import numpy as np
from scipy.optimize import fsolve, brentq

# ------------------------------------------------------------------ parameters
# knockout set, unit deposition stretches (uniax-MBE.feb)
c_e = 89.71
c_m = 261.4; d_m = 0.240
c_c = 234.9; d_c = 4.08
alpha = np.deg2rad(29.91)
beta_t = 0.056; beta_z = 0.067
beta_d = 0.5*(1.0 - beta_t - beta_z)
f_e = 0.34; f_c = 0.33
phi_e, phi_c, phi_m = f_e, f_c, 1.0 - f_e - f_c       # 0.34, 0.33, 0.33
Gm = Gc = Get = Gez = 1.0
kJ = 1.0e6; Jdep = 0.9999
F_o_total = 4*48.82                                    # homeostatic axial force (mN)

sa2, ca2 = np.sin(alpha)**2, np.cos(alpha)**2
# element-averaged referential structural tensors (diagonal, xyz)
# axes ordered (pull, free, free) = (axial, circumferential, radial)
T_ax   = np.array([1.0, 0.0, 0.0])                     # N2 (axial) = the pull direction
T_circ = np.array([0.0, 1.0, 0.0])                     # N1 (circumferential)
T_diag = np.array([ca2, sa2, 0.0])                     # Np, Nn (diagonal, axial-circ plane)
Ge2    = np.array([Gez**2, Get**2, 1.0/(Get*Gez)**2])  # Ge^2 (= I for unit stretches)

def Se_diag(lx, ly, lz, Jinv13=1.0):
    """Extra 2nd-PK stress (diagonal xyz).  Jinv13 = J^(-1/3): elastic stretch
    lei = li * J^(-1/3) (Fe = J^(-1/3) F).  Passive -> Jinv13 = 1."""
    lex2, ley2, lez2 = (lx*Jinv13)**2, (ly*Jinv13)**2, (lz*Jinv13)**2
    # uniform basis, axes (pull, free, free) = (axial, circumferential, radial):
    # index 0 = axial, 1 = circumferential, 2 = radial. No fibre lies along r.
    lto2 = ley2                               # circumferential fibre stretch^2
    lzo2 = lex2                               # axial (= the pull direction)
    lpo2 = ca2*lex2 + sa2*ley2                # diagonal, +-alpha from axial
    def Sh(c, d, l2, G): return c*(l2-1.0)*np.exp(d*(l2-1.0)**2)*G*G
    Sh_m  = Sh(c_m, d_m, Gm*Gm*lto2, Gm)
    Sh_ct = Sh(c_c, d_c, Gc*Gc*lto2, Gc)
    Sh_cz = Sh(c_c, d_c, Gc*Gc*lzo2, Gc)
    Sh_cp = Sh(c_c, d_c, Gc*Gc*lpo2, Gc)      # Sh_cn == Sh_cp (symmetric families)
    Se  = phi_e*c_e*Ge2.copy()
    Se += phi_m*Sh_m*T_circ
    Se += phi_c*(Sh_ct*beta_t*T_circ + Sh_cz*beta_z*T_ax + 2.0*Sh_cp*beta_d*T_diag)
    return Se

# ------------------------------------------------------------------- passive
def passive_stress(lx, ly, lz):
    J = lx*ly*lz
    Se = Se_diag(lx, ly, lz, 1.0)
    pen = kJ*np.log(Jdep*J)
    return (np.array([lx, ly, lz])**2 * Se + pen)/J   # sigma_xx, yy, zz

def passive_solve(lx, guess=None):
    g = guess if guess is not None else [lx**-0.5, lx**-0.5]
    ly, lz = fsolve(lambda v: passive_stress(lx, v[0], v[1])[1:], g, full_output=False)
    s = passive_stress(lx, ly, lz)
    J = lx*ly*lz
    return dict(lx=lx, ly=ly, lz=lz, J=J, sxx=s[0], syy=s[1], szz=s[2],
                force=s[0]*ly*lz)                       # Cauchy*current area = nominal force

def homeostatic_svo():
    """Passive state at the homeostatic force F_o -> sigma_vo = (1/3) tr(sigma)."""
    lx_o = brentq(lambda lx: passive_solve(lx)['force'] - F_o_total, 1.0, 3.0)
    st = passive_solve(lx_o)
    svo = (st['sxx'] + st['syy'] + st['szz'])/3.0
    return lx_o, svo, st

# ---------------------------------------------------------------------- MBE
def mbe_solve(lx, svo, guess=None):
    """Given axial stretch lx and homeostatic sigma_vo, solve lateral stretches."""
    def resid(v):
        ly, lz = v
        J = lx*ly*lz; Ji13 = J**(1.0/3.0)
        Se = Se_diag(lx, ly, lz, J**(-1.0/3.0))
        li2 = np.array([lx, ly, lz])**2
        # sigma_yy = sigma_zz = 0  <=>  ly^2 Se_yy = lz^2 Se_zz  (eq1)
        # p = sigma_v^x - svo, with p = (1/J) ly^2 Ji13 Se_yy    (eq2)
        eq1 = li2[1]*Se[1] - li2[2]*Se[2]
        svx = (Ji13/(3.0*J))*np.sum(li2*Se)
        p_lat = (Ji13/J)*li2[1]*Se[1]
        eq2 = p_lat - (svx - svo)
        return [eq1, eq2]
    g = guess if guess is not None else [lx**-0.4, lx**-0.4]
    ly, lz = fsolve(resid, g)
    J = lx*ly*lz; Ji13 = J**(1.0/3.0)
    Se = Se_diag(lx, ly, lz, J**(-1.0/3.0))
    li2 = np.array([lx, ly, lz])**2
    p = (Ji13/J)*li2[1]*Se[1]
    sxx = (Ji13/J)*li2[0]*Se[0] - p
    return dict(lx=lx, ly=ly, lz=lz, J=J, sxx=sxx, p=p)

# --------------------------------------------------- G&R time evolution (F(s))
# Growth law (main.tex Eq. stimulus-fc + p-vol-stress), backward Euler:
#     p^{n+1} = sigma_v^x - sigma_vo * (1 + (1/K_sig) * ln(J^{n+1}/J^n)/(k*ds))
# Kinematics enforce J = vartheta (Box 1), so growth is NOT a separate state
# variable: each step is a nonlinear ALGEBRAIC solve in (lx, ly, lz), with the
# history entering only through J^n. The MBE model is the same system with the
# ln(J/J^n) term dropped (rate-independent limit, ds -> inf).
K_sig = 1.0        # gain (no <K_sigma> tag in the .feb -> 1, as in Table 2)
k_gr  = 0.1        # 1/day  (<k_parameter>)

def _sig_grown(lx, ly, lz, svo, Jn=None, ds=None):
    """Total Cauchy stress (diagonal) in the grown state. Jn/ds given -> full KG
    (rate-dependent p); Jn=None -> MBE (equilibrated p = sigma_v^x - sigma_vo)."""
    J = lx*ly*lz
    Se = Se_diag(lx, ly, lz, J**(-1.0/3.0))
    li2 = np.array([lx, ly, lz])**2
    sx = (J**(1.0/3.0)/J) * li2 * Se            # extra Cauchy stress sigma^x
    svx = sx.sum()/3.0
    if Jn is None:
        p = svx - svo                                              # MBE
    else:
        p = svx - svo*(1.0 + np.log(J/Jn)/(K_sig*k_gr*ds))         # full KG
    return sx - p, J

def force_solve(F, svo, Jn=None, ds=None, guess=None):
    """Uniaxial force control: sigma_yy = sigma_zz = 0 and sigma_xx*(ly*lz*h0^2) = F.
    Three equations, three unknowns (lx, ly, lz)."""
    def resid(v):
        s, _ = _sig_grown(v[0], v[1], v[2], svo, Jn, ds)
        return [s[0]*v[1]*v[2] - F, s[1], s[2]]
    g = guess if guess is not None else [1.3, 0.9, 0.9]
    x = fsolve(resid, g)
    s, J = _sig_grown(x[0], x[1], x[2], svo, Jn, ds)
    return dict(lx=x[0], ly=x[1], lz=x[2], J=J, sxx=s[0], lg=J**(1.0/3.0))

def evolve(s_load, ds, svo, lx_o, mode="full", s_max=5.0):
    """March the G&R response for a force ramp F_o -> 3 F_o over s_load, then held.
    Returns columns matching the FEBio plot files *_c.txt:
        s/s_load , sigma_xx [kPa] , lambda/lambda_o , lambda_g
    (s is G&R time, i.e. FEBio time t minus gr_start = 1 day.)"""
    n = int(round(s_max*s_load/ds))
    # seed from the homeostatic passive state (s=0, F=F_o, J=1): without it the
    # first solve can land on a spurious root of the 3-equation system.
    st_o = passive_solve(lx_o)
    out, Jn, guess = [], 1.0, [lx_o, st_o['ly'], st_o['lz']]
    for i in range(n + 1):
        s = i*ds
        F = F_o_total*(1.0 + 2.0*min(s/s_load, 1.0))     # ramp to 3 F_o, then hold
        st = force_solve(F, svo, Jn if mode == "full" else None, ds, guess)
        guess = [st['lx'], st['ly'], st['lz']]
        if mode == "full":
            Jn = st['J']                                  # history update
        out.append((s/s_load, st['sxx'], st['lx']/lx_o, st['lg']))
    return np.array(out)

# ------------------------------------------------------------------ validation
def load_febio(path):
    import os
    if not os.path.exists(path): return None
    d = np.loadtxt(path)                                # lam1 sx J lam2 lam3
    return d

import warnings; warnings.filterwarnings("ignore")

def march(fb, solver, svo=None, lx_o=None):
    """Solve along the FEBio lam_x samples with continuation (previous
    lateral stretches seed the next solve), returning aligned rows.
    For MBE below the homeostatic stretch lx_o, growth is inactive and the
    response is passive (that is the regime FEBio is in before F_o)."""
    out, guess = [], None
    for row in fb:
        lx = row[0]
        if lx < 1.0001:
            guess = [1.0, 1.0]; continue
        if solver == "passive" or (lx_o is not None and lx <= lx_o):
            s = passive_solve(lx, guess)
        else:
            s = mbe_solve(lx, svo, guess)
        guess = [s['ly'], s['lz']]
        out.append((lx, s, row))
    return out

if __name__ == "__main__":
    import os
    base = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "aneurysm-knockout", "uniax-fractions"))
    lx_o, svo, st_o = homeostatic_svo()
    print(f"Homeostatic (F_o={F_o_total:.2f} mN): lam_x={lx_o:.4f}  "
          f"sigma_xx={st_o['sxx']:.2f}  J={st_o['J']:.5f}  sigma_vo={svo:.2f} kPa")
    print("(FEBio, AXIAL coupon: sigma_o=260.59 kPa -> agreement is now exact)\n")

    # The point-by-point F_o -> 8 F_o comparison lives in verify_coupon.py, which reads
    # the FEBio LOGS of the coupon decks (node7_*.txt / elem1_*.txt). It reports 0.08%
    # (passive) and 0.00% (MBE) over the whole branch.
    #
    # The old block here read figM_*.dat and applied an "MBE-cal" correction, solving at
    # 0.995 F_o because FEBio stored sigma_vo at the last sub-step before growth onset.
    # BOTH ARE GONE: the sub-step bug was fixed (must-points), so no calibration is needed,
    # and the .dat files are no longer written. It was silently comparing against
    # swirl-era data and reporting a 460% error.

    # ---------------------------------------------------------------- G&R runs
    # Reproduce the four curves of the paper's uniaxial figure and compare with
    # the FEBio plot data (*_c.txt from prep_uniax_consistent.py).
    print("=" * 78)
    print("G&R evolution: semi-analytical vs FEBio  (figure tension-4fib-all-again)")
    print("=" * 78)
    print(f"lambda_o = {lx_o:.4f}   sigma_vo = {svo:.2f} kPa   "
          f"(MBE plateau = 3*sigma_vo = {3*svo:.1f} kPa)\n")

    # ds = the FEBio step size of each run (backward Euler is ds-dependent)
    cases = [("MBE",   None, 10.0, 1.0,  "mbe"),
             ("sp1",   1.0,   1.0, 0.01, "full"),
             ("sp10",  10.0, 10.0, 0.10, "full"),
             ("sp100", 100.0, 100.0, 0.10, "full")]
    for tag, s_load, s_load_plot, ds, mode in cases:
        sl = s_load_plot
        d = evolve(sl, ds, svo, lx_o, mode=mode)
        np.savetxt(os.path.join(base, f"{tag}_a.txt"), d, fmt="%.6g",
                   header="s/s_load  sigma_xx  lam/lam_o  lambda_g  (semi-analytical)")
        fb = os.path.join(base, f"{tag}_c.txt")
        msg = ""
        if os.path.exists(fb):
            f = np.loadtxt(fb)
            f = f[f[:, 0] >= 0]                      # G&R phase only
            e = []
            for col in (1, 2, 3):
                a = np.interp(f[:, 0], d[:, 0], d[:, col])
                ok = np.abs(f[:, col]) > 1e-6
                e.append(np.max(np.abs(a[ok] - f[ok, col]) / np.abs(f[ok, col])) * 100)
            msg = (f"| vs FEBio, worst rel. err:  sigma_xx {e[0]:5.2f}%   "
                   f"lam/lam_o {e[1]:5.2f}%   lambda_g {e[2]:5.2f}%")
        print(f"  {tag:5s} ({mode:4s}, ds={ds:g}) end: sxx={d[-1,1]:7.1f} "
              f"lam/lamo={d[-1,2]:.3f} lg={d[-1,3]:.3f}  {msg}")
    print("\n  -> wrote <case>_a.txt (semi-analytical) next to <case>_c.txt (FEBio)")
