"""s1 copy (3 Oct 2026): reads the nu=* folders of s1/param-x (the 12 profiles of run_maps.py).
Pandas-free equivalent of read_all_fused.py (this python3 has numpy but not pandas).
Rebuilds MASTER_PARAMETRIC_STUDY.csv from the nu=* folders, identical logic/columns."""
import os, csv, numpy as np

PX = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(PX, "MASTER_PARAMETRIC_STUDY.csv")
REF = {'axis': (-2.31774149752e-15, 0.0, 0.605691058705, 0.0),
       'asym': ( 0.640459751827,   0.0, 0.101438859344, 0.0)}

def radial_stretch(x, u, y, v):
    r_cur = np.sqrt(x**2 + y**2); r_ref = np.sqrt((x-u)**2 + (y-v)**2)
    return np.where(r_ref > 1e-6, r_cur / r_ref, 1.0)

def main():
    rows = []
    carpetas = [d for d in os.listdir(PX)
                if os.path.isdir(os.path.join(PX, d)) and d.startswith('nu=')]
    for carpeta in carpetas:
        p = carpeta.split('_')
        try:
            nu_val = float(p[0].split('=')[1]); alph_val = p[1].split('=')[1]
            zod_str = p[2].split('=')[1]
            if 'unit-frac' in zod_str: zod_val, nota = 1.0, 'unit-frac'
            else:                       zod_val, nota = float(zod_str), 'normal'
        except Exception as e:
            print(f"skip malformed {carpeta}: {e}"); continue
        for type_ in ('axis', 'asym'):
            fEA = os.path.join(PX, carpeta, f"summary_{type_}_elem_A.csv")
            fNA = os.path.join(PX, carpeta, f"summary_{type_}_node_A.csv")
            if not (os.path.exists(fEA) and os.path.exists(fNA)): continue
            try:
                dEA = np.loadtxt(fEA, delimiter=',', comments='#')
                dNA = np.loadtxt(fNA, delimiter=',', comments='#')
                if dEA.ndim == 1: dEA = dEA[None, :]
                if dNA.ndim == 1: dNA = dNA[None, :]
                n = min(len(dEA), len(dNA)); dEA, dNA = dEA[:n], dNA[:n]
                JA = dEA[:, 6]; X = 1.0 - dEA[:, 0]; Y = 1.0 - dEA[:, 1]
                xo, uo, yo, vo = REF[type_]
                lro = radial_stretch(xo, uo, yo, vo)
                lr  = radial_stretch(dNA[:, 3], dNA[:, 4], dNA[:, 5], dNA[:, 6]) / lro
                eps = 1e-02
                m = (JA <= 99.0 + eps) & (JA >= 1.0 - eps)
                for i in np.where(m)[0]:
                    rows.append([carpeta, type_, nu_val, alph_val, zod_val,
                                 X[i], Y[i], lr[i], JA[i], nota])
            except Exception as e:
                print(f"error {carpeta} {type_}: {e}")
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['folder','type','z3_nu','alph','z4_zod','x_elastin_deg',
                    'y_collagen_deg','z1_lamb_r','z2_J','note'])
        w.writerows(rows)
    print(f"MASTER rebuilt: {len(rows)} data rows -> {OUT}")

if __name__ == "__main__":
    main()
