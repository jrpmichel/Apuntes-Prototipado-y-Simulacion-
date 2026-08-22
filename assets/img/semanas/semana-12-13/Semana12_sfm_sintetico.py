# -*- coding: utf-8 -*-
"""
Semana 12 - Prototipado y Simulacion (DITA, Universidad La Salle Bajio)
Codigo 2. Experimento sintetico de Structure from Motion.

Objetivo: cuantificar, sobre un banco sintetico con verdad de terreno exacta,
las dos propiedades del pipeline SfM que condicionan su uso metrologico.

  (A) Escalamiento del error de triangulacion con el numero de vistas.
      Se proyecta una pieza sintetica sobre N camaras dispuestas en dos
      anillos orbitales, se contamina la deteccion con ruido gaussiano de
      s_px [px], se triangula por DLT multivista (SVD) y se compara contra
      la verdad de terreno tras alineamiento de Umeyama. El escalamiento
      empirico se ajusta a  RMS ~ N^(-q)  y se contrasta con q = 1/2.

  (B) Propagacion de la escala metrica. La reconstruccion SfM queda definida
      salvo una similitud (7 gdl). La escala se fija con una barra calibrada
      de longitud L. Si cada extremo de la barra tiene incertidumbre 3D
      isotropa u_p, la escala hereda
              u_s / s = sqrt(2) * u_p / L,
      es decir, el error relativo del modelo completo es inversamente
      proporcional a la longitud de la barra. Se verifica por Monte Carlo.

Salida: Figuras/S12_sfm_sintetico.png (300 dpi) y resumen en consola.
Backend Agg: ejecutable en Google Colab sin display.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.ticker import NullFormatter, NullLocator

os.makedirs('Figuras', exist_ok=True)
rng = np.random.default_rng(2024)

# ----------------------------------------------------------------------
# 1. Camara: mismos intrinsecos del Codigo 1 (APS-C 24 MP, f = 35 mm)
# ----------------------------------------------------------------------
NPX_H, NPX_V = 6000, 4000
F_PX = 9417.0
CX, CY = NPX_H / 2.0, NPX_V / 2.0
K = np.array([[F_PX, 0.0, CX],
              [0.0, F_PX, CY],
              [0.0, 0.0, 1.0]])
S_PX = 0.30        # ruido de deteccion de features [px]
R_ORB = 600.0      # radio orbital [mm]
ELEV = [0.0, 30.0]  # anillos de captura [deg]


# ----------------------------------------------------------------------
# 2. Pieza sintetica: buje cilindrico sobre brida circular
# ----------------------------------------------------------------------
def pieza_sintetica(n=4000):
    """Devuelve (puntos [n,3] en mm, normales unitarias [n,3])."""
    r_buje, h_buje = 25.0, 80.0
    r_brida, t_brida = 55.0, 12.0
    n1 = int(0.45 * n)   # lateral del buje
    n2 = int(0.30 * n)   # canto y cara superior de la brida
    n3 = n - n1 - n2     # tapa del buje

    th = rng.uniform(0, 2 * np.pi, n1)
    z = rng.uniform(t_brida, t_brida + h_buje, n1)
    p1 = np.c_[r_buje * np.cos(th), r_buje * np.sin(th), z]
    nn1 = np.c_[np.cos(th), np.sin(th), np.zeros(n1)]

    th = rng.uniform(0, 2 * np.pi, n2)
    rad = np.sqrt(rng.uniform(r_buje ** 2, r_brida ** 2, n2))
    p2 = np.c_[rad * np.cos(th), rad * np.sin(th),
               np.full(n2, t_brida)]
    nn2 = np.tile(np.array([0.0, 0.0, 1.0]), (n2, 1))

    th = rng.uniform(0, 2 * np.pi, n3)
    rad = r_buje * np.sqrt(rng.uniform(0, 1, n3))
    p3 = np.c_[rad * np.cos(th), rad * np.sin(th),
               np.full(n3, t_brida + h_buje)]
    nn3 = np.tile(np.array([0.0, 0.0, 1.0]), (n3, 1))

    P = np.vstack([p1, p2, p3])
    Nrm = np.vstack([nn1, nn2, nn3])
    P[:, 2] -= (t_brida + h_buje) / 2.0        # centrar en Z
    return P, Nrm / np.linalg.norm(Nrm, axis=1, keepdims=True)


def poses_orbitales(n_cam, radio=R_ORB, elevaciones=ELEV):
    """Genera matrices [R|t] mirando al origen desde anillos orbitales."""
    poses, centros = [], []
    por_anillo = int(np.ceil(n_cam / len(elevaciones)))
    for el in elevaciones:
        for k in range(por_anillo):
            az = 2 * np.pi * k / por_anillo
            e = np.deg2rad(el)
            C = radio * np.array([np.cos(e) * np.cos(az),
                                  np.cos(e) * np.sin(az),
                                  np.sin(e)])
            zc = -C / np.linalg.norm(C)                    # eje optico
            up = np.array([0.0, 0.0, 1.0])
            if abs(np.dot(zc, up)) > 0.98:
                up = np.array([0.0, 1.0, 0.0])
            xc = np.cross(up, zc)
            xc /= np.linalg.norm(xc)
            yc = np.cross(zc, xc)
            R = np.vstack([xc, yc, zc])
            poses.append((R, -R @ C))
            centros.append(C)
    return poses[:n_cam], np.array(centros[:n_cam])


def proyectar(P, R, t):
    Pc = (R @ P.T).T + t
    uv = (K @ Pc.T).T
    return uv[:, :2] / uv[:, 2:3], Pc[:, 2]


def triangular_dlt(obs, Ps):
    """obs: lista de (u,v); Ps: lista de matrices 3x4. Devuelve X en mm."""
    A = []
    for (u, v), Pm in zip(obs, Ps):
        A.append(u * Pm[2] - Pm[0])
        A.append(v * Pm[2] - Pm[1])
    _, _, Vt = np.linalg.svd(np.array(A))
    X = Vt[-1]
    return X[:3] / X[3]


def umeyama(X, Y):
    """Similitud (s,R,t) que lleva X sobre Y en minimos cuadrados."""
    mx, my = X.mean(0), Y.mean(0)
    Xc, Yc = X - mx, Y - my
    S = (Yc.T @ Xc) / len(X)
    U, D, Vt = np.linalg.svd(S)
    d = np.sign(np.linalg.det(U @ Vt))
    W = np.diag([1.0, 1.0, d])
    R = U @ W @ Vt
    s = np.trace(np.diag(D) @ W) / (Xc ** 2).sum(1).mean()
    return s, R, my - s * (R @ mx)


def reconstruir(P_gt, Nrm, n_cam, s_px=S_PX, min_vistas=3):
    """Devuelve (X_rec, idx_validos, rms_reproyeccion)."""
    poses, _ = poses_orbitales(n_cam)
    Pmats = [K @ np.hstack([R, t.reshape(3, 1)]) for R, t in poses]
    obs = {i: [] for i in range(len(P_gt))}
    res_rep = []
    for j, (R, t) in enumerate(poses):
        C = -R.T @ t
        uv, zc = proyectar(P_gt, R, t)
        vis = ((Nrm * (C - P_gt)).sum(1) > 0) & (zc > 0)
        vis &= (uv[:, 0] > 0) & (uv[:, 0] < NPX_H)
        vis &= (uv[:, 1] > 0) & (uv[:, 1] < NPX_V)
        ruido = rng.normal(0.0, s_px, size=uv.shape)
        for i in np.flatnonzero(vis):
            obs[i].append((uv[i] + ruido[i], j))
    X, idx = [], []
    for i, lst in obs.items():
        if len(lst) >= min_vistas:
            X.append(triangular_dlt([o[0] for o in lst],
                                    [Pmats[o[1]] for o in lst]))
            idx.append(i)
            for uvo, j in lst:
                R, t = poses[j]
                uvp, _ = proyectar(X[-1][None, :], R, t)
                res_rep.extend((uvp[0] - uvo).tolist())
    return np.array(X), np.array(idx), float(np.sqrt(np.mean(
        np.array(res_rep) ** 2)))


# ----------------------------------------------------------------------
# 3. (A) Error de triangulacion frente al numero de vistas
# ----------------------------------------------------------------------
P_gt, N_gt = pieza_sintetica(4000)
N_LIST = [6, 8, 12, 16, 24, 36, 48]
rms_3d, rms_rep, n_pts, cob = [], [], [], []

for n_cam in N_LIST:
    Xr, idx, rrep = reconstruir(P_gt, N_gt, n_cam)
    s, R, t = umeyama(Xr, P_gt[idx])
    Xa = s * (R @ Xr.T).T + t
    d = np.linalg.norm(Xa - P_gt[idx], axis=1)
    rms_3d.append(float(np.sqrt(np.mean(d ** 2))))
    rms_rep.append(rrep)
    n_pts.append(len(idx))
    cob.append(100.0 * len(idx) / len(P_gt))

rms_3d = np.array(rms_3d)
coef = np.polyfit(np.log(N_LIST), np.log(rms_3d), 1)
Q_EMP = -coef[0]

# Reconstruccion de referencia con 24 camaras (punto de diseno del curso)
X24, idx24, rrep24 = reconstruir(P_gt, N_gt, 24)
s24, R24, t24 = umeyama(X24, P_gt[idx24])
X24a = s24 * (R24 @ X24.T).T + t24
d24 = np.linalg.norm(X24a - P_gt[idx24], axis=1)
U_P = float(np.sqrt(np.mean(d24 ** 2)))     # incertidumbre 3D por punto [mm]

print('=' * 70)
print('(A) TRIANGULACION MULTIVISTA')
print('=' * 70)
print('{:>6} {:>10} {:>12} {:>14} {:>14}'
      .format('N cam', 'pts', 'cobertura', 'RMS 3D [mm]', 'RMS rep [px/c]'))
for n_cam, npt, c, r3, rr in zip(N_LIST, n_pts, cob, rms_3d, rms_rep):
    print('{:>6d} {:>10d} {:>11.1f}% {:>14.4f} {:>14.3f}'
          .format(n_cam, npt, c, r3, rr))
print('-' * 70)
print('Exponente empirico  RMS ~ N^(-q):  q = {:.3f}  (teorico 0.5)'
      .format(Q_EMP))
print('Punto de diseno N = 24: RMS 3D = {:.4f} mm, '
      'residual de reproyeccion = {:.3f} px/coord'.format(U_P, rrep24))
print('Relacion RMS 3D / GSD(600 mm, 63.7 um) = {:.2f}'.format(U_P / 0.0637))

# ----------------------------------------------------------------------
# 4. (B) Escala metrica: Monte Carlo sobre la barra de calibracion
# ----------------------------------------------------------------------
L_LIST = np.array([25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0])
N_MC = 4000
err_s = []
for L in L_LIST:
    # dos targets separados L, cada extremo con incertidumbre 3D isotropa U_P
    e1 = rng.normal(0.0, U_P / np.sqrt(3.0), size=(N_MC, 3))
    e2 = rng.normal(0.0, U_P / np.sqrt(3.0), size=(N_MC, 3))
    p1 = np.zeros((N_MC, 3)) + e1
    p2 = np.tile([L, 0.0, 0.0], (N_MC, 1)) + e2
    L_med = np.linalg.norm(p2 - p1, axis=1)
    err_s.append(float(np.std(L_med / L)))       # u_s/s adimensional
err_s = np.array(err_s)
err_teo = np.sqrt(2.0) * (U_P / np.sqrt(3.0)) / L_LIST

L_REF = 150.0
i_ref = int(np.argmin(np.abs(L_LIST - L_REF)))
DIM_NOM = 110.0        # cota nominal de la pieza [mm]
print('=' * 70)
print('(B) ESCALA METRICA A PARTIR DE BARRA CALIBRADA')
print('=' * 70)
print('{:>10} {:>16} {:>18} {:>18}'
      .format('L [mm]', 'u_s/s [ppm]', 'u_s/s MC [ppm]',
              'error en 110 mm [um]'))
for L, em, et in zip(L_LIST, err_s, err_teo):
    print('{:>10.0f} {:>16.0f} {:>18.0f} {:>18.1f}'
          .format(L, et * 1e6, em * 1e6, em * DIM_NOM * 1000))
print('-' * 70)
print('Barra de {:.0f} mm -> u_s/s = {:.0f} ppm; sobre una cota de {:.0f} mm '
      'aporta {:.1f} um'.format(L_LIST[i_ref], err_s[i_ref] * 1e6, DIM_NOM,
                                err_s[i_ref] * DIM_NOM * 1000))
print('=' * 70)

# ----------------------------------------------------------------------
# 5. Figura de tres paneles
# ----------------------------------------------------------------------
AZUL = '#003366'
ROJO = '#B22222'
VERDE = '#2E7D32'

fig = plt.figure(figsize=(15.0, 4.5))

ax0 = fig.add_subplot(1, 3, 1, projection='3d')
sub = rng.choice(len(P_gt), 1200, replace=False)
ax0.scatter(P_gt[sub, 0], P_gt[sub, 1], P_gt[sub, 2], s=1.5,
            c=P_gt[sub, 2], cmap='viridis', alpha=0.85)
_, C24 = poses_orbitales(24)
ax0.scatter(C24[:, 0], C24[:, 1], C24[:, 2], marker='^', s=26,
            color=ROJO, depthshade=False, label='estaciones (N = 24)')
for C in C24:
    ax0.plot([C[0], 0], [C[1], 0], [C[2], 0], color=ROJO, lw=0.35, alpha=0.35)
ax0.set_xlabel('X [mm]', fontsize=8)
ax0.set_ylabel('Y [mm]', fontsize=8)
ax0.set_zlabel('Z [mm]', fontsize=8)
ax0.tick_params(labelsize=7)
ax0.set_title('(a) Banco sintetico: pieza y red de estaciones',
              fontsize=10, loc='left')
ax0.legend(fontsize=7.5, loc='upper right', frameon=False)
ax0.view_init(elev=20, azim=40)

ax1 = fig.add_subplot(1, 3, 2)
ax1.loglog(N_LIST, rms_3d * 1000, 'o-', color=AZUL, lw=1.9, ms=6,
           label='RMS 3D simulado')
ref = rms_3d[0] * (np.array(N_LIST, dtype=float) / N_LIST[0]) ** (-0.5)
ax1.loglog(N_LIST, ref * 1000, '--', color='gray', lw=1.4,
           label='referencia $N^{-1/2}$')
ax1.loglog([24], [U_P * 1000], 's', ms=9, mfc='none', mec=ROJO, mew=1.8)
ax1.annotate('punto de diseno\nN = 24: {:.0f} $\\mu$m'.format(U_P * 1000),
             xy=(24, U_P * 1000), xytext=(25, rms_3d[0] * 1000 * 0.95),
             fontsize=8, color=ROJO,
             arrowprops=dict(arrowstyle='->', color=ROJO, lw=0.9))
ax1.set_xlabel('Numero de estaciones $N$')
ax1.set_ylabel('RMS del error 3D [$\\mu$m]')
ax1.set_title('(b) Escalamiento del error de triangulacion ($q$ = {:.2f})'
              .format(Q_EMP), fontsize=10, loc='left')
ax1.xaxis.set_minor_locator(NullLocator())
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.set_xticks(N_LIST)
ax1.set_xticklabels([str(v) for v in N_LIST])
ax1.legend(fontsize=8, frameon=False)
ax1.grid(alpha=0.3, which='both')

ax2 = fig.add_subplot(1, 3, 3)
ax2.loglog(L_LIST, err_s * 1e6, 'o', color=VERDE, ms=6,
           label='Monte Carlo ({:d} realizaciones)'.format(N_MC))
ax2.loglog(L_LIST, err_teo * 1e6, '-', color=AZUL, lw=1.8,
           label='$\\sqrt{2}\\,u_p/L$')
ax2.axhline(1000, color=ROJO, ls=':', lw=1.2)
ax2.set_ylim(10, 2500)
ax2.text(26, 1150, 'criterio $10^{3}$ ppm (0.1 por ciento)',
         color=ROJO, fontsize=7.5)
ax2.set_xlabel('Longitud de la barra de escala $L$ [mm]')
ax2.set_ylabel('$u_s/s$ [ppm]')
ax2.set_title('(c) Herencia de la escala metrica', fontsize=10, loc='left')
ax2.legend(fontsize=8, frameon=False)
ax2.grid(alpha=0.3, which='both')

fig.tight_layout()
fig.savefig('Figuras/S12_sfm_sintetico.png', dpi=300, bbox_inches='tight')
# plt.show()
print('Figura guardada: Figuras/S12_sfm_sintetico.png')
