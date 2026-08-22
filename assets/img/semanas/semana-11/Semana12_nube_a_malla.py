# -*- coding: utf-8 -*-
"""
Semana 12 - Prototipado y Simulacion (DITA, Universidad La Salle Bajio)
Codigo 3. De la nube de puntos a la malla STL apta para simulacion.

Pipeline completo, auditable y con verdad de terreno analitica:

  1. Solido de referencia definido por una funcion de distancia con signo
     (SDF): brida circular + buje cilindrico con barreno pasante.
  2. Muestreo superficial y contaminacion controlada: ruido gaussiano normal
     de desviacion s_n y una fraccion f_out de outliers volumetricos.
  3. Filtrado estadistico de outliers (SOR): distancia media a los k vecinos
     mas proximos; se descarta el punto si  d_i > mu_d + n_sigma * sigma_d.
     Se reportan precision y exhaustividad frente al etiquetado real.
  4. Submuestreo por voxeles (rejilla de arista a_v) y estimacion de normales
     por analisis de componentes principales de la matriz de covarianza local
     (el autovector del autovalor minimo es la normal).
  5. Orientacion global de normales por propagacion de signo sobre el arbol
     de expansion minima del grafo kNN, con peso  w_ij = 1 - |n_i . n_j|.
  6. Reconstruccion implicita: campo escalar por minimos cuadrados moviles
     sobre planos tangentes locales, evaluado unicamente en la cascara de
     radio r_max en torno a la nube. El interior y el exterior se etiquetan
     por componentes conexas, evitando la cascara falsa caracteristica de
     Poisson sin recorte. Extraccion de la isosuperficie cero por Marching
     Cubes.
  7. Auditoria topologica: V, E, F, caracteristica de Euler, genero, aristas
     de frontera, aristas no-manifold, caras degeneradas, estanqueidad.
  8. Suavizado de Taubin (lambda / mu) que preserva volumen, decimacion por
     agrupamiento de vertices y exportacion a STL binario sin dependencias.
  9. Verificacion dimensional contra el SDF analitico: RMS, P95 y maximo de
     la desviacion, y error volumetrico.

Requisitos: numpy, scipy, scikit-image, matplotlib. Sin dependencias de GPU.
Salida: Figuras/S12_nube_a_malla.png (300 dpi), Figuras/S12_pieza_final.stl.
Backend Agg: ejecutable en Google Colab sin display.
"""

import os
import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import (minimum_spanning_tree,
                                  breadth_first_order,
                                  connected_components)
from scipy import ndimage
from skimage import measure

os.makedirs('Figuras', exist_ok=True)
rng = np.random.default_rng(1204)

# ----------------------------------------------------------------------
# Parametros del experimento (declarados)
# ----------------------------------------------------------------------
R_BRIDA, T_BRIDA = 55.0, 12.0     # radio y espesor de la brida [mm]
R_BUJE, H_BUJE = 25.0, 80.0       # radio y altura del buje [mm]
R_HOLE = 10.0                     # radio del barreno pasante [mm]
N_PTS = 120000                    # puntos muestreados en la superficie
S_N = 0.10                        # ruido normal, desviacion estandar [mm]
F_OUT = 0.02                      # fraccion de outliers
K_SOR, NSIG_SOR, IT_SOR = 20, 3.0, 2   # filtro estadistico robusto
A_VOX = 0.60                      # arista del voxel [mm]
K_NRM = 25                        # vecinos para PCA de normales
H_GRID = 0.70                     # paso de la rejilla implicita [mm]
R_MAX = 2.20                      # radio de la cascara de evaluacion [mm]
LAM_T, MU_T, IT_T = 0.50, -0.53, 3   # Taubin
A_CLUS = 0.90                     # arista de agrupamiento para decimacion


# ----------------------------------------------------------------------
# 1. Verdad de terreno: SDF analitica del solido
# ----------------------------------------------------------------------
def sdf_cilindro(P, r, z0, z1):
    """SDF exacta de un cilindro recto con tapas, eje Z."""
    zc, hz = 0.5 * (z0 + z1), 0.5 * (z1 - z0)
    dxy = np.hypot(P[..., 0], P[..., 1]) - r
    dz = np.abs(P[..., 2] - zc) - hz
    fuera = np.sqrt(np.maximum(dxy, 0.0) ** 2 + np.maximum(dz, 0.0) ** 2)
    return np.minimum(np.maximum(dxy, dz), 0.0) + fuera


def sdf_pieza(P):
    """Union de brida y buje, menos el barreno pasante. Centrada en Z."""
    P = np.asarray(P, dtype=float)
    z_off = (T_BRIDA + H_BUJE) / 2.0
    Q = P.copy()
    Q[..., 2] = Q[..., 2] + z_off
    d_bri = sdf_cilindro(Q, R_BRIDA, 0.0, T_BRIDA)
    d_buj = sdf_cilindro(Q, R_BUJE, T_BRIDA, T_BRIDA + H_BUJE)
    d_hol = sdf_cilindro(Q, R_HOLE, -1.0, T_BRIDA + H_BUJE + 1.0)
    return np.maximum(np.minimum(d_bri, d_buj), -d_hol)


VOL_GT = (np.pi * R_BRIDA ** 2 * T_BRIDA
          + np.pi * R_BUJE ** 2 * H_BUJE
          - np.pi * R_HOLE ** 2 * (T_BRIDA + H_BUJE))


def malla_de_sdf(func, lo, hi, paso):
    """Isosuperficie cero de un campo escalar por Marching Cubes."""
    ejes = [np.arange(lo[i], hi[i] + paso, paso) for i in range(3)]
    G = np.stack(np.meshgrid(*ejes, indexing='ij'), axis=-1)
    F = func(G)
    v, f, _, _ = measure.marching_cubes(F, level=0.0,
                                        spacing=(paso, paso, paso))
    v = v + np.array(lo)
    return v, f


def componente_mayor(V, F):
    """Conserva la pieza conexa de mayor numero de vertices."""
    e = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])),
                   shape=(len(V),) * 2)
    n_c, et = connected_components(A, directed=False)
    keep = np.argmax(np.bincount(et))
    mv = et == keep
    remap = -np.ones(len(V), dtype=np.int64)
    remap[mv] = np.arange(mv.sum())
    mf = mv[F].all(1)
    return V[mv], remap[F[mf]], n_c


def orientar_hacia_afuera(V, F):
    """Fuerza volumen positivo por el teorema de la divergencia."""
    tri = V[F]
    vol = np.einsum('ij,ij->i', np.cross(tri[:, 0], tri[:, 1]),
                    tri[:, 2]).sum() / 6.0
    return F if vol > 0 else F[:, ::-1]


# --- muestreo superficial ponderado por area de las caras ---------------
V_GT, F_GT = malla_de_sdf(sdf_pieza, (-62, -62, -55), (62, 62, 55), 0.45)
tri = V_GT[F_GT]
nrm_f = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
area = 0.5 * np.linalg.norm(nrm_f, axis=1)
val = area > 1e-9                       # descartar caras degeneradas
F_GT, tri, nrm_f, area = F_GT[val], tri[val], nrm_f[val], area[val]
nrm_f = nrm_f / np.linalg.norm(nrm_f, axis=1, keepdims=True)
prob = area / area.sum()
sel = rng.choice(len(F_GT), N_PTS, p=prob)
u, v_ = rng.random(N_PTS), rng.random(N_PTS)
fix = u + v_ > 1.0
u[fix], v_[fix] = 1.0 - u[fix], 1.0 - v_[fix]
P_surf = (tri[sel, 0] + u[:, None] * (tri[sel, 1] - tri[sel, 0])
          + v_[:, None] * (tri[sel, 2] - tri[sel, 0]))
N_surf = nrm_f[sel]

# --- contaminacion: ruido normal + outliers volumetricos ---------------
n_out = int(F_OUT * N_PTS)
P_ruido = P_surf + (rng.normal(0.0, S_N, N_PTS))[:, None] * N_surf
lo_b, hi_b = P_surf.min(0) - 8.0, P_surf.max(0) + 8.0
P_out = rng.uniform(lo_b, hi_b, size=(n_out, 3))
P_all = np.vstack([P_ruido, P_out])
es_outlier = np.zeros(len(P_all), dtype=bool)
es_outlier[N_PTS:] = True
orden = rng.permutation(len(P_all))
P_all, es_outlier = P_all[orden], es_outlier[orden]

print('=' * 72)
print('PIPELINE NUBE -> MALLA -> STL')
print('=' * 72)
print('Nube cruda: {:d} puntos ({:d} outliers, {:.1f} por ciento); '
      'ruido normal s_n = {:.2f} mm'
      .format(len(P_all), n_out, 100.0 * n_out / len(P_all), S_N))

# ----------------------------------------------------------------------
# 2. Filtrado estadistico de outliers (SOR)
# ----------------------------------------------------------------------
mant = np.ones(len(P_all), dtype=bool)
for _ in range(IT_SOR):
    kd = cKDTree(P_all[mant])
    d_knn, _ = kd.query(P_all[mant], k=K_SOR + 1)
    d_med = d_knn[:, 1:].mean(1)
    med = np.median(d_med)
    mad = 1.4826 * np.median(np.abs(d_med - med))   # escala robusta
    umbral = med + NSIG_SOR * mad
    sub_ok = d_med <= umbral
    idx_act = np.flatnonzero(mant)
    mant[idx_act[~sub_ok]] = False
TP = int(np.sum(es_outlier & ~mant))
FP = int(np.sum(~es_outlier & ~mant))
FN = int(np.sum(es_outlier & mant))
PREC = 100.0 * TP / max(TP + FP, 1)
REC = 100.0 * TP / max(TP + FN, 1)
P_f = P_all[mant]
print('SOR robusto (k = {:d}, n_sigma = {:.1f}, {:d} pasadas): umbral final '
      'd = {:.4f} mm; elimina {:d} puntos'
      .format(K_SOR, NSIG_SOR, IT_SOR, umbral, int(np.sum(~mant))))
print('     exhaustividad sobre outliers = {:.1f} por ciento; '
      'precision = {:.1f} por ciento; supervivientes = {:d}'
      .format(REC, PREC, len(P_f)))

# ----------------------------------------------------------------------
# 3. Submuestreo por voxeles
# ----------------------------------------------------------------------
idx_v = np.floor((P_f - P_f.min(0)) / A_VOX).astype(np.int64)
clave = (idx_v[:, 0] * 100003 + idx_v[:, 1]) * 100003 + idx_v[:, 2]
orden = np.argsort(clave, kind='stable')
clave_o, P_o = clave[orden], P_f[orden]
corte = np.flatnonzero(np.diff(clave_o)) + 1
grupos = np.split(np.arange(len(P_o)), corte)
P_v = np.array([P_o[g].mean(0) for g in grupos])
print('Voxel grid (a_v = {:.2f} mm): {:d} -> {:d} puntos '
      '(factor {:.2f})'.format(A_VOX, len(P_f), len(P_v),
                               len(P_v) / len(P_f)))

# ----------------------------------------------------------------------
# 4. Normales por PCA local
# ----------------------------------------------------------------------
kd_v = cKDTree(P_v)
_, nb = kd_v.query(P_v, k=K_NRM)
Q = P_v[nb] - P_v[nb].mean(1, keepdims=True)
C = np.einsum('nki,nkj->nij', Q, Q) / K_NRM
w, U = np.linalg.eigh(C)
N_v = U[:, :, 0]
curv = w[:, 0] / np.maximum(w.sum(1), 1e-12)     # variacion de superficie

# ----------------------------------------------------------------------
# 5. Orientacion global por MST del grafo kNN
# ----------------------------------------------------------------------
K_G = 12
_, nb_g = kd_v.query(P_v, k=K_G)
fil = np.repeat(np.arange(len(P_v)), K_G - 1)
col = nb_g[:, 1:].ravel()
peso = 1.0 - np.abs(np.einsum('ij,ij->i', N_v[fil], N_v[col]))
G = coo_matrix((peso + 1e-6, (fil, col)), shape=(len(P_v),) * 2)
T = minimum_spanning_tree(G)
T = T + T.T
semilla = int(np.argmax(P_v[:, 2]))
if N_v[semilla, 2] < 0:
    N_v[semilla] *= -1.0
orden_bfs, padres = breadth_first_order(T, semilla, directed=False)
for i in orden_bfs[1:]:
    p = padres[i]
    if p >= 0 and np.dot(N_v[i], N_v[p]) < 0:
        N_v[i] *= -1.0
coher = float(np.mean(np.einsum('ij,ij->i', N_v[fil], N_v[col]) > 0))
print('Normales: PCA con k = {:d}; coherencia de orientacion en el grafo '
      'kNN = {:.1f} por ciento'.format(K_NRM, 100.0 * coher))

# ----------------------------------------------------------------------
# 6. Campo implicito y Marching Cubes
# ----------------------------------------------------------------------
lo = P_v.min(0) - 3.0
hi = P_v.max(0) + 3.0
ejes = [np.arange(lo[i], hi[i] + H_GRID, H_GRID) for i in range(3)]
forma = tuple(len(e) for e in ejes)
G3 = np.stack(np.meshgrid(*ejes, indexing='ij'), axis=-1).reshape(-1, 3)

d1, _ = kd_v.query(G3, k=1)
cerca = d1 <= R_MAX
K_IM = 12
_, nb_im = kd_v.query(G3[cerca], k=K_IM)
dif = G3[cerca][:, None, :] - P_v[nb_im]
dist = np.linalg.norm(dif, axis=2)
w_im = np.exp(-(dist / (0.6 * R_MAX)) ** 2)
proy = np.einsum('nkj,nkj->nk', dif, N_v[nb_im])
f_shell = (w_im * proy).sum(1) / np.maximum(w_im.sum(1), 1e-12)

F3 = np.full(len(G3), np.nan)
F3[cerca] = f_shell
F3 = F3.reshape(forma)

lejos = np.isnan(F3)
lab, nlab = ndimage.label(lejos)
borde = set(np.unique(np.concatenate([
    lab[0].ravel(), lab[-1].ravel(), lab[:, 0].ravel(), lab[:, -1].ravel(),
    lab[:, :, 0].ravel(), lab[:, :, -1].ravel()])))
borde.discard(0)
ext = np.isin(lab, list(borde))
F3[lejos & ext] = +10.0                 # exterior conexo con la frontera
F3[lejos & ~ext] = -10.0                # cavidades internas del solido
print('Campo implicito: rejilla {} = {:.2f} M nodos; cascara evaluada = '
      '{:.1f} por ciento; componentes lejanas = {:d}'
      .format(forma, F3.size / 1e6, 100.0 * cerca.mean(), nlab))

V_m, F_m, _, _ = measure.marching_cubes(F3, level=0.0,
                                        spacing=(H_GRID,) * 3)
V_m = V_m + lo
F_MC0 = len(F_m)
V_m, F_m, N_COMP = componente_mayor(V_m, F_m)
F_m = orientar_hacia_afuera(V_m, F_m)
print('Marching Cubes bruto: {:d} caras en {:d} piezas conexas; se conserva '
      'la mayor con {:d} caras ({:.1f} por ciento)'
      .format(F_MC0, N_COMP, len(F_m), 100.0 * len(F_m) / F_MC0))


# ----------------------------------------------------------------------
# 7. Auditoria topologica
# ----------------------------------------------------------------------
def auditar(V, F):
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    tri = V[F]
    ar = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0],
                                       tri[:, 2] - tri[:, 0]), axis=1)
    c = V.mean(0)
    t0 = tri - c                       # origen en el centroide (robustez)
    vol = np.abs(np.einsum('ij,ij->i', np.cross(t0[:, 0], t0[:, 1]),
                           t0[:, 2]).sum() / 6.0)
    return dict(V=len(V), E=len(cnt), F=len(F),
                chi=len(V) - len(cnt) + len(F),
                frontera=int(np.sum(cnt == 1)),
                nomanifold=int(np.sum(cnt > 2)),
                degeneradas=int(np.sum(ar < 1e-9)),
                area=float(ar.sum()), volumen=float(vol))


aud0 = auditar(V_m, F_m)
print('-' * 72)
print('Marching Cubes: V = {V:d}, E = {E:d}, F = {F:d}, chi = {chi:d}, '
      'genero = {g:d}'.format(g=(2 - aud0['chi']) // 2, **aud0))
print('     aristas de frontera = {frontera:d}, no-manifold = {nomanifold:d}, '
      'degeneradas = {degeneradas:d}'.format(**aud0))
print('     estanca = {}'.format(aud0['frontera'] == 0
                                 and aud0['nomanifold'] == 0))


# ----------------------------------------------------------------------
# 8. Suavizado de Taubin, decimacion y exportacion
# ----------------------------------------------------------------------
def laplaciano_uniforme(V, F):
    e = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.vstack([e, e[:, ::-1]])
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])),
                   shape=(len(V),) * 2).tocsr()
    A.data[:] = 1.0
    gr = np.asarray(A.sum(1)).ravel()
    return (A @ V) / np.maximum(gr, 1)[:, None] - V


def taubin(V, F, lam=LAM_T, mu=MU_T, it=IT_T):
    Vs = V.copy()
    for _ in range(it):
        Vs = Vs + lam * laplaciano_uniforme(Vs, F)
        Vs = Vs + mu * laplaciano_uniforme(Vs, F)
    return Vs


def decimar_por_agrupamiento(V, F, a):
    idx = np.floor((V - V.min(0)) / a).astype(np.int64)
    _, inv = np.unique(idx, axis=0, return_inverse=True)
    n = inv.max() + 1
    Vn = np.zeros((n, 3))
    cnt = np.bincount(inv, minlength=n)[:, None]
    np.add.at(Vn, inv, V)
    Vn /= np.maximum(cnt, 1)
    Fn = inv[F]
    ok = ((Fn[:, 0] != Fn[:, 1]) & (Fn[:, 1] != Fn[:, 2])
          & (Fn[:, 0] != Fn[:, 2]))
    Fn = Fn[ok]
    # rotacion canonica: conserva la orientacion y permite deduplicar
    r = np.argmin(Fn, axis=1)
    Fn = np.take_along_axis(Fn, (r[:, None] + np.arange(3)[None, :]) % 3, 1)
    _, uniq = np.unique(Fn, axis=0, return_index=True)
    Fn = Fn[np.sort(uniq)]
    usados, Fn = np.unique(Fn, return_inverse=True)
    return Vn[usados], Fn.reshape(-1, 3)


def escribir_stl_binario(ruta, V, F):
    tri = V[F].astype(np.float32)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    with open(ruta, 'wb') as fh:
        fh.write(b'Semana12 DITA malla fotogrametrica'.ljust(80, b' '))
        fh.write(struct.pack('<I', len(F)))
        for k in range(len(F)):
            fh.write(struct.pack('<12fH', *n[k], *tri[k, 0], *tri[k, 1],
                                 *tri[k, 2], 0))
    return os.path.getsize(ruta)


V_s = taubin(V_m, F_m)
aud1 = auditar(V_s, F_m)
V_d, F_d = decimar_por_agrupamiento(V_s, F_m, A_CLUS)
V_d, F_d, _ = componente_mayor(V_d, F_d)
F_d = orientar_hacia_afuera(V_d, F_d)
aud2 = auditar(V_d, F_d)
BYTES = escribir_stl_binario('Figuras/S12_pieza_final.stl', V_s, F_m)

print('-' * 72)
print('Taubin (lambda = {:.2f}, mu = {:.2f}, {:d} pasadas): cambio de '
      'volumen = {:+.3f} por ciento'
      .format(LAM_T, MU_T, IT_T,
              100.0 * (aud1['volumen'] / aud0['volumen'] - 1.0)))
print('Decimacion por agrupamiento (a = {:.2f} mm): {:d} -> {:d} caras '
      '({:.1f} por ciento)'.format(A_CLUS, aud0['F'], aud2['F'],
                                   100.0 * aud2['F'] / aud0['F']))
print('     coste topologico: chi = {chi:d}, frontera = {frontera:d}, '
      'aristas no-manifold = {nomanifold:d}'.format(**aud2))
print('     ADVERTENCIA: el agrupamiento de vertices no preserva la variedad;')
print('     para el entregable se decima con QEM (MeshLab / Open3D) y se')
print('     reaudita. Se exporta la malla suavizada, estanca y de genero 1.')
print('STL binario exportado (malla suavizada): {:.2f} MB, {:d} caras'
      .format(BYTES / 1e6, aud1['F']))

# ----------------------------------------------------------------------
# 9. Verificacion dimensional contra el SDF analitico
# ----------------------------------------------------------------------
dev_mc = sdf_pieza(V_m)
dev_fin = sdf_pieza(V_s)
RMS_MC = float(np.sqrt(np.mean(dev_mc ** 2)))
RMS_FIN = float(np.sqrt(np.mean(dev_fin ** 2)))
P95 = float(np.percentile(np.abs(dev_fin), 95))
MAXD = float(np.max(np.abs(dev_fin)))
BIAS = float(np.mean(dev_fin))
E_VOL = 100.0 * (aud1['volumen'] / VOL_GT - 1.0)

print('-' * 72)
print('VERIFICACION CONTRA LA VERDAD DE TERRENO')
print('     RMS tras Marching Cubes      = {:.4f} mm'.format(RMS_MC))
print('     RMS de la malla entregable   = {:.4f} mm'.format(RMS_FIN))
print('     sesgo medio                  = {:+.4f} mm'.format(BIAS))
print('     P95 de |desviacion|          = {:.4f} mm'.format(P95))
print('     maximo de |desviacion|       = {:.4f} mm'.format(MAXD))
print('     volumen GT = {:.1f} mm3; malla = {:.1f} mm3; error = '
      '{:+.3f} por ciento'.format(VOL_GT, aud1['volumen'], E_VOL))
print('     area superficial de la malla = {:.1f} mm2'.format(aud1['area']))
print('     s_n de entrada = {:.3f} mm -> reduccion efectiva = {:.2f}x'
      .format(S_N, S_N / RMS_FIN))
print('=' * 72)

# ----------------------------------------------------------------------
# 10. Figura de tres paneles
# ----------------------------------------------------------------------
AZUL = '#003366'
ROJO = '#B22222'
VERDE = '#2E7D32'

fig = plt.figure(figsize=(15.0, 4.6))

ax0 = fig.add_subplot(1, 3, 1, projection='3d')
sub = rng.choice(len(P_all), 5000, replace=False)
inl = sub[~es_outlier[sub]]
out = sub[es_outlier[sub]]
ax0.scatter(P_all[inl, 0], P_all[inl, 1], P_all[inl, 2], s=0.8,
            c=P_all[inl, 2], cmap='viridis', alpha=0.55)
ax0.scatter(P_all[out, 0], P_all[out, 1], P_all[out, 2], s=9,
            color=ROJO, marker='x', lw=0.8, label='outliers inyectados')
ax0.set_title('(a) Nube cruda: ruido $\\sigma_n$ = {:.2f} mm + '
              '{:.0f} por ciento de outliers'.format(S_N, 100 * F_OUT),
              fontsize=9.5, loc='left')
ax0.set_xlabel('X [mm]', fontsize=8)
ax0.set_ylabel('Y [mm]', fontsize=8)
ax0.set_zlabel('Z [mm]', fontsize=8)
ax0.tick_params(labelsize=7)
ax0.legend(fontsize=7.5, frameon=False, loc='upper right')
ax0.view_init(elev=18, azim=35)

ax1 = fig.add_subplot(1, 3, 2, projection='3d')
V_p, F_p = decimar_por_agrupamiento(V_s, F_m, 2.20)   # solo para render
V_p, F_p, _ = componente_mayor(V_p, F_p)
tri_p = V_p[F_p]
col = plt.cm.viridis((tri_p[:, :, 2].mean(1) - V_s[:, 2].min())
                     / np.ptp(V_s[:, 2]))
pc = Poly3DCollection(tri_p, facecolors=col, edgecolors='k',
                      linewidths=0.12, alpha=1.0)
pc.set_sort_zpos(0)
ax1.add_collection3d(pc)
ax1.set_xlim(V_s[:, 0].min(), V_s[:, 0].max())
ax1.set_ylim(V_s[:, 1].min(), V_s[:, 1].max())
ax1.set_zlim(V_s[:, 2].min(), V_s[:, 2].max())
ax1.set_box_aspect((1, 1, np.ptp(V_s[:, 2]) / np.ptp(V_s[:, 0])))
ax1.set_title('(b) Malla entregable: {:d} caras, estanca, $\\chi$ = {:d} '
              '(genero 1)\n(vista simplificada a {:d} caras)'
              .format(aud1['F'], aud1['chi'], len(F_p)),
              fontsize=9.5, loc='left')
ax1.set_xlabel('X [mm]', fontsize=8)
ax1.set_ylabel('Y [mm]', fontsize=8)
ax1.set_zlabel('Z [mm]', fontsize=8)
ax1.tick_params(labelsize=7)
ax1.view_init(elev=18, azim=35)

ax2 = fig.add_subplot(1, 3, 3)
ax2.hist(dev_fin * 1000, bins=90, color=AZUL, alpha=0.82,
         edgecolor='none', density=True)
ax2.axvline(0.0, color='k', lw=0.8)
for x, c, lab in [(RMS_FIN * 1000, ROJO, 'RMS = {:.0f} $\\mu$m'
                   .format(RMS_FIN * 1000)),
                  (-RMS_FIN * 1000, ROJO, None),
                  (P95 * 1000, VERDE, 'P95 = {:.0f} $\\mu$m'
                   .format(P95 * 1000)),
                  (-P95 * 1000, VERDE, None)]:
    ax2.axvline(x, color=c, ls='--', lw=1.3, label=lab)
ax2.set_xlim(-300, 300)
ax2.set_xlabel('Desviacion malla - solido analitico [$\\mu$m]')
ax2.set_ylabel('Densidad [1/$\\mu$m]')
ax2.set_title('(c) Verificacion dimensional (sesgo = {:+.0f} $\\mu$m)'
              .format(BIAS * 1000), fontsize=9.5, loc='left')
ax2.legend(fontsize=8, frameon=False)
ax2.grid(alpha=0.3)

fig.tight_layout()
fig.savefig('Figuras/S12_nube_a_malla.png', dpi=300, bbox_inches='tight')
# plt.show()
print('Figura guardada: Figuras/S12_nube_a_malla.png')
