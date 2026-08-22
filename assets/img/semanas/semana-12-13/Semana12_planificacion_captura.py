# -*- coding: utf-8 -*-
"""
Semana 12 - Prototipado y Simulacion (DITA, Universidad La Salle Bajio)
Codigo 1. Presupuesto optico de una sesion fotogrametrica.

Calcula, para una configuracion camara-objeto declarada:
  (1) GSD  = p * Z / f                      (resolucion lateral en el objeto)
  (2) sZ   = Z^2 * sd / (fpx * B)           (incertidumbre de profundidad por
                                             triangulacion, par convergente)
  (3) B    = 2 * Z * sin(da/2)              (base efectiva en captura orbital)
  (4) o    = 1 - B / W,  W = 2*Z*tan(FOV/2) (solapamiento entre estaciones)
  (5) N    = 360/da * n_anillos             (numero de imagenes)

El modelo de matching degrada la precision de correspondencia con la
separacion angular (distorsion perspectiva del parche). Se adopta el modelo
fenomenologico declarado:
      sd(da) = sd0 * [ 1 + (da/da_c)^2 ]
con da_c calibrable por el usuario (valor por defecto 20 deg, coherente con
el rango util reportado para descriptores tipo DSP-SIFT). La competencia
entre (2) -que mejora con da- y este termino -que empeora con da- produce
un optimo interior que el script localiza numericamente.

Salida: Figuras/S12_planificacion_captura.png (300 dpi) y resumen en consola.
Backend Agg: ejecutable en Google Colab sin display.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

os.makedirs('Figuras', exist_ok=True)
rng = np.random.default_rng(12)

# ----------------------------------------------------------------------
# 1. Configuracion nominal declarada (APS-C 24 MP + objetivo de 35 mm)
# ----------------------------------------------------------------------
W_SENS = 22.3          # ancho del sensor [mm]
NPX_H = 6000           # pixeles en la dimension ancha
P_PIX = W_SENS / NPX_H  # pitch de pixel [mm]
F_MM = 35.0            # distancia focal [mm]
Z_NOM = 600.0          # distancia de trabajo nominal objeto-camara [mm]
SD0 = 0.30             # ruido de correspondencia a angulo nulo [px]
DA_C = 20.0            # angulo caracteristico de degradacion del matching [deg]
N_ANILLOS = 2          # anillos de captura (0 deg y +30 deg)

F_PX = F_MM / P_PIX                     # focal en pixeles
FOV = 2.0 * np.arctan(W_SENS / (2.0 * F_MM))   # campo de vision horizontal [rad]


def gsd(z, f_mm=F_MM, p=P_PIX):
    """Ground sampling distance [mm/px] a distancia z [mm]."""
    return p * z / f_mm


def base_orbital(z, da_deg):
    """Base estereo efectiva [mm] entre dos estaciones separadas da_deg."""
    return 2.0 * z * np.sin(np.deg2rad(da_deg) / 2.0)


def sigma_d(da_deg, sd0=SD0, da_c=DA_C):
    """Ruido de correspondencia [px] en funcion de la separacion angular."""
    return sd0 * (1.0 + (da_deg / da_c) ** 2)


def sigma_z(da_deg, z=Z_NOM, degradacion=True):
    """Incertidumbre estandar de profundidad [mm] de un par convergente."""
    b = base_orbital(z, da_deg)
    sd = sigma_d(da_deg) if degradacion else SD0
    return (z ** 2) * sd / (F_PX * b)


def solapamiento(da_deg, z=Z_NOM):
    """Fraccion de solapamiento entre estaciones consecutivas (0-1)."""
    w = 2.0 * z * np.tan(FOV / 2.0)      # huella lateral del encuadre [mm]
    return 1.0 - base_orbital(z, da_deg) / w


def n_imagenes(da_deg, n_anillos=N_ANILLOS):
    return np.ceil(360.0 / da_deg) * n_anillos


# ----------------------------------------------------------------------
# 2. Optimo de separacion angular
# ----------------------------------------------------------------------
res = minimize_scalar(sigma_z, bounds=(2.0, 80.0), method='bounded')
DA_OPT = float(res.x)
SZ_OPT = float(res.fun)

GSD_NOM = gsd(Z_NOM)
SZ_15 = float(sigma_z(15.0))
SZ_15_ideal = float(sigma_z(15.0, degradacion=False))
OVL_OPT = float(solapamiento(DA_OPT))
N_OPT = float(n_imagenes(DA_OPT))

# separacion angular maxima admisible por criterio de solapamiento >= 0.60
da_grid = np.linspace(2.0, 80.0, 4000)
ovl = solapamiento(da_grid)
DA_OVL60 = float(da_grid[np.argmin(np.abs(ovl - 0.60))])
DA_OVL80 = float(da_grid[np.argmin(np.abs(ovl - 0.80))])

# relacion de aspecto del presupuesto: profundidad frente a lateral
RATIO = SZ_OPT / GSD_NOM

print('=' * 68)
print('PRESUPUESTO OPTICO DE LA SESION FOTOGRAMETRICA')
print('=' * 68)
print('Pitch de pixel                p     = {:.4f} mm'.format(P_PIX))
print('Focal en pixeles              f_px  = {:.0f} px'.format(F_PX))
print('Campo de vision horizontal    FOV   = {:.1f} deg'
      .format(np.rad2deg(FOV)))
print('GSD a Z = {:.0f} mm             GSD   = {:.4f} mm/px ({:.1f} um)'
      .format(Z_NOM, GSD_NOM, GSD_NOM * 1000))
print('-' * 68)
print('Separacion angular optima     da*   = {:.1f} deg'.format(DA_OPT))
print('Incertidumbre de profundidad  sZ*   = {:.4f} mm ({:.1f} um)'
      .format(SZ_OPT, SZ_OPT * 1000))
print('Relacion sZ*/GSD                    = {:.2f}'.format(RATIO))
print('Solapamiento en da*                 = {:.3f}'.format(OVL_OPT))
print('Imagenes requeridas en da*    N     = {:.0f} ({} anillos)'
      .format(N_OPT, N_ANILLOS))
print('-' * 68)
print('sZ a da = 15 deg (con degradacion)  = {:.4f} mm'.format(SZ_15))
print('sZ a da = 15 deg (sin degradacion)  = {:.4f} mm'.format(SZ_15_ideal))
print('da con solapamiento = 0.80          = {:.1f} deg'.format(DA_OVL80))
print('da con solapamiento = 0.60          = {:.1f} deg'.format(DA_OVL60))
print('-' * 68)
# Optimo restringido: el solapamiento minimo del pipeline SfM domina el diseno
SZ_OVL60 = float(sigma_z(DA_OVL60))
SZ_OVL80 = float(sigma_z(DA_OVL80))
PENAL60 = 100.0 * (SZ_OVL60 / SZ_OPT - 1.0)
PENAL80 = 100.0 * (SZ_OVL80 / SZ_OPT - 1.0)
print('OPTIMO RESTRINGIDO (la restriccion activa es el solapamiento)')
print('  o >= 0.60 -> da = {:.1f} deg, sZ = {:.1f} um, N = {:.0f}, '
      'penalizacion = {:.1f} por ciento'
      .format(DA_OVL60, SZ_OVL60 * 1000, n_imagenes(DA_OVL60), PENAL60))
print('  o >= 0.80 -> da = {:.1f} deg, sZ = {:.1f} um, N = {:.0f}, '
      'penalizacion = {:.1f} por ciento'
      .format(DA_OVL80, SZ_OVL80 * 1000, n_imagenes(DA_OVL80), PENAL80))
print('=' * 68)

# ----------------------------------------------------------------------
# 3. Figura de tres paneles
# ----------------------------------------------------------------------
AZUL = '#003366'
ROJO = '#B22222'
VERDE = '#2E7D32'
NARANJA = '#E07B00'

fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.4))

# --- (a) GSD frente a distancia de trabajo -----------------------------
z_grid = np.linspace(200.0, 1500.0, 500)
for f_i, est in zip([24.0, 35.0, 50.0], ['-.', '-', '--']):
    ax[0].plot(z_grid, gsd(z_grid, f_mm=f_i) * 1000.0, est, lw=1.8,
               color=AZUL if f_i == 35.0 else 'gray',
               label='f = {:.0f} mm'.format(f_i))
ax[0].axhline(100.0, color=ROJO, lw=1.2, ls=':')
ax[0].text(215, 106, 'criterio GSD = 100 um', color=ROJO, fontsize=8)
ax[0].plot([Z_NOM], [GSD_NOM * 1000.0], 'o', ms=7, color=ROJO, zorder=5)
ax[0].annotate('config. nominal\n{:.0f} um'.format(GSD_NOM * 1000.0),
               xy=(Z_NOM, GSD_NOM * 1000.0), xytext=(760, 40),
               fontsize=8, color=ROJO,
               arrowprops=dict(arrowstyle='->', color=ROJO, lw=0.9))
ax[0].set_xlabel('Distancia de trabajo $Z$ [mm]')
ax[0].set_ylabel('GSD [$\\mu$m/px]')
ax[0].set_title('(a) Resolucion lateral', fontsize=10, loc='left')
ax[0].legend(fontsize=8, frameon=False)
ax[0].grid(alpha=0.3)

# --- (b) Incertidumbre de profundidad frente a separacion angular ------
sz_deg = sigma_z(da_grid) * 1000.0
sz_ide = sigma_z(da_grid, degradacion=False) * 1000.0
ax[1].plot(da_grid, sz_ide, '--', color='gray', lw=1.5,
           label='geometria pura $\\sigma_d=\\sigma_{d0}$')
ax[1].plot(da_grid, sz_deg, '-', color=AZUL, lw=2.0,
           label='con degradacion de matching')
ax[1].plot([DA_OPT], [SZ_OPT * 1000.0], 'o', ms=8, color=ROJO, zorder=5)
ax[1].annotate('$\\Delta\\alpha^{{*}}$ = {:.1f}$^\\circ$\n$\\sigma_Z$ = {:.0f} $\\mu$m'
               .format(DA_OPT, SZ_OPT * 1000.0),
               xy=(DA_OPT, SZ_OPT * 1000.0), xytext=(DA_OPT + 12, SZ_OPT * 1000 + 95),
               fontsize=8.5, color=ROJO,
               arrowprops=dict(arrowstyle='->', color=ROJO, lw=0.9))
ax[1].axvspan(DA_OVL60, 80.0, color=NARANJA, alpha=0.12)
Y1MAX = float(np.nanmax(sz_deg[da_grid < 70]))
ax[1].text(DA_OVL60 + 1.0, 0.58 * Y1MAX,
           'solapamiento < 0.60\n(SfM inestable)', fontsize=7.5,
           color=NARANJA)
ax[1].set_xlim(2, 70)
ax[1].set_ylim(0, Y1MAX)
ax[1].set_xlabel('Separacion angular $\\Delta\\alpha$ [$^\\circ$]')
ax[1].set_ylabel('$\\sigma_Z$ [$\\mu$m]')
ax[1].set_title('(b) Presupuesto de profundidad', fontsize=10, loc='left')
ax[1].legend(fontsize=7.5, frameon=False, loc='upper right')
ax[1].grid(alpha=0.3)

# --- (c) Solapamiento y numero de imagenes -----------------------------
ax[2].plot(da_grid, ovl, '-', color=VERDE, lw=2.0)
ax[2].axhline(0.60, color=ROJO, lw=1.1, ls=':')
ax[2].axhline(0.80, color=NARANJA, lw=1.1, ls=':')
ax[2].text(38, 0.625, 'minimo admisible 0.60', color=ROJO, fontsize=7.5)
ax[2].text(38, 0.825, 'recomendado 0.80', color=NARANJA, fontsize=7.5)
ax[2].axvline(DA_OPT, color=AZUL, lw=1.2, ls='--')
ax[2].set_xlim(2, 70)
ax[2].set_ylim(0, 1.0)
ax[2].set_xlabel('Separacion angular $\\Delta\\alpha$ [$^\\circ$]')
ax[2].set_ylabel('Solapamiento $o$ [-]', color=VERDE)
ax[2].tick_params(axis='y', labelcolor=VERDE)
ax[2].grid(alpha=0.3)

ax2b = ax[2].twinx()
ax2b.plot(da_grid, n_imagenes(da_grid), '-', color=AZUL, lw=1.6, alpha=0.85)
ax2b.set_ylabel('$N$ imagenes ({} anillos)'.format(N_ANILLOS), color=AZUL)
ax2b.tick_params(axis='y', labelcolor=AZUL)
ax2b.set_yscale('log')
ax[2].set_title('(c) Solapamiento y volumen de datos', fontsize=10, loc='left')

fig.tight_layout()
fig.savefig('Figuras/S12_planificacion_captura.png', dpi=300,
            bbox_inches='tight')
# plt.show()
print('Figura guardada: Figuras/S12_planificacion_captura.png')
