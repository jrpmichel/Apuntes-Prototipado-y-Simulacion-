# ============================================================================
# Semana11_metrica_area.py
# Curso: Prototipado y Simulacion (DITA) -- Universidad La Salle Bajio
# Semana 11: Validacion Experimental y Metrologia Computacional
#
# Metrica de area de Oberkampf & Roy evaluada sobre funciones de
# distribucion acumulada (CDF): ensamble Monte Carlo del modelo contra
# replicas experimentales. Caso: carga hidraulica H de una bomba
# centrifuga en su punto de operacion nominal (familia de proyectos:
# banco de pruebas de bombas centrifugas).
#
# Modelo:      H(Q) = a - b Q^2, con a y b inciertos (hoja de datos +
#              ajuste de curva caracteristica).
# Experimento: n replicas de H medidas con transductor de presion;
#              el sistema real presenta perdidas adicionales no
#              modeladas (deficiencia de modelo).
#
# Salida: Figuras/S11_metrica_cdf.png
# Reproducibilidad: numpy.random.default_rng(semilla fija)
# ============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('Figuras', exist_ok=True)
rng = np.random.default_rng(1164)

# ----------------------------------------------------------------------------
# 1. ENSAMBLE MONTE CARLO DEL MODELO (propagacion de u_input)
# ----------------------------------------------------------------------------
Q_op = 0.090          # [m3/s] caudal nominal de operacion
M = 20000             # muestras Monte Carlo

a_mu, a_u = 32.0, 0.60    # [m]        ordenada de la curva caracteristica
b_mu, b_u = 1250.0, 80.0  # [m s2/m6]  coeficiente cuadratico

a_s = rng.normal(a_mu, a_u, M)
b_s = rng.normal(b_mu, b_u, M)
H_sim = a_s - b_s * Q_op**2

# ----------------------------------------------------------------------------
# 2. REPLICAS EXPERIMENTALES (fisica real con perdidas no modeladas)
# ----------------------------------------------------------------------------
n_exp = 40
H_true_mu = 21.35     # [m] media real: perdidas adicionales de succion
s_exp = 0.55          # [m] repetibilidad del transductor + proceso
H_exp = rng.normal(H_true_mu, s_exp, n_exp)

# ----------------------------------------------------------------------------
# 3. CDFs EMPIRICAS Y METRICA DE AREA
# ----------------------------------------------------------------------------
def ecdf(sample, grid):
    """CDF empirica evaluada sobre 'grid'."""
    sample_sorted = np.sort(sample)
    return np.searchsorted(sample_sorted, grid, side='right') / len(sample)

y_lo = min(H_sim.min(), H_exp.min()) - 0.5
y_hi = max(H_sim.max(), H_exp.max()) + 0.5
y_grid = np.linspace(y_lo, y_hi, 6000)

F_S = ecdf(H_sim, y_grid)
F_D = ecdf(H_exp, y_grid)

d_area = np.trapezoid(np.abs(F_S - F_D), y_grid)   # [m], unidades de H

# Estadisticos de referencia para el reporte
E_medias = H_sim.mean() - H_exp.mean()
u_input_H = H_sim.std(ddof=1)              # dispersion inducida por entradas
u_D_H = H_exp.std(ddof=1) / np.sqrt(n_exp) # incertidumbre de la media exp.

print('=' * 72)
print('  METRICA DE AREA (OBERKAMPF & ROY) -- CARGA DE BOMBA CENTRIFUGA')
print('=' * 72)
print('Modelo (MC, M=%d):  media = %.3f m;  desv = %.3f m'
      % (M, H_sim.mean(), u_input_H))
print('Experimento (n=%d): media = %.3f m;  desv = %.3f m'
      % (n_exp, H_exp.mean(), H_exp.std(ddof=1)))
print('Error de comparacion de medias E = %.3f m' % E_medias)
print('Metrica de area d_area = %.3f m' % d_area)
print('Fraccion de d_area respecto a la media experimental: %.2f por ciento'
      % (100.0 * d_area / H_exp.mean()))

# ----------------------------------------------------------------------------
# 4. FIGURA: CDFs Y AREA ENCERRADA
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.0, 5.8), dpi=300)
ax.fill_betweenx(F_S, y_grid, y_grid, color='none')  # noop para orden de capas
ax.fill_between(y_grid, F_S, F_D, color='gold', alpha=0.55,
                label=r'Area $= d_{area} = %.2f$ m' % d_area)
ax.plot(y_grid, F_S, color='tab:blue', lw=2,
        label=r'$F_S(y)$: ensamble Monte Carlo del modelo (M=%d)' % M)
ax.step(np.sort(H_exp), np.arange(1, n_exp + 1) / n_exp, where='post',
        color='k', lw=1.6,
        label=r'$F_D(y)$: CDF empirica experimental (n=%d)' % n_exp)
ax.set_xlabel('Carga hidraulica H [m]', fontsize=12)
ax.set_ylabel('Probabilidad acumulada', fontsize=12)
ax.set_title('Metrica de area sobre CDFs: desajuste estocastico '
             'modelo-experimento', fontsize=12.5)
ax.set_ylim(-0.02, 1.02)
ax.grid(True, ls='--', alpha=0.45)
ax.legend(loc='upper left', fontsize=10, frameon=True)
fig.tight_layout()
fig.savefig('Figuras/S11_metrica_cdf.png')
plt.close(fig)
# plt.show()
print('\nFigura escrita en Figuras/S11_metrica_cdf.png')
