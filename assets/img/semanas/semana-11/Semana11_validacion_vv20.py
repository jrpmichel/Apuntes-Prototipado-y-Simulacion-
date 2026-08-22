# ============================================================================
# Semana11_validacion_vv20.py
# Curso: Prototipado y Simulacion (DITA) -- Universidad La Salle Bajio
# Semana 11: Validacion Experimental y Metrologia Computacional
#
# Framework ASME V&V 20 aplicado a un modelo termico transitorio de
# calentamiento de un conductor (familia de proyectos: cable OPGW bajo
# corriente de falla). El modelo computacional S omite la perdida
# radiativa (deficiencia de modelo delta_model deliberada); el
# "experimento" D se sintetiza integrando la fisica completa y
# contaminando con sesgo sistematico y ruido aleatorio de termopar.
#
# Componentes calculados:
#   u_num   : GCI por extrapolacion de Richardson (3 pasos de tiempo)
#   u_input : propagacion GUM de primer orden de u_h (coef. conveccion)
#   u_D     : combinacion de sesgo sistematico y repetibilidad (n replicas)
#   u_val   : combinacion en cuadratura (ASME V&V 20)
#
# Reproducibilidad: numpy.random.default_rng(semilla fija)
# Salidas: Figuras/S11_validacion_continua.png
#          Figuras/S11_gci_convergencia.png
# ============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

os.makedirs('Figuras', exist_ok=True)
rng = np.random.default_rng(2026)

# ----------------------------------------------------------------------------
# 1. PARAMETROS FISICOS DEL SISTEMA
# ----------------------------------------------------------------------------
P_joule = 25.0     # [W]      potencia disipada por efecto Joule
h_conv  = 12.0     # [W/m2K]  coeficiente de conveccion (parametro incierto)
u_h     = 1.0      # [W/m2K]  incertidumbre estandar de h (correlacion empirica)
A_sup   = 0.05     # [m2]     area superficial del tramo de conductor
mc      = 180.0    # [J/K]    capacitancia termica (m*c)
T_amb   = 300.0    # [K]      temperatura ambiente
eps_rad = 0.8      # [-]      emisividad (fisica omitida por el modelo S)
sigma_SB = 5.670e-8  # [W/m2K4] constante de Stefan-Boltzmann

t_final = 900.0    # [s]      horizonte de comparacion
n_pts   = 15       # puntos de comparacion experimento-modelo
t_comp  = np.linspace(60.0, t_final, n_pts)

# ----------------------------------------------------------------------------
# 2. MODELO COMPUTACIONAL S (sin radiacion) RESUELTO POR EULER EXPLICITO
# ----------------------------------------------------------------------------
def solve_model_euler(dt, h_val, t_eval):
    """Integra dT/dt = (P - h A (T - Tamb)) / mc con Euler explicito y
    devuelve T interpolada en t_eval. La solucion numerica (no la
    analitica) es la que se verifica: su error de discretizacion es u_num."""
    n_steps = int(np.ceil(t_final / dt)) + 1
    t = np.linspace(0.0, dt * (n_steps - 1), n_steps)
    T = np.empty(n_steps)
    T[0] = T_amb
    for k in range(n_steps - 1):
        dTdt = (P_joule - h_val * A_sup * (T[k] - T_amb)) / mc
        T[k + 1] = T[k] + dt * dTdt
    return np.interp(t_eval, t, T)

# ----------------------------------------------------------------------------
# 3. VERIFICACION: GCI CON TRES PASOS DE TIEMPO (Richardson)
# ----------------------------------------------------------------------------
r_ref = 2.0                       # razon de refinamiento
dt3, dt2, dt1 = 8.0, 4.0, 2.0     # grueso, medio, fino
t_ref = 600.0                     # instante de referencia para el GCI
S3 = solve_model_euler(dt3, h_conv, np.array([t_ref]))[0]
S2 = solve_model_euler(dt2, h_conv, np.array([t_ref]))[0]
S1 = solve_model_euler(dt1, h_conv, np.array([t_ref]))[0]

p_obs = np.log(abs(S3 - S2) / abs(S2 - S1)) / np.log(r_ref)
Fs = 1.25                          # factor de seguridad (Roache)
GCI_fine = Fs * abs(S1 - S2) / (r_ref**p_obs - 1.0)   # [K] expandido
u_num_scalar = GCI_fine / 1.15     # GCI ~ intervalo uniforme -> estandar
u_num = u_num_scalar * np.ones_like(t_comp)

# Solucion analitica exacta del modelo (solo para la grafica de convergencia)
tau = mc / (h_conv * A_sup)
T_exact_ref = T_amb + (P_joule / (h_conv * A_sup)) * (1.0 - np.exp(-t_ref / tau))

# ----------------------------------------------------------------------------
# 4. PREDICCION S EN LA MALLA FINA Y SENSIBILIDAD u_input
# ----------------------------------------------------------------------------
S_model = solve_model_euler(dt1, h_conv, t_comp)

# Sensibilidad local dS/dh por diferencias centradas (recorriendo el solver)
delta_h = 0.05 * h_conv
S_hp = solve_model_euler(dt1, h_conv + delta_h, t_comp)
S_hm = solve_model_euler(dt1, h_conv - delta_h, t_comp)
dS_dh = (S_hp - S_hm) / (2.0 * delta_h)
u_input = np.abs(dS_dh) * u_h

# ----------------------------------------------------------------------------
# 5. EXPERIMENTO SINTETICO D (fisica completa + metrologia de termopar)
# ----------------------------------------------------------------------------
def rhs_true(t, T):
    q_conv = h_conv * A_sup * (T[0] - T_amb)
    q_rad = eps_rad * sigma_SB * A_sup * (T[0]**4 - T_amb**4)
    return [(P_joule - q_conv - q_rad) / mc]

sol = solve_ivp(rhs_true, [0.0, t_final], [T_amb], t_eval=t_comp,
                rtol=1e-10, atol=1e-10)
T_true = sol.y[0]

b_real = 1.2       # [K] sesgo sistematico realizado (termopar descalibrado)
u_b    = 1.0       # [K] incertidumbre estandar del sesgo (certificado)
s_rep  = 0.8       # [K] desviacion estandar de una lectura individual
n_rep  = 5         # replicas por punto de comparacion

lecturas = (T_true[:, None] + b_real
            + rng.normal(0.0, s_rep, size=(n_pts, n_rep)))
D_exp = lecturas.mean(axis=1)
s_muestral = lecturas.std(axis=1, ddof=1)
u_D = np.sqrt(u_b**2 + s_muestral**2 / n_rep)

# ----------------------------------------------------------------------------
# 6. ASME V&V 20: ERROR DE COMPARACION Y BANDA DE VALIDACION
# ----------------------------------------------------------------------------
E = S_model - D_exp
u_val = np.sqrt(u_num**2 + u_input**2 + u_D**2)
k_cov = 2.0
is_valid = np.abs(E) <= k_cov * u_val

# ----------------------------------------------------------------------------
# 7. REPORTE DE TERMINAL
# ----------------------------------------------------------------------------
print('=' * 78)
print('  REPORTE DE VALIDACION ASME V&V 20 -- MODELO TERMICO TRANSITORIO')
print('=' * 78)
print('Verificacion (GCI, Euler explicito, t_ref = %.0f s):' % t_ref)
print('  S(dt=%.0f)=%.4f K  S(dt=%.0f)=%.4f K  S(dt=%.0f)=%.4f K'
      % (dt3, S3, dt2, S2, dt1, S1))
print('  Orden observado p = %.3f (teorico Euler: 1)' % p_obs)
print('  GCI_fino = %.4f K  ->  u_num = %.4f K' % (GCI_fine, u_num_scalar))
print('  |S_fino - T_exacta| = %.4f K (verificacion contra solucion cerrada)'
      % abs(S1 - T_exact_ref))
print('-' * 78)
print('Sensibilidad de entrada: max|dS/dh| = %.3f K/(W/m2K); '
      'max u_input = %.3f K' % (np.max(np.abs(dS_dh)), np.max(u_input)))
print('Metrologia: u_b = %.2f K; s promedio = %.2f K; '
      'u_D en [%.2f, %.2f] K' % (u_b, s_muestral.mean(),
                                 u_D.min(), u_D.max()))
print('-' * 78)
hdr = '%-8s%-10s%-10s%-10s%-9s%-22s' % ('t [s]', 'S [K]', 'D [K]',
                                        'E [K]', '2u_val', 'Diagnostico')
print(hdr)
print('-' * 78)
for i in range(n_pts):
    status = 'no refutado' if is_valid[i] else 'DEFICIENCIA DE MODELO'
    print('%-8.0f%-10.2f%-10.2f%-+10.2f%-9.2f%-22s'
          % (t_comp[i], S_model[i], D_exp[i], E[i],
             k_cov * u_val[i], status))
n_fail = int(np.sum(~is_valid))
print('-' * 78)
print('Puntos con |E| > 2 u_val: %d de %d '
      '(primera deteccion en t = %.0f s)'
      % (n_fail, n_pts, t_comp[~is_valid][0] if n_fail else -1))

# ----------------------------------------------------------------------------
# 8. FIGURA 1: VALIDACION CONTINUA
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.5, 6.0), dpi=300)
ax.errorbar(t_comp, D_exp, yerr=k_cov * u_D, fmt='o', color='k',
            mfc='white', ecolor='gray', elinewidth=1.3, capsize=4,
            label=r'Experimento $D \pm 2u_D$ (n=%d replicas)' % n_rep, zorder=4)
ax.plot(t_comp, S_model, 'b-', lw=2, label=r'Modelo computacional $S$')
ax.fill_between(t_comp, S_model - k_cov * u_val, S_model + k_cov * u_val,
                color='tab:blue', alpha=0.15,
                label=r'Banda de validacion $S \pm 2u_{val}$')
ax.scatter(t_comp[~is_valid], D_exp[~is_valid], color='red', s=70,
           zorder=5, label=r'$|E| > 2u_{val}$: deficiencia de modelo')
ax.set_xlabel('Tiempo [s]', fontsize=12)
ax.set_ylabel('Temperatura [K]', fontsize=12)
ax.set_title('Validacion ASME V&V 20: modelo convectivo vs. sistema '
             'convectivo-radiativo', fontsize=12.5)
ax.grid(True, ls='--', alpha=0.45)
ax.legend(loc='lower right', fontsize=10, frameon=True)
fig.tight_layout()
fig.savefig('Figuras/S11_validacion_continua.png')
plt.close(fig)

# ----------------------------------------------------------------------------
# 9. FIGURA 2: ESTUDIO DE CONVERGENCIA (GCI)
# ----------------------------------------------------------------------------
dts = np.array([16.0, 8.0, 4.0, 2.0, 1.0, 0.5])
errs = np.array([abs(solve_model_euler(d, h_conv,
                 np.array([t_ref]))[0] - T_exact_ref) for d in dts])

fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=300)
ax.loglog(dts, errs, 'o-', color='tab:blue', lw=1.8,
          label='Error de discretizacion (Euler)')
ref = errs[0] * (dts / dts[0])**1.0
ax.loglog(dts, ref, 'k--', lw=1.2, label=r'Pendiente de referencia $p=1$')
ax.axhline(u_num_scalar, color='tab:red', ls=':', lw=1.5,
           label=r'$u_{num}$ estimada por GCI')
ax.set_xlabel(r'Paso de tiempo $\Delta t$ [s]', fontsize=12)
ax.set_ylabel(r'$|T_{\Delta t} - T_{exacta}|$ en $t = 600$ s [K]', fontsize=12)
ax.set_title('Verificacion: convergencia de malla temporal y GCI',
             fontsize=12.5)
ax.grid(True, which='both', ls='--', alpha=0.45)
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig('Figuras/S11_gci_convergencia.png')
plt.close(fig)
# plt.show()
print('\nFiguras escritas en Figuras/.')
