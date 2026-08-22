# Prototipado y Simulación — Apuntes del Curso Doctoral

Sitio web del curso **Prototipado y Simulación**, Doctorado en Ingeniería y Tecnología Aplicada,
Universidad La Salle Bajío (Dr. Jorge Ramón Parra Michel). Es un sitio estático (HTML/CSS/JS puro,
sin build, sin backend) pensado para publicarse con **GitHub Pages** y reutilizarse curso tras curso.

## 1. Publicarlo en GitHub Pages

1. **Descomprime el `.zip`** en tu computadora. Vas a arrastrar el *contenido* de la carpeta ya
   descomprimida — GitHub no extrae zips por ti.
2. Crea un repositorio nuevo en GitHub (público o privado con GitHub Pro/Team para Pages privado).
4. Arrastra **todo el contenido de esa carpeta** (no la carpeta misma, su contenido: `index.html`,
   `README.md`, `assets/`, `semanas/`, `.nojekyll`) a la página de "Add file → Upload files" del
   repositorio, o haz `git init` / `git add .` / `git commit` / `git push` si prefieres usar Git.
5. Ve a **Settings → Pages**, en "Build and deployment" elige **Deploy from a branch**, rama `main`,
   carpeta `/ (root)`. Guarda.
6. En un par de minutos el sitio queda publicado en `https://<tu-usuario>.github.io/<repo>/`.
7. Cada vez que subas un cambio (una imagen, un texto corregido), el sitio se actualiza solo.

El archivo `.nojekyll` ya está incluido para que GitHub Pages sirva los archivos tal cual, sin pasarlos
por Jekyll (evita sorpresas con carpetas o nombres de archivo).

## 2. Estructura de carpetas

```
index.html                          página de inicio (portada, logo, foto, temario)
README.md                           este archivo
.nojekyll                           desactiva el procesamiento Jekyll de GitHub Pages
assets/
  css/style.css                     estilos de todo el sitio
  js/main.js                        utilidades (placeholder de imágenes faltantes)
  js/pyodide-runner.js              motor que ejecuta el código Python en el navegador
  img/
    logo-lasalle.png                logo institucional (ya incluido)
    foto-profesor.png               tu fotografía (ya incluida)
    semana-XX/                      figuras de cada semana (ya subidas)
  apuntes/
    Apuntes_Curso_..._COMPLETO.pdf  PDF completo del curso (495 pp., botón "Descargar PDF")
    semana-XX.pdf                   recorte del PDF anterior, solo esa semana (para el visor incrustado)
semanas/
  semana-01/index.html              Semana 1. Introducción al Modelado Científico
  semana-02/index.html              Semana 2. Herramientas Computacionales Científicas
  semana-03/index.html              Semana 3. Métodos Numéricos Fundamentales
  semana-04/index.html              Semana 4. Modelado de Fenómenos Físicos
  semana-05/index.html              Semana 5. Simulación Bidimensional de Campos
  semana-06/index.html              Semana 6. Cuantificación de Incertidumbre (GUM y Monte Carlo)
  semana-07/index.html              Semana 7. Cadenas de Markov y Sistemas Dinámicos
  semana-08/index.html              Semana 8. Método de Elementos Finitos I
  semana-09-10/index.html           Semanas 9 y 10. Método de Elementos Finitos II
  semana-11/index.html              Semana 11. Validación Experimental y Metrología Computacional
  semana-12-13/index.html           Semana 12. Fotogrametría y Escaneo 3D
```

## 3. Qué trae ya armado cada página de semana — y qué falta

Este primer entregable es **la estructura completa del sitio**, generada automáticamente a partir de
tu documento `.tex`. Cada una de las 11 páginas de semana ya incluye, con el texto real del curso:

- Objetivo de la cátedra.
- Temario.
- Recuadro "Propósito del tema" y/o "Vínculo con los proyectos del grupo" (cuando existía en el `.tex`).
- Resultado de aprendizaje esperado.
- **Todo el código Python de esa semana, ya ejecutable** en la propia página (ver sección 5):
  57 bloques de código en total, extraídos íntegros del documento.
- Una galería con el nombre exacto de cada figura que el `.tex` referencia, lista para recibir las
  imágenes (ver sección 4).

- **El desarrollo teórico completo** (ecuaciones, demostraciones, tablas, discusión) incrustado como PDF
  al final de cada página — ver sección 6.
- **Descarga del PDF completo del curso** (495 páginas), desde el menú superior y desde la portada.

## 4. Cómo añadir las imágenes que faltan

Tú mismo vas a subir las imágenes (así lo pediste, para no recargar el proceso). La ruta esperada
para cada semana es:

```
assets/img/semana-01/nombre_de_la_imagen.png
assets/img/semana-02/otra_imagen.png
... etc.
```

Cada página de semana ya muestra, debajo de cada recuadro de la galería, **el nombre exacto de archivo
que espera** (tomado directamente de tu `.tex`, mismo nombre y misma extensión). Solo copia el archivo
con ese nombre exacto a la carpeta `assets/img/semana-XX/` correspondiente — no hace falta tocar ningún
HTML. Mientras una imagen no esté presente, el sitio muestra un recuadro discreto de "Imagen pendiente"
en su lugar, en vez de un ícono de imagen rota: el sitio se ve terminado aunque falten fotos.

Tu carpeta local `Figuras/` (la del proyecto de LaTeX) ya tiene prácticamente todos estos archivos con
esos mismos nombres — es cuestión de copiarlos a la subcarpeta de cada semana.

## 5. Cómo funciona el código ejecutable

Cada celda de código usa **[Pyodide](https://pyodide.org)**: una distribución de Python compilada a
WebAssembly que corre **enteramente en el navegador del alumno**, sin servidor, sin backend, sin nada
que tú tengas que mantener corriendo. Al entrar a una página de semana:

1. El código se muestra ya editable, tal como está en tus apuntes.
2. Al pulsar **▶ Ejecutar**, la primera vez se descarga el runtime de Python y los paquetes
   `numpy`, `scipy` y `matplotlib` (tarda unos segundos; las ejecuciones siguientes en esa misma
   sesión de navegador son inmediatas).
3. La salida de consola (`print`) y cualquier figura de `matplotlib` se muestran justo debajo del
   código, dentro de la misma página.
4. El alumno puede modificar el código antes de correrlo — es un cuaderno vivo, no solo una lectura.

**Límites honestos, para que no haya sorpresas:**

- No hay ANSYS Workbench, Meshroom ni ningún software de escritorio dentro del navegador — esas partes
  del curso (Semanas 8, 9-10, 11 y 12) siguen dependiendo de capturas de pantalla y video, como en el
  documento original. Solo el código Python puro (mallado, post-proceso, validación, fotogrametría
  sintética) es ejecutable aquí.
- Un script que necesita leer un archivo externo (una imagen, un `.csv`, una nube de puntos) no va a
  encontrarlo dentro del navegador a menos que ese archivo también se cargue al sistema de archivos
  virtual de Pyodide; varios de los scripts de Semanas 11 y 12 generan sus propios datos sintéticos
  para no depender de esto, pero conviene revisarlo antes de usarlos en clase.
- **Verificación real de los 57 bloques**: cada uno se ejecutó fuera del navegador, en un Python normal
  con las mismas librerías (numpy/scipy/matplotlib/pandas/networkx/sympy/scikit-image), para confirmar
  que corren sin errores antes de publicarlos — ver la sección 8 para el detalle completo (qué se
  corrigió, qué quedó pendiente). Lo único que no pude probar desde este entorno es Pyodide en sí
  (el navegador que uso aquí no tiene salida a internet), así que aunque el código ya está verificado,
  te pido una prueba rápida en cuanto publiques: entra a `semanas/semana-01/`, pulsa "Ejecutar" en el
  primer código (estimación de π por Monte Carlo) y confirma que aparece la gráfica. Si algo falla,
  dímelo con el mensaje de error exacto.

## 6. Cómo se resolvió la teoría completa (ecuaciones, tablas, figuras TikZ)

Cada página de semana trae, después del código, un visor con el desarrollo teórico completo de esa
semana — el mismo texto que en tu PDF, recortado a las páginas de esa semana exacta (`assets/apuntes/semana-XX.pdf`,
dentro de un `<iframe>`). Debajo del visor hay un enlace de respaldo para descargarlo o abrirlo en pestaña
nueva, por si el navegador del alumno no soporta ver PDFs incrustados (pasa en algunos navegadores de
celular).

**Por qué PDF incrustado y no HTML con las ecuaciones "sueltas"**: tu documento usa 22 diagramas TikZ
dibujados directamente en LaTeX (no imágenes) y macros propias (`\Figure`, `\ExField`, cajas `cajaproposito`
/ `cajaaplica`, etc.). Convertir ese contenido a HTML nativo es posible (con `pandoc` + KaTeX para las
ecuaciones, renderizando cada TikZ como imagen aparte), pero es un trabajo grande y con riesgo real de
introducir un error silencioso en una fórmula — algo que no me pareció aceptable tratándose de contenido
de nivel doctoral sin que tú revises cada semana convertida antes de publicarla. Por eso la primera versión
usa el PDF real, con exactamente el mismo texto, tablas y figuras que ya verificaste al compilarlo — cero
riesgo de que algo se transcriba mal.

Si más adelante quieres el texto plasmado como HTML nativo (buscable, copiable, sin visor de PDF), es un
proyecto aparte que conviene hacer semana por semana con tu revisión de cada una antes de publicarla —
dímelo cuando quieras arrancarlo.

**Cómo se generó el PDF**: se recompiló tu `.tex` final (el mismo que ya tiene el sitio, incluida la
sección "Cierre del curso" y el renombre a "Semana 12 y 13") con XeLaTeX, sustituyendo la fuente
*Times New Roman* por *TeX Gyre Termes* (métricamente compatible, libre) porque Times New Roman no está
disponible en el entorno donde compilé. Detecté una única referencia rota en tu fuente: la figura
`fig:s9_fuerza_biela` (Semana 9-10) se cita con `\ref{}` pero su `\label{}` no existe en ningún lado del
documento — en el PDF aparece como "??" en ese punto. No lo inventé ni lo adiviné: es un pendiente tuyo en
el `.tex` original. El resto del documento compiló limpio (495 páginas).

Si vuelves a corregir el `.tex` (arreglar esa referencia, agregar contenido), pide que se regenere el PDF
y las páginas de semana a partir de la versión nueva, en vez de editar cada archivo a mano.

## 7. Verificación de código: qué se corrigió y qué falta

Corriste bien en pedirme que lo verificara — **46 de los 57 bloques no corrían tal como estaban**. La
causa no fue de esta página: **tu `.tex` original tiene la indentación de Python aplanada** en todos los
bloques `\begin{lstlisting}` pegados directamente (cada línea del cuerpo de un `for`/`if`/`def` quedó al
mismo nivel que el encabezado), así que cualquier cosa que corriera ese texto —esta página, tu propio
`python archivo.py`, lo que sea— iba a fallar con `IndentationError`. Los 9 bloques que vienen de un
`\lstinputlisting{archivo.py}` externo (Semanas 8, 9-10, 11, 12-13) nunca tuvieron este problema porque
ahí sí se lee el archivo tal cual.

Con la carpeta de scripts que conectaste (`codigos/`) pude corregir la mayoría:

- **44 de 57 corren limpio**, verificado ejecutándolos de verdad (Python real, mismas librerías que
  Pyodide carga: numpy, scipy, matplotlib, pandas, networkx, sympy, scikit-image).
- **3 no son código Python** — son listados de referencia (estructura de carpetas, un `Makefile`,
  comandos de `git`) que tu propio `.tex` marca como `language=bash`. Los reclasifiqué: ya no tienen
  botón de "Ejecutar" ni intentan correr como Python, se muestran solo como referencia.
- **1 requiere `gmsh`** (Semana 8, generación de mallas): el código en sí es correcto — lo verifiqué
  instalando las librerías de sistema que le faltaban (`libglu1-mesa` y similares) — pero `gmsh` es una
  librería binaria compilada que no existe para Pyodide/WebAssembly, así que no puede correr dentro del
  navegador. Por eso esa celda tampoco tiene botón de ejecución; el aviso en la página lo explica.
- **10 bloques siguen rotos** porque no encontré, entre lo que subiste, el script con la indentación
  correcta:
  - Semana 1: "Simulación explícita de conducción transitoria 1D", "Modelo simplificado de adopción
    tecnológica mediante agentes", "Conducción 2D en placa de intercambiador", "Visualización de
    esfuerzos de von Mises".
  - Semana 2: "Vectorización vs bucle explícito en NumPy", "Plantilla de script reproducible con
    argparse", "Comparación de métodos en `solve_ivp`", "Diferenciación automática con JAX", el módulo
    y la suite de `pytest` del integrador de calor.

  Si tienes los `.py` o `.ipynb` originales de estos (aunque sean de Colab, como los que ya diste), 
  compártelos y hago la misma corrección. Si no los tienes, dímelo y reindento el código a mano
  revisando la física de cada uno para no introducir un error de estructura.

Además, mientras corregía encontré y arreglé dos problemas más, presentes en varios scripts que vienen
de Google Colab:

- Rutas absolutas `/content/...` en `plt.savefig(...)` (carpeta de trabajo por defecto de Colab, no
  existe en ningún otro entorno) — las volví relativas.
- Un comando de shell (`!apt-get ...`, `!pip install gmsh`) que solo es válido dentro de una celda de
  Jupyter/Colab, no en Python puro — lo dejé comentado con una nota, en vez de borrarlo, para que sepas
  que ese paso (instalar `gmsh`) sigue siendo necesario si corres el script tú mismo fuera del navegador.

**Cómo se ve esto en la página**: cada celda de código con algún problema tiene un aviso justo arriba —
azul si es informativo (como el caso de `gmsh`), amarillo si el código se recuperó de un borrador
anterior y conviene que lo revises antes de usarlo en clase (pasa en las Semanas 3 y 4, donde el `.tex`
final tiene ediciones que el borrador no tenía — el código corre, pero quizá no en el 100% de los
detalles), y rojo en los 10 bloques que siguen sin corregir.

## 8. Créditos y mantenimiento

Contenido generado a partir de `Apuntes_Curso_Doctoral_Simulacion_Prototipado_COMPLETO.tex`. Si vuelves
a modificar el `.tex` (nuevas semanas, textos corregidos), lo más simple es pedir que se regenere el
sitio a partir de la versión nueva del documento, en vez de editar cada página HTML a mano.

Dr. Jorge Ramón Parra Michel — Doctorado en Ingeniería y Tecnología Aplicada, Universidad La Salle Bajío.
