/* ============================================================
   pyodide-runner.js
   Ejecuta las celdas de código Python directamente en el navegador
   del alumno, sin servidor ni backend, usando Pyodide (Python
   compilado a WebAssembly). Pensado para los scripts numpy /
   scipy / matplotlib de los Apuntes de Prototipado y Simulación.

   Limitaciones honestas (documentar también en el README):
   - No hay acceso a ANSYS Workbench, Meshroom ni a ningún software
     de escritorio: esas partes del curso siguen siendo capturas
     de pantalla / video, no código ejecutable.
   - gmsh y JAX tampoco existen como paquetes de Pyodide (WebAssembly):
     las celdas que los requieren no tienen botón de ejecución aquí,
     aunque el código en sí ya se verificó correcto por fuera del
     navegador (ver README, sección de verificación de código).
   - Un script que lee archivos externos (imágenes, CSV, .stl) no
     funcionará aquí a menos que ese archivo se cargue primero al
     sistema de archivos virtual de Pyodide.
   - La primera ejecución en cada página tarda unos segundos: debe
     descargar el runtime de Python (WASM) y los paquetes.
   ============================================================ */

(function () {
  const PYODIDE_CDN = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
  const PACKAGES = ["numpy", "scipy", "matplotlib"];

  // Paquetes adicionales que SI existen precompilados para Pyodide, cargados
  // solo si la celda los importa (para no alargar la primera descarga en
  // celdas que no los necesitan). El nombre de la izquierda es el que
  // aparece en el "import X" de Python; el de la derecha es el nombre del
  // paquete en el indice de Pyodide.
  const OPTIONAL_PACKAGES = {
    pandas: "pandas",
    networkx: "networkx",
    sympy: "sympy",
    skimage: "scikit-image",
  };
  const loadedPackages = new Set(PACKAGES);

  function detectOptionalPackages(code) {
    const found = [];
    for (const modName of Object.keys(OPTIONAL_PACKAGES)) {
      const re = new RegExp("^\\s*(?:import|from)\\s+" + modName + "\\b", "m");
      if (re.test(code) && !loadedPackages.has(OPTIONAL_PACKAGES[modName])) {
        found.push(OPTIONAL_PACKAGES[modName]);
      }
    }
    return found;
  }

  let pyodideReadyPromise = null;
  let pyodideScriptPromise = null;

  function loadPyodideScript() {
    if (pyodideScriptPromise) return pyodideScriptPromise;
    pyodideScriptPromise = new Promise((resolve, reject) => {
      const s = document.createElement("script");
      s.src = PYODIDE_CDN + "pyodide.js";
      s.onload = resolve;
      s.onerror = () => reject(new Error("No se pudo cargar Pyodide desde el CDN."));
      document.head.appendChild(s);
    });
    return pyodideScriptPromise;
  }

  async function getPyodide(onStatus) {
    if (pyodideReadyPromise) return pyodideReadyPromise;
    pyodideReadyPromise = (async () => {
      onStatus("Descargando runtime de Python (WebAssembly)…");
      await loadPyodideScript();
      const pyodide = await window.loadPyodide({ indexURL: PYODIDE_CDN });
      onStatus("Instalando numpy, scipy y matplotlib…");
      await pyodide.loadPackage(PACKAGES);
      // Backend no interactivo + captura de figuras a memoria.
      await pyodide.runPythonAsync(`
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as _plt_bootstrap
`);
      onStatus("Listo.");
      return pyodide;
    })();
    return pyodideReadyPromise;
  }

  function b64ImagesFromMatplotlib(pyodide) {
    // Devuelve una lista de imágenes PNG en base64, una por figura abierta,
    // y cierra las figuras para no arrastrarlas a la siguiente ejecución.
    const code = `
import base64, io
import matplotlib.pyplot as plt
_imgs = []
for _n in plt.get_fignums():
    _fig = plt.figure(_n)
    _buf = io.BytesIO()
    _fig.savefig(_buf, format="png", dpi=130, bbox_inches="tight")
    _buf.seek(0)
    _imgs.append(base64.b64encode(_buf.read()).decode("ascii"))
plt.close("all")
_imgs
`;
    const result = pyodide.runPython(code);
    return result.toJs();
  }

  async function runCell(cell) {
    const btn = cell.querySelector(".run-btn");
    const textarea = cell.querySelector(".code-input");
    const outputWrap = cell.querySelector(".code-output-wrap");
    const consoleEl = cell.querySelector(".code-output-console");
    const figuresEl = cell.querySelector(".code-output-figures");
    const statusEl = cell.querySelector(".code-status");

    btn.disabled = true;
    outputWrap.hidden = false;
    consoleEl.classList.remove("has-error");
    consoleEl.textContent = "";
    figuresEl.innerHTML = "";
    const originalLabel = btn.textContent;
    btn.textContent = "Cargando…";

    try {
      const pyodide = await getPyodide((msg) => {
        statusEl.textContent = msg;
        btn.textContent = "Cargando…";
      });

      const code = textarea.value;
      const extras = detectOptionalPackages(code);
      if (extras.length) {
        btn.textContent = "Cargando…";
        statusEl.textContent = "Instalando " + extras.join(", ") + "…";
        await pyodide.loadPackage(extras);
        extras.forEach((p) => loadedPackages.add(p));
      }

      btn.textContent = "Ejecutando…";
      statusEl.textContent = "Ejecutando en tu navegador (WebAssembly, sin servidor)…";

      let stdout = "";
      pyodide.setStdout({ batched: (s) => { stdout += s + "\n"; } });
      pyodide.setStderr({ batched: (s) => { stdout += s + "\n"; } });

      await pyodide.runPythonAsync(code);

      const imgs = b64ImagesFromMatplotlib(pyodide);
      imgs.forEach((b64) => {
        const img = document.createElement("img");
        img.src = "data:image/png;base64," + b64;
        img.alt = "Figura generada por el script";
        figuresEl.appendChild(img);
      });

      consoleEl.textContent = stdout.trim() || "(el script se ejecutó sin imprimir nada en consola)";
      statusEl.textContent = "Ejecución terminada.";
    } catch (err) {
      consoleEl.classList.add("has-error");
      consoleEl.textContent = String(err && err.message ? err.message : err);
      statusEl.textContent = "La ejecución terminó con un error — revisa el código o si depende de archivos externos.";
    } finally {
      btn.disabled = false;
      btn.textContent = originalLabel;
    }
  }

  function initCodeCells() {
    document.querySelectorAll(".code-cell").forEach((cell) => {
      const btn = cell.querySelector(".run-btn");
      if (!btn) return;
      btn.addEventListener("click", () => runCell(cell));
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initCodeCells);
  } else {
    initCodeCells();
  }
})();
