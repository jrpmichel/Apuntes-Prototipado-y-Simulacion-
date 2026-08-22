/* main.js — utilidades generales del sitio */
(function () {
  // Reemplaza cualquier <img> rota (imagen aún no subida) por un
  // recuadro de "pendiente" prolijo en vez del ícono roto del navegador.
  function handleBrokenImages() {
    document.querySelectorAll("img[data-fallback-name]").forEach((img) => {
      img.addEventListener("error", function onErr() {
        img.removeEventListener("error", onErr);
        const box = document.createElement("div");
        box.className = "img-pendiente";
        box.textContent = "Imagen pendiente: " + img.dataset.fallbackName;
        img.replaceWith(box);
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", handleBrokenImages);
  } else {
    handleBrokenImages();
  }
})();
