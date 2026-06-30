# PÓRTICO-UACh — Laboratorio Virtual de Análisis Estructural 3D

Capa **académica** de la **Universidad Austral de Chile** (Instituto de Obras Civiles,
Facultad de Ciencias de la Ingeniería) construida sobre el motor open source
**[portico-core](https://github.com/jpreyes/portico-core)**.

Material docente desarrollado por **Dr. Juan Patricio Reyes C.**, para estudiantes de
**Arquitectura**, **Construcción** e **Ingeniería Civil en Obras Civiles**.
Aplicación web de **análisis y diseño estructural 3D (FEM)** que corre íntegramente en el
navegador, sin instalación.

> Herramienta de carácter **docente**: para proyectos reales use software profesional
> validado y el criterio de un ingeniero calculista competente.

## Arquitectura — overlay sobre portico-core

PÓRTICO-UACh **no forkea** el motor: lo consume como **submódulo git** y sólo aporta lo
específico de la UACh. El código del motor/UI vive en `vendor/portico-core/`; la raíz del
overlay aporta marca, datos servidos, asistente y despliegue.

```
portico-uach/
├── index.html              # del overlay: carga el core desde vendor/ + branding + bootstrap
├── branding.default.json   # marca UACh (escudo, logos, enlaces) — white-label del core
├── overlay/bootstrap.js    # registra contribuciones UACh por los seams del core
│                           #   (badge ACADÉMICO, traducciones EN del bloque legal)
├── manifest.webmanifest    # PWA UACh
├── sw.js                   # service worker (offline)
├── icons/                  # escudo/logos UACh + iconos PWA
├── assistant/              # DATOS del asistente (preset Chile: NCh) servidos en la raíz
├── examples/               # ejemplos guiados (del core) servidos en la raíz
├── worker/asistente.js     # Cloudflare Worker: API del asistente IA (BYO LLM)
├── wrangler.jsonc          # despliegue Cloudflare
└── vendor/portico-core/    # SUBMÓDULO (motor + UI), pinneado a un commit. NO se edita.
```

El core es **white-label por configuración** y expone *seams* de extensión
(`js/branding.js`, `js/ext/extensions.js`, `js/solver/backend.js`, `js/i18n/`); el overlay
los usa sin tocar el código del motor. Detalle del contrato en
[`docs/EXTENDING.md`](https://github.com/jpreyes/portico-core/blob/main/docs/EXTENDING.md)
del core.

## Cómo correr (local)

Requiere `git` (con submódulos) y Python 3.

```bash
# Clonar CON el submódulo del core:
git clone --recurse-submodules https://github.com/jpreyes/portico-uach
cd portico-uach
# Si ya clonaste sin submódulos:
git submodule update --init --recursive

# Servir como sitio estático (servidor sin caché, MIME correctos):
python serve.py 8765
```

Luego abrir **http://localhost:8765**. (La app requiere HTTP: no funciona abriendo
`index.html` como `file://`.) El asistente IA queda inactivo en local salvo que configures
un endpoint LLM (ver abajo); el resto de la app funciona completo.

> Este repositorio se llamaba `structweb3d`; GitHub redirige la URL antigua.

## Asistente IA (endpoint LLM «trae tu propio servicio»)

El core es agnóstico del proveedor: en el diálogo **Asistente** se configura la URL del
endpoint (se guarda en `localStorage`). El **Cloudflare Worker** de este overlay
(`worker/asistente.js`) implementa ese endpoint:

- `POST /api/assistant` · `{ message }` → `{ spec }` (spec según `spec.schema.json` del core)
- `POST /api/assistant/modificar` · `{ message, model, selection }` → `{ ops }`
- `POST /api/assistant/feedback` · `{ id, comentario? }`
- `GET  /api/assistant/log?token=…` (corpus para revisión)

El LLM **sólo traduce** lenguaje natural → spec; el modelo lo construye el cliente de forma
determinista. Las API keys viven como **secretos del Worker** (nunca en el código ni en el
navegador):

```bash
npx wrangler secret put OPENROUTER_API_KEY   # u OPENAI_API_KEY
npx wrangler deploy
```

## Mantenimiento

- **Actualizar el motor:** `cd vendor/portico-core && git pull && cd ../.. && git add vendor/portico-core`.
  Tras subir el pin, **re-sincronizar** `assistant/` (datos) y `examples/` desde
  `vendor/portico-core/` y re-aplicar el preset Chile (ver [`CONTRIBUTING.md`](CONTRIBUTING.md)).
- **Normativa:** los datos NCh provienen de `presets/chile/` del core, copiados sobre
  `assistant/`. El motor es agnóstico; sólo estos datos definen la jurisdicción.

## Licencia

**GNU AGPL-3.0** (ver [`LICENSE`](LICENSE)), igual que portico-core. Usted puede usar,
estudiar, modificar y redistribuir el software, y tiene derecho a obtener el código fuente
correspondiente, incluso al usar la aplicación a través de la red (AGPL §13).

© 2026 Instituto de Obras Civiles, Universidad Austral de Chile · Dr. Juan Patricio Reyes C.

Bibliotecas de terceros: [Three.js](https://github.com/mrdoob/three.js) (MIT) ·
[numeric.js](https://github.com/sloisel/numeric) (MIT).
