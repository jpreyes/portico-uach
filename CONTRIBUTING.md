# Contribuir a PÓRTICO-UACh

PÓRTICO-UACh es la **capa académica (overlay)** de la Universidad Austral de Chile sobre el
motor open source **[portico-core](https://github.com/jpreyes/portico-core)**. Antes de
contribuir, decide **dónde** va tu cambio.

## ¿Dónde va el cambio?

| Tipo de cambio | Repositorio |
|---|---|
| Motor de cálculo, solver, UI genérica, formato `.s3d`, diseño, IO, asistente (generador determinista) | **portico-core** (PR upstream) |
| Marca/branding UACh, logos, textos institucionales, bloque legal, traducciones de esos textos | este repo (`branding.default.json`, `overlay/`, `index.html`) |
| Cloudflare Worker del asistente IA, secretos, despliegue | este repo (`worker/`, `wrangler.jsonc`) |
| Datos de normativa (NCh) | preset `presets/chile/` en **portico-core**; aquí sólo se copian sobre `assistant/` |

> **Regla de oro (dependencia unidireccional):** el overlay depende del core; el core
> **nunca** depende del overlay. **No edites `vendor/portico-core/`** — es un submódulo de
> sólo lectura. Si algo genérico necesita cambiar, hazlo como PR a portico-core y luego sube
> el pin del submódulo aquí.

## Cómo aporta el overlay (sin forkear)

El core es white-label y expone *seams* (ver
[`docs/EXTENDING.md`](https://github.com/jpreyes/portico-core/blob/main/docs/EXTENDING.md)).
Las contribuciones UACh se registran en `overlay/bootstrap.js`, que corre **antes** de
`app.js`:

- **Branding:** edita `branding.default.json` (nombre, tagline, logo, enlaces). El core lo
  aplica a los elementos `data-brand` de `index.html`.
- **Badges / config / análisis:** `extensions.registerBadge(…)`, `registerConfigSection(…)`,
  `setFlag(…)`, `registerAnalysis(…)` desde `overlay/bootstrap.js`.
- **i18n de textos UACh:** el motor traduce por *source string* (español = fuente). Para que
  un texto UACh tenga versión en inglés, escribe cada frase como **un único nodo de texto**
  (sin `<b>` a mitad de frase) e inyecta la traducción con `Object.assign(EN, { … })` en
  `overlay/bootstrap.js`.

## Convenciones

- **Idioma:** UI, comentarios y mensajes de commit en **español**. Cierra los commits con el
  trailer `Co-Authored-By: Claude`.
- **Sin build:** no hay bundler ni `package.json`; todo carga vía importmap y se sirve
  estático. Verifica con `python serve.py 8765`.
- **Cache-busting `?v=NNN`:** cada import/URL de los archivos **del overlay** lleva una
  versión global para invalidar cachés del navegador/SW; al publicar un cambio, súbela en
  todos los archivos del overlay a la vez (incl. `sw.js` → `CACHE_VERSION`). Los `?v=` del
  **core** vienen pinneados por el submódulo; no los tocas tú.
- **Nunca `git add -A`:** hay rutas intencionalmente sin versionar (`vendor/portico-core/` es
  un submódulo; `excel/`, `referencias/`, `node_modules/`, `CLAUDE.md` y demás archivos de
  agentes están en `.gitignore`). Stagea rutas explícitas.
- **Secretos:** las API keys (`OPENROUTER_API_KEY`/`OPENAI_API_KEY`, `ASIS_LOG_TOKEN`) viven
  en el dashboard de Cloudflare, **nunca** en el código.
- **Sintaxis ESM:** `node --input-type=module --check < archivo.js` (no uses `node --check`,
  trata el `.js` como CommonJS).

## Actualizar el motor (submódulo)

```bash
cd vendor/portico-core && git pull origin main && cd ../..
git add vendor/portico-core            # sube el pin del submódulo
```

Como el core hace `fetch` **relativo a la página** de sus datos, el overlay sirve copias en
la raíz. Tras subir el pin, **re-sincroniza**:

```bash
V=vendor/portico-core
# Datos del asistente (sin .js: el código se importa del submódulo):
cp $V/assistant/{rules.json,profiles.csv,materials.csv,live_loads.csv,design_params.json,spec.schema.json,examples.json,example_spec.json} assistant/
# Re-aplicar el preset Chile (NCh):
cp $V/presets/chile/{rules.json,design_params.json,live_loads.csv} assistant/
# Ejemplos guiados del core:
cp $V/examples/* examples/
```

Si cambió el contrato del spec (`spec.schema.json`), revisa el SYSTEM prompt de
`worker/asistente.js`. Corre los tests de verificación afectados del core
(`cd vendor/portico-core && node test_*.mjs`) antes de subir el pin.

## Despliegue

Producción auto-despliega desde `main` (Cloudflare). Para probar:

```bash
npx wrangler deploy
```

## Licencia

Al contribuir aceptas que tu aporte se distribuya bajo **GNU AGPL-3.0** (ver
[`LICENSE`](LICENSE)).
