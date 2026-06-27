# Cambios — hormigón (NCh170), textos de UI y copy de capacidades

Registro de los cambios aplicados en esta tanda. Versión de la app: **v212**.

---

## 1. Hormigón: nomenclatura H → G (NCh170:2016)

La nomenclatura chilena del hormigón cambió de **«H»** (designación por resistencia
**cúbica**) a **«G»** (grado por resistencia **cilíndrica** característica `f'c`, NCh170:2016).
En PÓRTICO el **número del grado ya correspondía a `f'c`** (p.ej. el material «H30» usaba
`f'c = 30 MPa`, `E = 4700·√f'c`), por lo que el cambio es un **renombre directo** sin alterar
propiedades: `G30 = f'c 30 MPa`, mismos `E`, `G`, `ν`, `ρ`.

### Qué se cambió
| Lugar | Antes | Ahora |
|---|---|---|
| Material por defecto del modelo (`js/model/model.js`) | `Concreto H30` | `Hormigón G30` |
| Catálogo de materiales (`js/design/materials_catalog.js`) | Hormigón H25/H30/H40 | **Hormigón G20/G25/G30/G40** (+ familia) |
| Catálogo del asistente (`asistente/materiales.csv`) | H20…H50 | **G20…G50** |
| Plantilla CSV (`js/model/serializer.js`) | `Concreto H30` | `Hormigón G30` |
| Macromodelo muro de relleno (comentario) | col 30×30 H25 | col 30×30 G25 |
| Esquema de la ficha (`asistente/ficha.schema.json`) | H20/H25/H30/H40, `H{fc}` | **G20/G25/G30/G40, `G{fc}`** |
| Prompt del asistente (`worker/asistente.js`, `probar_pipeline.mjs`, `n8n_flujo.json`) | «hormigon Hxx» | «hormigón grado **G** (NCh170:2016)» |
| Ejemplos y tests del asistente | H30/H40/H50 | G30/G40/G50 |
| Placeholders del asistente (index.html / app.js) | «hormigón H30» | «hormigón **G30**» |

### Compatibilidad hacia atrás
El generador del asistente (`asistente/generador.js`) **acepta tanto «Gxx» como el legado
«Hxx»**: el regex de contexto reconoce `[hg]\d`, y la búsqueda hace `byName.get('g'+fc) ||
byName.get('h'+fc)`. Así, si alguien escribe «H30» o «G30», ambos resuelven al material **G30**.
Los archivos `.s3d` antiguos **no se modifican** (el nombre del material es sólo un string que
se conserva tal cual al abrir).

Verificado: `G30`, `hormigón G30`, `H30` (legado), `fc=30` → todos resuelven a **G30** con
`E = 2.87·10⁷ kN/m²`.

---

## 2. Textos de UI (propuestas del docente en `docs/textos_ui.csv`)

Se aplicaron las 9 propuestas escritas en la columna `texto_propuesto`:

| Ubicación | Antes | Ahora |
|---|---|---|
| Logo del header | `PÓRTICO  IOC · UACh` | `PÓRTICO` (se quitó «IOC · UACh») |
| Menú Archivo | `📝 Memoria de Cálculo (Word .docx)…` | `📄 Memoria / Bases de Cálculo (Word)…` |
| Menú Editar | `Suavizar malla (calidad)…` | `Suavizar malla…` |
| Menú Análisis | `▶ Time-history modal — acelerograma…` | `▶ Time-history modal…` |
| Menú Análisis | `▶ Rótulas Plásticas — colapso…` | `▶ Pushover — colapso…` |
| Menú Análisis | `▶ Time-history NO LINEAL — rótulas (edificio de corte)…` | `▶ Time-history NO LINEAL — Modal…` |
| Placeholders del asistente | `…hormigón H30…` | `…hormigón G30…` |

> Nota: `docs/textos_ui.csv` queda como tu documento de trabajo. Tras aplicar estas
> propuestas puedes limpiar esas filas; si re-guardas y corres `node docs/refresh_ui_texts.mjs`
> se refrescan los números de línea (cambió `index.html`).

---

## 3. Pestañas de resultados «Pushover» diferenciadas

Al renombrar la pestaña **«Rótulas» → «Pushover»** quedaba duplicada con la pestaña de
pushover por control de desplazamiento que ya existía. Se diferenciaron:

- `data-rtab="plastico"` (colapso incremental por formación de rótulas, control de carga) →
  **«Pushover (colapso)»**.
- `data-rtab="dc"` (control de desplazamiento) → **«Pushover (DC)»**.

Los ítems de menú ya estaban diferenciados: **«Pushover — colapso…»** vs
**«Pushover — control de desplazamiento…»**.

---

## 4. Copy de capacidades (la app dejaba de parecer «solo barras»)

El texto de presentación sub-vendía la app y las **«Limitaciones conocidas»** del *Acerca de*
estaban **obsoletas y eran falsas** (decían «sin verificación de diseño», «solo elementos
barra», «solo cargas uniformes», «análisis lineal»). Se reescribió todo el copy:

- **Landing (splash):** tagline → «análisis **y diseño** estructural 3D»; descripción y chips
  con la gama real (barras + áreas, espectro NCh433, pandeo·P-Δ, no lineal·pushover, diseño
  acero/hormigón/madera/aluminio, derivas, IFC/BIM, asistente IA).
- **Acerca de:** «Capacidades» reescritas y completas (modelado, apoyos/conexiones, cargas,
  análisis, **diseño multinorma**, interoperabilidad, productividad). «Limitaciones» →
  **«Alcances y limitaciones»** con las limitaciones REALES (herramienta docente; diseño
  orientativo; áreas lineal-elásticas; dinámica NL de alta fidelidad en integración).
- **`<head>`:** título → «Análisis y Diseño»; se añadió `meta description` + Open Graph
  (`og:title`/`og:description`/`og:type`) para compartir/buscar.
- **`README.md`:** nueva sección «Capacidades» + nota docente; ejemplo de material → Hormigón G30.

---

## Verificación

- Batería de tests verde: IFC (40), solver `#2…#6` (masa rotacional, resortes acoplados,
  Guyan, uplift, suelo p-y), y tests del asistente (`test_generador`, `test_primitivas`,
  `test_puente_galpon`).
- `node --input-type=module --check` OK en los módulos tocados; HTML del *Acerca de* balanceado.
- Bump de versión global a **v212** (cache-busting).
