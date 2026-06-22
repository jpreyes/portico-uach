// build_barqueta.mjs — Puente de la BARQUETA (Sevilla, Arenas & Pantaleón, 1992).
// Modelo 3D: arco central único, RED de péndolas inclinadas centradas (network),
// pórticos triangulares de extremo, y TABLERO modelado con elementos de ÁREA
// (shell = membrana + placa) que además actúa de tirante.
//   node tools/build_barqueta.mjs
import fs from 'fs';
import path from 'path';
import { execFileSync } from 'child_process';
import { Model } from '../js/model/model.js';
import { Serializer } from '../js/model/serializer.js';
import { runStatic } from './verif/runners.mjs';
import { renderModelSVG } from './verif/figure.mjs';

const ROOT = process.cwd();
const LOGOS = 'icons/UACh-color-negro.svg,icons/Facultad-color-negro.svg,icons/IOC-color.svg';

const L = 168, W = 18, NX = 21, NY = 2, hp = 5, rise = 24;   // luz, ancho, malla, apoyo del arco, flecha
const dx = L / NX, dyc = W / NY;
const zArch = (x) => hp + rise * (1 - ((x - L / 2) / (L / 2)) ** 2);

const m = new Model(); m.mode = '3D'; m.materials.clear(); m.sections.clear();
const steel = m.addMaterial({ name: 'Acero', E: 2.0e8, G: 7.7e7, nu: 0.3, rho: 7.85 }).id;
const conc  = m.addMaterial({ name: 'Hormigón losa', E: 3.0e7, G: 1.25e7, nu: 0.2, rho: 2.5 }).id;
const sArch = m.addSection({ name: 'Arco (cajón acero)', A: 0.6, Iy: 0.30, Iz: 0.30, J: 0.4, Avy: 0.3, Avz: 0.3, kappay: 0.6, kappaz: 0.6 }).id;
const sLeg  = m.addSection({ name: 'Pata de pórtico', A: 0.4, Iy: 0.15, Iz: 0.15, J: 0.2, Avy: 0.2, Avz: 0.2, kappay: 0.6, kappaz: 0.6 }).id;
const sHang = m.addSection({ name: 'Péndola', A: 0.012, Iy: 1e-4, Iz: 1e-4, J: 1e-5, Avy: 0.006, Avz: 0.006, kappay: 0.9, kappaz: 0.9 }).id;
const tDeck = 0.4;   // espesor equivalente del tablero (losa ortótropa)

// ── Tablero: malla de nodos (NX+1)×(NY+1) en z=0 → áreas shell ───────────────────
const dn = [];   // dn[i][j] = nodeId
for (let i = 0; i <= NX; i++) {
  dn.push([]);
  for (let j = 0; j <= NY; j++) {
    const x = i * dx, y = j * dyc;
    // apoyos en las 4 esquinas (uno fijo en planta, los demás liberan según dirección)
    let r = {};
    if ((i === 0 || i === NX) && (j === 0 || j === NY)) {
      if (i === 0 && j === 0) r = { ux: 1, uy: 1, uz: 1 };         // articulado fijo
      else if (i === 0) r = { uy: 1, uz: 1 };                       // libera X (dilatación long.)
      else if (j === 0) r = { ux: 1, uz: 1 };                       // libera Y
      else r = { uz: 1 };                                           // rodillo
    }
    dn[i].push(m.addNode(x, y, 0, r).id);
  }
}
// áreas shell (QUAD) del tablero
for (let i = 0; i < NX; i++) for (let j = 0; j < NY; j++)
  m.addArea([dn[i][j], dn[i + 1][j], dn[i + 1][j + 1], dn[i][j + 1]], conc, { thickness: tDeck, behavior: 'shell' });

// ── Arco central (plano y=W/2) ───────────────────────────────────────────────────
const jc = NY / 2;   // fila central (centerline)
const arch = [];
for (let i = 0; i <= NX; i++) {
  const x = i * dx;
  // arranques (i=0, NX) en el apoyo del pórtico (z=hp); el resto sobre la parábola
  arch.push(m.addNode(x, W / 2, zArch(x)).id);
}
for (let i = 0; i < NX; i++) m.addElement(arch[i], arch[i + 1], steel, sArch);

// ── Pórticos triangulares de extremo: del arranque del arco a las 2 esquinas ────
for (const i of [0, NX]) {
  m.addElement(arch[i], dn[i][0], steel, sLeg);
  m.addElement(arch[i], dn[i][NY], steel, sLeg);
  m.addElement(arch[i], dn[i][jc], steel, sLeg);   // y al eje del tablero
}

// ── RED de péndolas centradas (network): inclinadas y cruzadas ───────────────────
const off = 2;   // desfase → péndolas inclinadas que se cruzan
for (let i = 1; i < NX; i++) {
  for (const t of [i - off, i + off]) if (t >= 1 && t <= NX - 1) m.addElement(arch[i], dn[t][jc], steel, sHang);
}

// ── Cargas: peso propio (barras) + carga de tablero como fuerzas nodales ─────────
const lc = m.addLoadCase('PP + tablero + tránsito', true);
const qDeck = 12;   // kN/m² (losa + pavimento + tránsito)
for (let i = 0; i <= NX; i++) for (let j = 0; j <= NY; j++) {
  const n = m.nodes.get(dn[i][j]); if (Object.values(n.restraints).some(v => v)) continue;
  const trib = (i === 0 || i === NX ? dx / 2 : dx) * (j === 0 || j === NY ? dyc / 2 : dyc);
  m.addLoad(lc.id, { type: 'nodal', nodeId: dn[i][j], F: [0, 0, -qDeck * trib, 0, 0, 0] });
}

const res = await runStatic(m, lc.id, true);
let Rz = 0; for (const nd of m.nodes.values()) { const r = nd.restraints; if (r && (r.ux || r.uy || r.uz)) Rz += res.getReaction(nd.id)[2]; }
const sum = res.getSummary();

// figura 3D (áreas + barras + deformada)
const nodes = new Map(), elements = [], supports = new Set();
for (const nd of m.nodes.values()) { nodes.set(nd.id, [nd.x, nd.y, nd.z]); const r = nd.restraints; if (r && (r.ux || r.uy || r.uz)) supports.add(nd.id); }
for (const e of m.elements.values()) elements.push({ n1: e.n1, n2: e.n2 });
const areas = [...m.areas.values()].map(a => a.nodes);
let mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
for (const c of nodes.values()) for (let k = 0; k < 3; k++) { mn[k] = Math.min(mn[k], c[k]); mx[k] = Math.max(mx[k], c[k]); }
const diag = Math.hypot(mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]) || 1;
let maxT = 0; const defo = new Map();
for (const id of nodes.keys()) { const d = res.getNodeDisp(id); defo.set(id, d); maxT = Math.max(maxT, Math.hypot(d[0], d[1], d[2])); }
const amp = maxT > 0 ? 0.10 * diag / maxT : 0; const deformed = new Map();
for (const [id, c] of nodes) { const d = defo.get(id); deformed.set(id, [c[0] + amp * d[0], c[1] + amp * d[1], c[2] + amp * d[2]]); }
fs.mkdirSync(path.join(ROOT, 'docs/ejemplos/img'), { recursive: true });
fs.writeFileSync(path.join(ROOT, 'docs/ejemplos/img/puente_barqueta.svg'), renderModelSVG({ nodes, elements, areas, deformed, supports, width: 980 }), 'utf8');
fs.writeFileSync(path.join(ROOT, 'examples', 'puente_barqueta.s3d'), new Serializer().toJSON(m), 'utf8');

const mdTable = (h, rows) => `| ${h.join(' | ')} |\n| ${h.map(() => '---').join(' | ')} |\n` + rows.map(r => `| ${r.join(' | ')} |`).join('\n') + '\n';
const md = `# Puente de la Barqueta (Sevilla, 1992) — arco con red de péndolas y tablero de áreas

**Tipo:** ejemplo 3D con **geometría real** y **tablero modelado con elementos de área (shell)** · **Modelo:** [\`examples/puente_barqueta.s3d\`](../../examples/puente_barqueta.s3d)

## Descripción

El **Puente de la Barqueta** (Sevilla, Juan José Arenas y Marcos J. Pantaleón, Expo'92) salva el Guadalquivir con **un solo vano de 168 m**. Es un **arco atirantado** con **un único arco central** del que cuelga el tablero mediante una **red de péndolas inclinadas centradas** (los característicos tirantes rojos cruzados), y con **pórticos triangulares** en cada extremo (las «puertas») que recogen el arranque del arco. El tablero (mixto acero-hormigón) actúa además de **tirante**, cerrando el empuje del arco.

${mdTable(['Propiedad', 'Valor'], [['Luz', '168 m (vano único)'], ['Ancho del tablero (modelo)', `${W} m`], ['Arco', 'único, central, flecha ~24 m'], ['Péndolas', 'red inclinada centrada (network)'], ['Extremos', 'pórticos triangulares («puertas»)'], ['Tablero', 'elementos de ÁREA (shell), tirante'], ['Autores / año', 'Arenas & Pantaleón / 1992']])}
## Modelo en Pórtico

- El **tablero** se modela con **${m.areas.size} elementos de área (QUAD shell)** — membrana (acción de tirante en su plano) + placa (flexión transversal). Es la novedad pedida: el tablero como áreas, no como viga.
- El **arco** central (plano longitudinal medio) es una cadena de elementos viga; los **pórticos triangulares** de extremo conectan el arranque del arco con las esquinas del tablero.
- La **red de péndolas** se arma con elementos inclinados **cruzados** (cada nodo del arco baja a nodos del eje del tablero desfasados ±${off} → patrón network), reproduciendo los tirantes inclinados característicos.
- Apoyos en las **4 esquinas** del tablero (uno fijo en planta, los demás liberan la dilatación); el tablero-tirante absorbe el **empuje** del arco.

![Puente de la Barqueta](img/puente_barqueta.svg)

*Figura. Vista 3D: tablero (áreas), arco central, pórticos triangulares y red de péndolas, con la deformada (×escala) bajo peso propio + carga de tablero.*

## Resultados (peso propio + carga de tablero ${qDeck} kN/m²)

${mdTable(['Magnitud', 'Valor'], [
  ['Nodos · elementos · áreas', `${m.nodes.size} · ${m.elements.size} · ${m.areas.size}`],
  ['ΣReacciones verticales', `${Rz.toFixed(0)} kN`],
  ['Desplazamiento máx. |u|', `${(sum.maxU * 1000).toFixed(1)} mm`],
  ['Axial máx. |N| (arco/tirante)', `${sum.maxN.toFixed(0)} kN`],
])}
## Conclusión

El modelo reproduce la **forma real** de la Barqueta —arco único central, **red de péndolas inclinadas centradas** y pórticos triangulares de extremo— con el **tablero modelado por elementos de área (shell)** que trabaja de tirante. Resuelve en equilibrio bajo peso propio + carga de tablero. Ejemplo avanzado que combina **barras + áreas** en un puente real. *(El cálculo riguroso de las péndolas usa el análisis geométrico/no lineal de Pórtico — Kg/NL-lite.)*
`;
const mdPath = path.join(ROOT, 'docs/ejemplos/puente_barqueta.md');
fs.writeFileSync(mdPath, md, 'utf8');
execFileSync('node', ['tools/md2pdf.mjs', path.relative(ROOT, mdPath), '--logos', LOGOS], { cwd: ROOT, stdio: 'ignore' });
console.log(`✓ puente_barqueta  ·  ${m.nodes.size} nodos, ${m.elements.size} elem, ${m.areas.size} áreas  ·  ΣRz=${Rz.toFixed(0)} kN  ·  umax=${(sum.maxU * 1000).toFixed(1)} mm  ·  Nmax=${sum.maxN.toFixed(0)} kN`);
