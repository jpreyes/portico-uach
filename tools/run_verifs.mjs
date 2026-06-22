// run_verifs.mjs — pipeline BATCHED de verificaciones.
//   node tools/run_verifs.mjs            # todos los casos de tools/verif/cases/
//   node tools/run_verifs.mjs 1-014      # sólo los que matcheen
//
// Por cada caso: carga el .s3d (construido a mano) → corre el solver de Pórtico
// HEADLESS → genera la figura 3D (SVG isométrico) → compara contra la referencia
// → arma el .md → genera el .pdf con md2pdf (membrete IOC). Lo único caso-a-caso
// es construir el .s3d y validar; todo lo demás es automático.
import fs from 'fs';
import path from 'path';
import { execFileSync } from 'child_process';
import { Serializer } from '../js/model/serializer.js';
import { runAnalysis } from './verif/runners.mjs';
import { renderModelSVG } from './verif/figure.mjs';

const ROOT = process.cwd();
const LOGOS = 'icons/UACh-color-negro.svg,icons/Facultad-color-negro.svg,icons/IOC-color.svg';
const pat = process.argv[2] || '';

const pct = (v, ref) => { const d = (v - ref) / ref * 100; return Math.abs(d) < 0.005 ? '0 %' : `${d >= 0 ? '+' : ''}${d.toFixed(2)} %`; };
const esc = s => String(s);

function mdTable(header, rows) {
  return `| ${header.join(' | ')} |\n| ${header.map(() => '---').join(' | ')} |\n`
    + rows.map(r => `| ${r.join(' | ')} |`).join('\n') + '\n';
}

// Figura: nodos/elementos + deformada de la forma modal (amp auto).
function buildFigure(model, out, caseDef) {
  const nodes = new Map(), elements = [], supports = new Set();
  for (const n of model.nodes.values()) {
    nodes.set(n.id, [n.x, n.y, n.z]);
    // Apoyo "real" = ≥2 restricciones de traslación (pin/empotramiento). Evita
    // marcar restricciones artificiales de un solo GDL (p.ej. Ux para excluir axial).
    const r = n.restraints; if (r && ((r.ux ? 1 : 0) + (r.uy ? 1 : 0) + (r.uz ? 1 : 0)) >= 2) supports.add(n.id);
  }
  for (const e of model.elements.values()) elements.push({ n1: e.n1, n2: e.n2 });
  // diagonal del bbox para escalar la amplitud
  let mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
  for (const c of nodes.values()) for (let k = 0; k < 3; k++) { mn[k] = Math.min(mn[k], c[k]); mx[k] = Math.max(mx[k], c[k]); }
  const diag = Math.hypot(mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]) || 1;
  let deformed = null;
  if (caseDef.figure && out.type === 'modal') {
    const amp = 0.18 * diag;   // la forma modal viene normalizada a máx. traslación = 1
    const shape = out.res.getModeShape((caseDef.figure.mode || 1) - 1);
    deformed = new Map();
    for (const [id, c] of nodes) { const s = shape.get(id) || [0, 0, 0]; deformed.set(id, [c[0] + amp * s[0], c[1] + amp * s[1], c[2] + amp * s[2]]); }
  }
  return renderModelSVG({ nodes, elements, supports, deformed, width: 900 });
}

function buildComparison(cmp, out) {
  const pv = cmp.portico(out.res);
  const header = ['Modo', 'Descripción', `Independiente (${cmp.unit})`, `SAP2000 (${cmp.unit})`, 'dif. SAP', `**Pórtico (${cmp.unit})**`, '**dif. Pórtico**'];
  const rows = cmp.rows.map((r, i) => {
    const p = pv[i];
    return [String(i + 1), r.desc, r.indep.toFixed(cmp.decimals), r.sap.toFixed(cmp.decimals), pct(r.sap, r.indep),
      `**${p.toFixed(cmp.decimals)}**`, `**${pct(p, r.indep)}**`];
  });
  return { table: mdTable(header, rows), pv };
}

async function runCase(file) {
  const mod = (await import('./verif/cases/' + file)).default;
  const model = new Serializer().fromJSON(fs.readFileSync(path.join(ROOT, mod.s3d), 'utf8'));
  const out = await runAnalysis(model, mod);

  // figura
  const svg = buildFigure(model, out, mod);
  const imgRel = `img/${mod.slug}.svg`;
  fs.mkdirSync(path.join(ROOT, 'docs/verificaciones/img'), { recursive: true });
  fs.writeFileSync(path.join(ROOT, 'docs/verificaciones', imgRel), svg, 'utf8');

  // comparación + sustitución de placeholders {{Pi}}/{{Di}} en `extra`
  const { table, pv } = buildComparison(mod.compare, out);
  let extra = mod.extra || '';
  extra = extra.replace(/\{\{P(\d+)\}\}/g, (_, i) => pv[+i].toFixed(mod.compare.decimals))
               .replace(/\{\{D(\d+)\}\}/g, (_, i) => pct(pv[+i], mod.compare.rows[+i].indep));
  const caption = typeof mod.figure?.caption === 'function' ? mod.figure.caption(out.res) : (mod.figure?.caption || '');

  const md = `# Verificación ${mod.id} — ${mod.title}

**Capacidad verificada:** ${mod.capability}.
**Referencia:** ${mod.referenceText}
**Modelo Pórtico:** [\`${mod.s3d}\`](../../${mod.s3d})

## Descripción del problema

${mod.intro}

${mdTable(['Propiedad', 'Valor'], mod.props)}
## Modelo en Pórtico

${mod.modelNotes.map(n => `- ${n}`).join('\n')}

![${esc(caption)}](${imgRel})

*Figura 1. ${caption}*

## Resultados — comparación

${mod.compare.intro}

${table}
${extra ? extra + '\n\n' : ''}## Conclusión

${mod.conclusion}
`;
  const mdPath = path.join(ROOT, 'docs/verificaciones', mod.slug + '.md');
  fs.writeFileSync(mdPath, md, 'utf8');

  // PDF
  execFileSync('node', ['tools/md2pdf.mjs', path.relative(ROOT, mdPath), '--logos', LOGOS], { cwd: ROOT, stdio: 'ignore' });

  // resumen a consola
  const maxDiff = Math.max(...mod.compare.rows.map((r, i) => Math.abs((pv[i] - r.indep) / r.indep * 100)));
  console.log(`✓ ${mod.id}  ${mod.slug}  ·  máx |dif| = ${maxDiff.toFixed(3)} %  ·  ${mod.slug}.pdf`);
  return { id: mod.id, maxDiff };
}

const files = fs.readdirSync(path.join(ROOT, 'tools/verif/cases')).filter(f => f.endsWith('.mjs') && (!pat || f.includes(pat))).sort();
if (!files.length) { console.error('Sin casos en tools/verif/cases/ que matcheen', pat); process.exit(1); }
console.log(`Corriendo ${files.length} caso(s)…`);
for (const f of files) { try { await runCase(f); } catch (e) { console.error(`✗ ${f}: ${e.message}`); } }
console.log('Listo.');
