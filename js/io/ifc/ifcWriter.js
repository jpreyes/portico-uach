// ──────────────────────────────────────────────────────────────────────────────
// io/ifc/ifcWriter.js — exportador IFC (STEP / ISO-10303-21) · #75-#77, G19
//
// Escribe el modelo NEUTRO de `io/` como un .ifc (IFC4) de TEXTO, sin dependencias —
// la contraparte del parser de `ifcLoader.js`.  Cada barra se emite como IfcBeam (o
// IfcColumn si es casi vertical) con:
//   • la jerarquía espacial mínima  IfcProject → IfcSite → IfcBuilding → IfcBuildingStorey;
//   • su EJE como representación 'Axis' (IfcPolyline de los dos extremos, en metros
//     globales) → lo que nuestro propio importador vuelve a leer (round-trip);
//   • material (IfcMaterial + IfcMechanicalMaterialProperties con E en Pa) y sección
//     (IfcRectangleProfileDef dimensionada para REPRODUCIR A e Iz exactos).
// Unidades de salida: metros (factor 1).  AUTÓNOMO (Node + navegador).
// ──────────────────────────────────────────────────────────────────────────────

// GUID IFC (22 chars del alfabeto base64 de IFC; el 1.º acotado a 2 bits). No se exige
// unicidad estricta a los visores; basta que sean válidos y distintos.
const G64 = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$';
function ifcGuid() {
  let s = G64[Math.floor(Math.random() * 4)];
  for (let i = 1; i < 22; i++) s += G64[Math.floor(Math.random() * 64)];
  return s;
}

// álgebra vectorial mínima (para superficies de área)
const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const vmul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const vcross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const vlen = a => Math.hypot(a[0], a[1], a[2]);
const vunit = a => { const l = vlen(a); return l > 1e-12 ? vmul(a, 1 / l) : [0, 0, 0]; };

// real IFC: siempre con punto decimal (1 → "1.")
function num(v) { let s = (+v || 0).toString(); if (!/[.eE]/.test(s)) s += '.'; return s; }
// string IFC: comilla escapada '' y no-ASCII como \X2\HHHH…\X0\ (compatibilidad amplia)
function str(s) {
  s = String(s ?? '');
  let out = '', i = 0;
  while (i < s.length) {
    const c = s.charCodeAt(i);
    if (c > 0x7e || c < 0x20) {
      let hex = '';
      while (i < s.length && (s.charCodeAt(i) > 0x7e || s.charCodeAt(i) < 0x20)) { hex += s.charCodeAt(i).toString(16).toUpperCase().padStart(4, '0'); i++; }
      out += `\\X2\\${hex}\\X0\\`;
    } else { out += s[i] === "'" ? "''" : s[i]; i++; }
  }
  return `'${out}'`;
}

/**
 * Modelo neutro de `io/` → texto IFC4 (.ifc).
 * @param {object} neutral  salida de `modelToNeutral`
 * @param {object} [opts]   { name }
 * @returns {string}
 */
export function neutralToIFC(neutral, opts = {}) {
  const name = opts.name || neutral.meta?.name || 'PORTICO';
  const lines = [];
  let id = 0;
  const e = (body) => { const ref = '#' + (++id); lines.push(`${ref}=${body};`); return ref; };

  // ── contexto, unidades, proyecto ──────────────────────────────────────────────
  const origin = e('IFCCARTESIANPOINT((0.,0.,0.))');
  const cs = e(`IFCAXIS2PLACEMENT3D(${origin},$,$)`);
  const lenUnit = e('IFCSIUNIT(*,.LENGTHUNIT.,$,.METRE.)');
  const angUnit = e('IFCSIUNIT(*,.PLANEANGLEUNIT.,$,.RADIAN.)');
  const unitAssign = e(`IFCUNITASSIGNMENT((${lenUnit},${angUnit}))`);
  const ctx = e(`IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,${cs},$)`);
  const project = e(`IFCPROJECT(${str(ifcGuid())},$,${str(name)},$,$,$,$,(${ctx}),${unitAssign})`);

  // ── jerarquía espacial: Site → Building → Storey (placements identidad) ──
  const sitePl = e(`IFCLOCALPLACEMENT($,${cs})`);
  const bldgPl = e(`IFCLOCALPLACEMENT(${sitePl},${cs})`);
  const storeyPl = e(`IFCLOCALPLACEMENT(${bldgPl},${cs})`);
  const site = e(`IFCSITE(${str(ifcGuid())},$,'Site',$,$,${sitePl},$,$,.ELEMENT.,$,$,$,$,$)`);
  const building = e(`IFCBUILDING(${str(ifcGuid())},$,'Building',$,$,${bldgPl},$,$,.ELEMENT.,$,$,$)`);
  const storey = e(`IFCBUILDINGSTOREY(${str(ifcGuid())},$,'Planta 1',$,$,${storeyPl},$,$,.ELEMENT.,0.)`);
  e(`IFCRELAGGREGATES(${str(ifcGuid())},$,$,$,${project},(${site}))`);
  e(`IFCRELAGGREGATES(${str(ifcGuid())},$,$,$,${site},(${building}))`);
  e(`IFCRELAGGREGATES(${str(ifcGuid())},$,$,$,${building},(${storey}))`);

  // ── materiales (IfcMaterial + propiedades mecánicas: E en Pa) ──
  const matRef = new Map();
  for (const m of (neutral.materials || [])) {
    const mr = e(`IFCMATERIAL(${str(m.name || 'Material')})`);
    matRef.set(m.id, mr);
    const E = (m.E || 0) * 1000, G = (m.G || 0) * 1000;   // kN/m² → Pa
    if (E > 0) e(`IFCMECHANICALMATERIALPROPERTIES(${mr},$,${num(E)},${num(G)},${num(m.nu ?? 0.3)},${num(m.alpha ?? 1.2e-5)})`);
  }
  const secById = new Map((neutral.sections || []).map(s => [s.id, s]));

  // sección → rectángulo (b,h) que reproduce A e Iz exactos:  h=√(12·Iz/A), b=A/h
  const rectOf = (s) => {
    const A = s.A > 0 ? s.A : 1e-3;
    let h = (s.Iz > 0 && A > 0) ? Math.sqrt(12 * s.Iz / A) : Math.sqrt(A);
    if (!isFinite(h) || h <= 0) h = Math.sqrt(A);
    let b = A / h; if (!isFinite(b) || b <= 0) b = Math.sqrt(A);
    return { b, h };
  };

  // ── barras: punto×2 → polilínea 'Axis' → IfcBeam/IfcColumn ──
  const memPl = e(`IFCLOCALPLACEMENT(${storeyPl},${cs})`);
  const nodeById = new Map((neutral.nodes || []).map(n => [n.id, n]));
  const elemsByMatSec = new Map();   // "mat|sec" → [elemRef, …]
  const allElems = [];
  for (const mb of (neutral.members || [])) {
    const a = nodeById.get(mb.ni), b = nodeById.get(mb.nj);
    if (!a || !b) continue;
    const p1 = e(`IFCCARTESIANPOINT((${num(a.x)},${num(a.y)},${num(a.z)}))`);
    const p2 = e(`IFCCARTESIANPOINT((${num(b.x)},${num(b.y)},${num(b.z)}))`);
    const pl = e(`IFCPOLYLINE((${p1},${p2}))`);
    const sr = e(`IFCSHAPEREPRESENTATION(${ctx},'Axis','Curve3D',(${pl}))`);
    const ps = e(`IFCPRODUCTDEFINITIONSHAPE($,$,(${sr}))`);
    const L = Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z) || 1;
    const vertical = Math.abs(b.z - a.z) / L > 0.85;
    const tipo = vertical ? 'IFCCOLUMN' : 'IFCBEAM';
    const pdt = vertical ? '.COLUMN.' : '.BEAM.';
    const elName = `${vertical ? 'Pilar' : 'Viga'} ${mb.id}`;
    const el = e(`${tipo}(${str(ifcGuid())},$,${str(elName)},$,$,${memPl},${ps},$,${pdt})`);
    allElems.push(el);
    const key = `${mb.mat || 1}|${mb.sec || 1}`;
    if (!elemsByMatSec.has(key)) elemsByMatSec.set(key, []);
    elemsByMatSec.get(key).push(el);
  }

  // ── asociación material+perfil por combinación (mat,sec) usada ──
  for (const [key, els] of elemsByMatSec) {
    const [mi, si] = key.split('|').map(Number);
    const sec = secById.get(si);
    const mr = matRef.get(mi) || matRef.values().next().value;
    const { b, h } = rectOf(sec || { A: 1e-3, Iz: 0 });
    const prof = e(`IFCRECTANGLEPROFILEDEF(.AREA.,${str(sec?.name || 'Sección')},$,${num(b)},${num(h)})`);
    const mProf = e(`IFCMATERIALPROFILE(${str(sec?.name || 'Sección')},$,${mr || '$'},${prof},$,$)`);
    const mSet = e(`IFCMATERIALPROFILESET(${str(sec?.name || 'Sección')},$,(${mProf}),$)`);
    const mUse = e(`IFCMATERIALPROFILESETUSAGE(${mSet},$,$)`);
    e(`IFCRELASSOCIATESMATERIAL(${str(ifcGuid())},$,$,$,(${els.join(',')}),${mUse})`);
  }

  // ── áreas: muro/losa como sólido extruido (contorno medio × espesor) ──
  const nById = new Map((neutral.nodes || []).map(n => [n.id, [n.x, n.y, n.z]]));
  for (const ar of (neutral.areas || [])) {
    const cs = (ar.nodes || []).map(i => nById.get(i)).filter(Boolean);
    if (cs.length < 3) continue;
    const t = ar.thickness || 0.2;
    let nrm = vunit(vcross(vsub(cs[1], cs[0]), vsub(cs[2], cs[0])));
    if (vlen(nrm) < 1e-9) continue;
    const vertical = Math.abs(nrm[2]) <= 0.7;
    const lx = vunit(vsub(cs[1], cs[0])), lz = nrm, ly = vcross(lz, lx);
    const base0 = vsub(cs[0], vmul(lz, t / 2));                 // plano base del sólido (centra en la superficie media)
    const poly2d = cs.map(c => [vdot(vsub(c, cs[0]), lx), vdot(vsub(c, cs[0]), ly)]);
    const ptRefs = poly2d.map(q => e(`IFCCARTESIANPOINT((${num(q[0])},${num(q[1])}))`));
    const pl = e(`IFCPOLYLINE((${ptRefs.join(',')},${ptRefs[0]}))`);   // contorno cerrado
    const prof = e(`IFCARBITRARYCLOSEDPROFILEDEF(.AREA.,${str(ar.name || 'Área')},${pl})`);
    const oPt = e(`IFCCARTESIANPOINT((${num(base0[0])},${num(base0[1])},${num(base0[2])}))`);
    const zD = e(`IFCDIRECTION((${num(lz[0])},${num(lz[1])},${num(lz[2])}))`);
    const xD = e(`IFCDIRECTION((${num(lx[0])},${num(lx[1])},${num(lx[2])}))`);
    const ap = e(`IFCAXIS2PLACEMENT3D(${oPt},${zD},${xD})`);
    const solid = e(`IFCEXTRUDEDAREASOLID(${prof},${ap},${e('IFCDIRECTION((0.,0.,1.))')},${num(t)})`);
    const sr = e(`IFCSHAPEREPRESENTATION(${ctx},'Body','SweptSolid',(${solid}))`);
    const ps = e(`IFCPRODUCTDEFINITIONSHAPE($,$,(${sr}))`);
    const tipo = vertical ? 'IFCWALL' : 'IFCSLAB';
    const pdt = vertical ? '.STANDARD.' : '.FLOOR.';
    const el = e(`${tipo}(${str(ifcGuid())},$,${str((vertical ? 'Muro ' : 'Losa ') + ar.id)},$,$,${memPl},${ps},$,${pdt})`);
    allElems.push(el);
    const mr = matRef.get(ar.mat) || matRef.values().next().value;
    if (mr) {
      const ml = e(`IFCMATERIALLAYER(${mr},${num(t)},$)`);
      const mls = e(`IFCMATERIALLAYERSET((${ml}),${str((vertical ? 'Muro ' : 'Losa ') + ar.id)})`);
      const mlu = e(`IFCMATERIALLAYERSETUSAGE(${mls},.AXIS2.,.POSITIVE.,0.)`);
      e(`IFCRELASSOCIATESMATERIAL(${str(ifcGuid())},$,$,$,(${el}),${mlu})`);
    }
  }

  // ── contención de TODOS los elementos (barras + áreas) en el piso ──
  if (allElems.length)
    e(`IFCRELCONTAINEDINSPATIALSTRUCTURE(${str(ifcGuid())},$,'Elementos',$,(${allElems.join(',')}),${storey})`);

  // ── ensamblar el archivo ──
  const ts = new Date().toISOString().slice(0, 19);
  const header =
    `ISO-10303-21;\nHEADER;\n` +
    `FILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1');\n` +
    `FILE_NAME(${str(name + '.ifc')},'${ts}',(''),(''),'PORTICO','PORTICO','');\n` +
    `FILE_SCHEMA(('IFC4'));\nENDSEC;\nDATA;\n`;
  return header + lines.join('\n') + '\nENDSEC;\nEND-ISO-10303-21;\n';
}
