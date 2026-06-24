// ──────────────────────────────────────────────────────────────────────────────
// io/ifc/ifcToPortico.js — IFC → modelo de PÓRTICO (vía modelo neutro) · #76, G19
//
// Orquesta loader + clasificador + geometría para producir, en DOS pasos separados:
//   1) analyzeIFC(text)  → `items`: una ficha por elemento (eje en global, sección y
//      material aproximados, nivel, estado y advertencias) SIN tocar el modelo.  Esto es
//      lo que alimenta la tabla y el preview side-by-side del diálogo (#77).
//   2) itemsToNeutral(items, seleccionados) → modelo NEUTRO de `io/` con los elementos
//      elegidos: crea nodos con SNAP por tolerancia, deduplica materiales/secciones y
//      arma las barras.  El modelo de PÓRTICO se construye con `neutralToModel`.
//
// Convención de ejes: IFC y PÓRTICO comparten Z VERTICAL, así que las coordenadas pasan
// directas (sólo conversión de unidades a metros).  AUTÓNOMO (Node + navegador).
// ──────────────────────────────────────────────────────────────────────────────
import { parseIFC, lengthUnit } from './ifcLoader.js?v=204';
import { classify, KIND_LABEL } from './ifcClassifier.js?v=204';
import { memberAxis, bodyProfile, profileProps, areaSurface } from './ifcGeometrySimplifier.js?v=204';
import { Warnings } from './ifcWarnings.js?v=204';
import { neutralToModel } from '../neutral.js?v=204';

const DEFAULT_TOL = 0.01;       // m — tolerancia de snap de nodos coincidentes
// material genérico (acero) cuando el IFC no trae propiedades mecánicas — kN/m²
const GENERIC = { E: 2.0e8, nu: 0.3, rho: 7.85, alpha: 1.2e-5 };

// ── resolución de material + perfil desde el RelatingMaterial asociado ────────────
function resolveMaterial(model, materialRef, mechE) {
  const out = { name: '', E: null, profile: null };
  let m = model.get(materialRef), guard = 0;
  while (m && guard++ < 8) {
    switch (m.type) {
      case 'IFCMATERIAL':
        out.name = out.name || (m.args[0] || '').toString();
        if (mechE.has(m.id)) out.E = mechE.get(m.id);
        return out;
      case 'IFCMATERIALLIST':
        m = model.get((m.args[0] || [])[0]); continue;
      case 'IFCMATERIALPROFILESETUSAGE':
        m = model.get(m.args[0]); continue;                 // ForProfileSet
      case 'IFCMATERIALPROFILESET':
        m = model.get((m.args[2] || [])[0]); continue;       // MaterialProfiles[0] (args: Name,Desc,Profiles,…)
      case 'IFCMATERIALPROFILE':                             // (Name, Desc, Material, Profile, …)
        if (m.args[3]) out.profile = m.args[3];              // el nombre de material lo da el IfcMaterial, no el perfil
        m = model.get(m.args[2]); continue;                  // → IfcMaterial (nombre/E)
      case 'IFCMATERIALLAYERSETUSAGE':                        // áreas: muros/losas por capas
        m = model.get(m.args[0]); continue;                  // ForLayerSet
      case 'IFCMATERIALLAYERSET':
        m = model.get((m.args[0] || [])[0]); continue;       // MaterialLayers[0]
      case 'IFCMATERIALLAYER':
        m = model.get(m.args[0]); continue;                  // → IfcMaterial
      default:
        return out;                                          // capas (muros) u otros → genérico
    }
  }
  return out;
}

// mapa  materialId → E (kN/m²)  desde IfcMechanicalMaterialProperties (IFC2x3)
function mechModulus(model) {
  const map = new Map();
  for (const p of model.ofType('IFCMECHANICALMATERIALPROPERTIES')) {
    // (Material, DynamicViscosity, YoungModulus[Pa], ShearModulus, PoissonRatio, ThermalExpansion)
    const mid = model.isRef(p.args[0]) ? p.args[0].ref : null;
    const Epa = +p.args[2] || 0;
    if (mid && Epa > 0) map.set(mid, Epa / 1000);            // Pa → kN/m²
  }
  return map;
}

/**
 * Analiza un .ifc y produce las fichas de todos los elementos (sin construir el modelo).
 * @param {string} text  contenido del .ifc
 * @param {object} [opts]  { tol }
 * @returns {{ schema, unit, levels, counts, items, warnings:Warnings }}
 */
export function analyzeIFC(text, opts = {}) {
  const model = parseIFC(text);
  const unit = lengthUnit(model);
  const { elements, levels, counts } = classify(model);
  const mechE = mechModulus(model);
  const levelName = new Map(levels.map(l => [l.id, l.name]));
  const W = new Warnings();

  const items = [];
  for (const el of elements) {
    const w = new Warnings();
    const item = {
      ifcId: el.id, ifcType: el.ifcType, kind: el.kind, kindLabel: KIND_LABEL[el.kind] || el.kind,
      supported: el.supported, isArea: !!el.isArea, name: el.name,
      levelName: el.storeyId != null ? (levelName.get(el.storeyId) || '—') : '—',
      segments: [], corners: null, thickness: 0, matName: '', E: null, secName: '', sec: null,
      status: 'unsupported', warnings: w,
    };

    if (!el.supported) {
      item.warnings.add(`${item.kindLabel} aún no soportado`);
      items.push(item); continue;
    }

    // material (común a barras y áreas)
    const setMaterial = () => {
      const mat = resolveMaterial(model, el.materialRef, mechE);
      item.matName = mat.name || 'Genérico (IFC)';
      if (mat.E && mat.E > 0) item.E = mat.E;
      else { item.E = GENERIC.E; w.add('Sin propiedades mecánicas: material genérico (acero)'); }
      return mat;
    };

    // ── ÁREA (muro/losa/placa): superficie de 3–4 esquinas + espesor ──
    if (el.isArea) {
      const surf = areaSurface(model, model.get(el.id), el.kind, unit.factor, w);
      if (!surf || surf.corners.length < 3) { item.status = 'no-geom'; w.add('Sin geometría de superficie reconocible: no se puede importar'); items.push(item); continue; }
      item.corners = surf.corners;
      item.thickness = surf.thickness;
      item.areaKind = surf.via;
      setMaterial();
      item.secName = `e = ${(surf.thickness * 1000).toFixed(0)} mm`;
      item.status = 'ok';
      items.push(item); continue;
    }

    // geometría → eje(s)
    const axis = memberAxis(model, model.get(el.id), unit.factor, w);
    if (!axis || !axis.segments.length) { item.status = 'no-geom'; w.add('Sin geometría de eje reconocible: no se puede importar'); items.push(item); continue; }
    item.segments = axis.segments;

    // material
    const mat = setMaterial();

    // sección: perfil del sólido extruido, o el del material, o genérica
    let prof = bodyProfile(model, model.get(el.id)) || mat.profile;
    const sp = prof ? profileProps(model, prof, unit.factor, w) : null;
    if (sp) { item.sec = { A: sp.A, Iy: sp.Iy, Iz: sp.Iz, J: sp.J }; item.secName = sp.name; if (sp.approx) item.secApprox = true; }
    else { item.secName = 'Genérica'; w.add('Sin sección reconocible: sección genérica'); }

    item.status = 'ok';
    items.push(item);
  }

  // resumen global de advertencias (agrupado)
  for (const it of items) for (const m of it.warnings.list) W.add(m);

  return { schema: model.schema, unit, levels, counts, items, warnings: W };
}

// ── snap espacial de nodos por tolerancia (hash de celdas) ────────────────────────
function makeSnapper(tol) {
  const nodes = [];                 // [ [x,y,z], … ]
  const grid = new Map();           // 'cx,cy,cz' → [idx, …]
  const cell = (p) => `${Math.round(p[0] / tol)},${Math.round(p[1] / tol)},${Math.round(p[2] / tol)}`;
  const t2 = tol * tol;
  return (p) => {
    const cx = Math.round(p[0] / tol), cy = Math.round(p[1] / tol), cz = Math.round(p[2] / tol);
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const bucket = grid.get(`${cx + dx},${cy + dy},${cz + dz}`);
      if (!bucket) continue;
      for (const i of bucket) { const q = nodes[i], d0 = q[0] - p[0], d1 = q[1] - p[1], d2 = q[2] - p[2]; if (d0 * d0 + d1 * d1 + d2 * d2 <= t2) return i; }
    }
    const idx = nodes.length; nodes.push([p[0], p[1], p[2]]);
    const k = cell(p); if (!grid.has(k)) grid.set(k, []); grid.get(k).push(idx);
    return idx;
  };
}

/**
 * Convierte las fichas SELECCIONADAS (y con estado 'ok') en un modelo neutro de `io/`.
 * @param {Array} items   fichas de `analyzeIFC`
 * @param {Set<number>} [selected]  ifcIds elegidos; si se omite, todos los 'ok'
 * @param {object} [opts]  { tol, name }
 * @returns {{ neutral, stats, warnings:string[] }}
 */
export function itemsToNeutral(items, selected = null, opts = {}) {
  const tol = opts.tol || DEFAULT_TOL;
  const snap = makeSnapper(tol);
  const warnings = [];

  const matKey = new Map(), materials = [];
  const secKey = new Map(), sections = [];
  const members = [], areas = [];
  const nodeCoords = [];   // se llena al final desde el snapper

  // acumulador de nodos: el snapper devuelve índices 0..N-1; reconstruimos coords después
  const usedNodes = new Map(); // idx → [x,y,z]
  const getNode = (p) => { const i = snap(p); if (!usedNodes.has(i)) usedNodes.set(i, p); return i; };

  let skipped = 0, skippedAreas = 0;
  for (const it of items) {
    if (it.status !== 'ok') continue;
    if (selected && !selected.has(it.ifcId)) continue;

    // material (dedupe por nombre+E) — común a barras y áreas
    const mk = `${it.matName}|${Math.round((it.E || GENERIC.E))}`;
    let mIdx = matKey.get(mk);
    if (mIdx == null) { mIdx = materials.length + 1; matKey.set(mk, mIdx); materials.push({ id: mIdx, name: it.matName, E: it.E || GENERIC.E, nu: GENERIC.nu, rho: GENERIC.rho, alpha: GENERIC.alpha }); }

    // ── ÁREA (muro/losa/placa) → cáscara de 3–4 nodos ──
    if (it.corners && it.corners.length >= 3) {
      const ids = it.corners.map(p => getNode(p) + 1);
      const uniq = []; for (const id of ids) if (!uniq.includes(id)) uniq.push(id);   // colapsa esquinas coincidentes
      if (uniq.length >= 3 && uniq.length <= 4) areas.push({ id: areas.length + 1, nodes: uniq, mat: mIdx, thickness: it.thickness || 0.2, behavior: 'shell' });
      else skippedAreas++;
      continue;
    }

    // ── BARRA → sección + miembro(s) ──
    const s = it.sec;
    const sk = s ? `${it.secName}|${s.A.toExponential(4)}|${s.Iy.toExponential(4)}|${s.Iz.toExponential(4)}` : `gen|${it.secName}`;
    let sIdx = secKey.get(sk);
    if (sIdx == null) { sIdx = sections.length + 1; secKey.set(sk, sIdx); sections.push(s ? { id: sIdx, name: it.secName, A: s.A, Iy: s.Iy, Iz: s.Iz, J: s.J } : { id: sIdx, name: it.secName }); }

    for (const [pa, pb] of it.segments) {
      if (Math.hypot(pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]) < tol) { skipped++; continue; }
      const ni = getNode(pa) + 1, nj = getNode(pb) + 1;   // ids 1..N
      if (ni === nj) { skipped++; continue; }
      members.push({ id: members.length + 1, ni, nj, mat: mIdx, sec: sIdx, releases: Array(12).fill(0), beta: 0 });
    }
  }

  // nodos en orden de índice
  const maxIdx = Math.max(-1, ...usedNodes.keys());
  for (let i = 0; i <= maxIdx; i++) {
    const c = usedNodes.get(i) || [0, 0, 0];
    nodeCoords.push({ id: i + 1, x: c[0], y: c[1], z: c[2], restraints: { ux: 0, uy: 0, uz: 0, rx: 0, ry: 0, rz: 0 }, mass: null });
  }

  if (skipped) warnings.push(`${skipped} tramo(s) de longitud ~0 descartado(s) tras el snap (tol ${tol} m)`);
  if (skippedAreas) warnings.push(`${skippedAreas} área(s) descartada(s): tras el snap no quedaron 3–4 esquinas distintas`);
  if (!members.length && !areas.length) warnings.push('No se generó ninguna barra ni área con la selección actual');

  const neutral = {
    units: { length: 'm', force: 'kN' },
    meta: { name: opts.name || 'IFC', source: 'ifc', warnings },
    nodes: nodeCoords, materials, sections, members, areas, loadCases: [],
  };
  return { neutral, stats: { nodes: nodeCoords.length, members: members.length, areas: areas.length, materials: materials.length, sections: sections.length }, warnings };
}

/** Atajo: texto IFC → `Model` de PÓRTICO con TODOS los elementos soportados (para Node/tests). */
export function ifcToModel(text, opts = {}) {
  const { items, warnings } = analyzeIFC(text, opts);
  const { neutral, stats, warnings: w2 } = itemsToNeutral(items, null, opts);
  return { model: neutralToModel(neutral), stats, warnings: [...warnings.summary(), ...w2] };
}
