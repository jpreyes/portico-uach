// ──────────────────────────────────────────────────────────────────────────────
// io/ifc/ifcClassifier.js — clasificación de elementos IFC · #75, G19
//
// Recorre un `IfcModel` y clasifica los elementos constructivos según lo que PÓRTICO
// sabe importar HOY (sólo BARRAS):
//   • barras  → IfcBeam / IfcColumn / IfcMember (y sus *StandardCase)   ── soportado
//   • muros   → IfcWall                                                 ── no soportado aún
//   • losas   → IfcSlab                                                 ── no soportado aún
//   • placas  → IfcPlate                                               ── no soportado aún
// (Los muros/losas/placas se LISTAN como «no soportado» para que el usuario vea qué hay
//  en el archivo; se activarán cuando exista importación real de áreas.)
//
// Además resuelve las RELACIONES útiles: a qué NIVEL pertenece cada elemento
// (IfcRelContainedInSpatialStructure) y qué MATERIAL/PERFIL tiene asociado
// (IfcRelAssociatesMaterial).  AUTÓNOMO (Node + navegador).
// ──────────────────────────────────────────────────────────────────────────────

// tipo IFC → clase de PÓRTICO + si es importable hoy
const KIND = new Map([
  ['IFCBEAM', 'beam'], ['IFCBEAMSTANDARDCASE', 'beam'],
  ['IFCCOLUMN', 'column'], ['IFCCOLUMNSTANDARDCASE', 'column'],
  ['IFCMEMBER', 'member'], ['IFCMEMBERSTANDARDCASE', 'member'],
  ['IFCWALL', 'wall'], ['IFCWALLSTANDARDCASE', 'wall'], ['IFCWALLELEMENTEDCASE', 'wall'],
  ['IFCSLAB', 'slab'], ['IFCSLABSTANDARDCASE', 'slab'], ['IFCSLABELEMENTEDCASE', 'slab'],
  ['IFCPLATE', 'plate'], ['IFCPLATESTANDARDCASE', 'plate'],
  ['IFCFOOTING', 'footing'], ['IFCPILE', 'pile'],
]);
// barras → IfcBeam/Column/Member; áreas → IfcWall/Slab/Plate (3–4 nodos en PÓRTICO)
const SUPPORTED = new Set(['beam', 'column', 'member', 'wall', 'slab', 'plate']);
const AREA_KINDS = new Set(['wall', 'slab', 'plate']);

/** Etiqueta legible (es) por clase, para la UI. */
export const KIND_LABEL = { beam: 'Viga', column: 'Pilar', member: 'Barra', wall: 'Muro', slab: 'Losa', plate: 'Placa', footing: 'Zapata', pile: 'Pilote' };

/**
 * Clasifica los elementos del `IfcModel`.
 * @returns {{
 *   elements: Array<{id,ifcType,kind,supported,name,predefined,storeyId,materialRef}>,
 *   levels:   Array<{id,name,elevation}>,
 *   counts:   Record<string,number>
 * }}
 */
export function classify(model) {
  // ── niveles (IfcBuildingStorey) ──
  const levels = model.ofType('IFCBUILDINGSTOREY').map(s => ({
    id: s.id,
    name: (s.args[2] || s.args[7] || `Nivel ${s.id}`).toString(),
    elevation: +s.args[9] || 0,
  })).sort((a, b) => a.elevation - b.elevation);

  // ── contención espacial: elemento → nivel ──
  const storeyOf = new Map();
  for (const rel of model.ofType('IFCRELCONTAINEDINSPATIALSTRUCTURE')) {
    const structure = rel.args[5];                 // RelatingStructure
    const sid = model.isRef(structure) ? structure.ref : null;
    for (const o of (rel.args[4] || [])) if (model.isRef(o)) storeyOf.set(o.ref, sid);
  }

  // ── asociación de material/perfil: elemento → RelatingMaterial ──
  const matOf = new Map();
  for (const rel of model.ofType('IFCRELASSOCIATESMATERIAL')) {
    const mat = rel.args[5];                        // RelatingMaterial
    for (const o of (rel.args[4] || [])) if (model.isRef(o)) matOf.set(o.ref, mat);
  }

  // ── elementos ──
  const elements = [];
  const counts = {};
  for (const [type, kind] of KIND) {
    for (const e of model.ofType(type)) {
      counts[kind] = (counts[kind] || 0) + 1;
      elements.push({
        id: e.id, ifcType: type, kind, supported: SUPPORTED.has(kind), isArea: AREA_KINDS.has(kind),
        name: (e.args[2] || `${KIND_LABEL[kind] || kind} ${e.id}`).toString(),
        predefined: (typeof e.args[8] === 'string' ? e.args[8] : '') || '',
        storeyId: storeyOf.get(e.id) ?? null,
        materialRef: matOf.get(e.id) ?? null,
      });
    }
  }

  return { elements, levels, counts };
}

export { SUPPORTED, AREA_KINDS, KIND };
