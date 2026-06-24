// ──────────────────────────────────────────────────────────────────────────────
// macro_registry.js — REGISTRO de MACROMODELOS (#86)
//
// Un macromodelo es un subsistema estructural complejo que se resuelve con POCOS
// elementos calibrados (barras/cables/resortes/links no lineales) en vez de un mallado
// fino — p.ej. un muro de relleno → puntal diagonal equivalente.  Este registro hace que
// AGREGAR un macromodelo nuevo sea registrar un descriptor `{ id, name, nodes, params,
// expand }`; no hay que tocar el solver ni (con la UI genérica) escribir diálogo nuevo.
//
// El autor desarrolla la TEORÍA del macromodelo aparte (calibración del puntal/resorte/
// histéresis) y aquí sólo se «enchufa»: `expand(model, nodeIds, props)` construye en el
// modelo la red de elementos ya calibrada y los marca con `el.macro`/`el.macroType`.
//
//   registerMacro(def)                    → registra/sobre-escribe un macromodelo
//   getMacro(id) / listMacros()           → consulta (la UI lee esto para el menú)
//   insertMacro(model, id, nodeIds, props)→ ejecuta el expand; { error } | { …, macroId }
// ──────────────────────────────────────────────────────────────────────────────

const _macros = new Map();

/**
 * Registra un macromodelo conectable.
 * @param {object} def
 *   id        {string}   identificador único ('infill', 'shearwall', 'bracing', …)
 *   name      {string}   nombre legible (menú/diálogo)
 *   desc      {string}   descripción corta (1 línea) — referencia teórica
 *   nodes     {number}   nº de nodos que el usuario selecciona (p.ej. 4 esquinas)
 *   nodesHint {string}   texto de ayuda sobre qué nodos seleccionar
 *   dims      {'2D'|'3D'|null}  restricción de modo (null = ambos)
 *   params    {Array}    descriptores de parámetros para auto-generar el diálogo:
 *                        [{ key, label, default, step?, min? }]
 *   expand    {(model, nodeIds:number[], props:object) => ({error}|{macroId,...})}
 *             construye la red interna en el modelo (elementos calibrados); debe
 *             marcar cada elemento creado con `el.macro=<id num>` y `el.macroType=def.id`,
 *             y registrar el macro en `model.macros` (ver `insertInfill` de ejemplo).
 */
export function registerMacro(def) {
  if (!def || !def.id) throw new Error('registerMacro: falta id');
  if (typeof def.expand !== 'function') throw new Error(`registerMacro «${def.id}»: falta expand()`);
  _macros.set(def.id, {
    id: def.id, name: def.name || def.id, desc: def.desc || '',
    nodes: def.nodes || 0, nodesHint: def.nodesHint || '', dims: def.dims || null,
    params: def.params || [], expand: def.expand,
  });
  return def.id;
}

export function getMacro(id) { return _macros.get(id) || null; }

/** Descriptores de los macromodelos registrados (sin las funciones) para la UI. */
export function listMacros() {
  return [..._macros.values()].map(({ id, name, desc, nodes, nodesHint, dims, params }) => ({ id, name, desc, nodes, nodesHint, dims, params }));
}

/**
 * Inserta un macromodelo en el modelo.
 * @returns {{error:string} | {macroId:number, ...}}
 */
export function insertMacro(model, id, nodeIds, props = {}) {
  const def = getMacro(id);
  if (!def) return { error: `Macromodelo desconocido: ${id}` };
  if (def.dims && model.mode !== def.dims) return { error: `«${def.name}» requiere un modelo ${def.dims}.` };
  if (def.nodes && (!Array.isArray(nodeIds) || nodeIds.length !== def.nodes)) return { error: `Seleccione ${def.nodes} nodo(s): ${def.nodesHint || ''}` };
  return def.expand(model, nodeIds, props);
}
