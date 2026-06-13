// ──────────────────────────────────────────────────────────────────────────────
// Discretización de elementos — dividir y unir
//
// splitElement(model, elemId, nParts)   → divide un elemento en N sub-elementos
// splitByLength(model, elemId, len)     → divide cada ≈len metros
// discretizeAll(model, opts)            → divide todos los elementos
// joinElements(model, elemIds)          → une una cadena colineal en 1 elemento
//
// Reglas al dividir:
//  - Las liberaciones del extremo 1 quedan en el primer sub-elemento;
//    las del extremo 2 en el último (los nodos intermedios son continuos).
//  - Las cargas distribuidas (UDL) se replican en cada sub-elemento,
//    en TODOS los casos de carga.
// ──────────────────────────────────────────────────────────────────────────────

const EPS = 1e-9;

/** Divide un elemento en nParts sub-elementos. Devuelve los IDs nuevos (o null). */
export function splitElement(model, elemId, nParts) {
  nParts = Math.max(2, Math.floor(nParts));
  const elem = model.elements.get(elemId);
  if (!elem) return null;
  const n1 = model.nodes.get(elem.n1);
  const n2 = model.nodes.get(elem.n2);
  if (!n1 || !n2) return null;

  // Nodos intermedios
  const chain = [elem.n1];
  for (let k = 1; k < nParts; k++) {
    const t = k / nParts;
    const nd = model.addNode(
      n1.x + t * (n2.x - n1.x),
      n1.y + t * (n2.y - n1.y),
      n1.z + t * (n2.z - n1.z)
    );
    chain.push(nd.id);
  }
  chain.push(elem.n2);

  // Sub-elementos
  const newIds = [];
  for (let k = 0; k < nParts; k++) {
    const e = model.addElement(chain[k], chain[k + 1], elem.matId, elem.secId);
    if (!e) continue;
    newIds.push(e.id);
  }
  // Liberaciones: extremo 1 → primer sub-elemento; extremo 2 → último
  if (elem.releases && newIds.length) {
    const first = model.elements.get(newIds[0]);
    const last  = model.elements.get(newIds[newIds.length - 1]);
    for (let i = 0; i < 6; i++) {
      first.releases[i]    = elem.releases[i]     || 0;
      last.releases[6 + i] = elem.releases[6 + i] || 0;
    }
  }

  // Cargas distribuidas: replicar en cada sub-elemento (todos los casos)
  for (const lc of model.loadCases.values()) {
    const distLoads = lc.loads.filter(l => l.type === 'dist' && l.elemId === elemId);
    if (!distLoads.length) continue;
    lc.loads = lc.loads.filter(l => !(l.type === 'dist' && l.elemId === elemId));
    for (const dl of distLoads) {
      for (const id of newIds) {
        lc.loads.push({ ...dl, elemId: id });
      }
    }
  }

  model.elements.delete(elemId);
  return newIds;
}

/** Divide un elemento en tramos de ≈targetLen metros. */
export function splitByLength(model, elemId, targetLen) {
  const elem = model.elements.get(elemId);
  if (!elem || !(targetLen > EPS)) return null;
  const n1 = model.nodes.get(elem.n1);
  const n2 = model.nodes.get(elem.n2);
  if (!n1 || !n2) return null;
  const L = Math.hypot(n2.x - n1.x, n2.y - n1.y, n2.z - n1.z);
  const nParts = Math.max(1, Math.round(L / targetLen));
  if (nParts < 2) return [elemId];   // ya es más corto que el objetivo
  return splitElement(model, elemId, nParts);
}

/**
 * Discretiza todos los elementos del modelo.
 * opts: { parts: N }  ó  { length: metros }
 * Devuelve el nº de elementos resultantes.
 */
export function discretizeAll(model, opts = {}) {
  const ids = [...model.elements.keys()];
  for (const id of ids) {
    if (opts.length > EPS)      splitByLength(model, id, opts.length);
    else if (opts.parts >= 2)   splitElement(model, id, opts.parts);
  }
  return model.elements.size;
}

/**
 * Une una cadena de elementos colineales en un solo elemento.
 * Requisitos: mismo material y sección, cadena continua y colineal,
 * nodos intermedios sin apoyos, cargas nodales, masas, diafragmas
 * ni otros elementos conectados, y sin rótulas intermedias.
 * Devuelve { ok, reason?, elemId?, removedNodes? }.
 */
export function joinElements(model, elemIds) {
  const elems = elemIds.map(id => model.elements.get(id)).filter(Boolean);
  if (elems.length < 2) return { ok: false, reason: 'Seleccione al menos 2 elementos' };

  const ref = elems[0];
  if (!elems.every(e => e.matId === ref.matId && e.secId === ref.secId)) {
    return { ok: false, reason: 'Los elementos deben tener el mismo material y sección' };
  }

  // Contar ocurrencias de nodos: extremos aparecen 1 vez, intermedios 2
  const count = new Map();
  for (const e of elems) {
    count.set(e.n1, (count.get(e.n1) || 0) + 1);
    count.set(e.n2, (count.get(e.n2) || 0) + 1);
  }
  const endNodes = [...count.entries()].filter(([, c]) => c === 1).map(([id]) => id);
  const midNodes = [...count.entries()].filter(([, c]) => c === 2).map(([id]) => id);
  if (endNodes.length !== 2 || [...count.values()].some(c => c > 2)) {
    return { ok: false, reason: 'Los elementos no forman una cadena simple' };
  }

  // Ordenar la cadena desde un extremo
  const [startId, endId] = endNodes;
  const remaining = new Set(elems.map(e => e.id));
  const ordered = [];   // [{elem, flip}] — flip=true si el elemento apunta contra la cadena
  let cur = startId;
  while (remaining.size) {
    let found = null;
    for (const id of remaining) {
      const e = model.elements.get(id);
      if (e.n1 === cur)      { found = { elem: e, flip: false }; cur = e.n2; break; }
      else if (e.n2 === cur) { found = { elem: e, flip: true  }; cur = e.n1; break; }
    }
    if (!found) return { ok: false, reason: 'Cadena discontinua' };
    remaining.delete(found.elem.id);
    ordered.push(found);
  }
  if (cur !== endId) return { ok: false, reason: 'Cadena discontinua' };

  // Colinealidad: todos los nodos sobre la recta start→end
  const nS = model.nodes.get(startId), nE = model.nodes.get(endId);
  const d  = [nE.x - nS.x, nE.y - nS.y, nE.z - nS.z];
  const Lt = Math.hypot(...d);
  if (Lt < EPS) return { ok: false, reason: 'Longitud nula' };
  for (const mid of midNodes) {
    const n = model.nodes.get(mid);
    const v = [n.x - nS.x, n.y - nS.y, n.z - nS.z];
    const cx = v[1] * d[2] - v[2] * d[1];
    const cy = v[2] * d[0] - v[0] * d[2];
    const cz = v[0] * d[1] - v[1] * d[0];
    if (Math.hypot(cx, cy, cz) / Lt > 1e-6) {
      return { ok: false, reason: 'Los elementos no son colineales' };
    }
  }

  // Nodos intermedios deben estar "limpios"
  for (const mid of midNodes) {
    const n = model.nodes.get(mid);
    if (Object.values(n.restraints || {}).some(v => v)) {
      return { ok: false, reason: `Nodo ${mid} tiene apoyos` };
    }
    const nm = n.nodeMass;
    if (nm && (nm.mx || nm.my || nm.mz)) {
      return { ok: false, reason: `Nodo ${mid} tiene masa nodal` };
    }
    for (const lc of model.loadCases.values()) {
      if (lc.loads.some(l => l.type === 'nodal' && l.nodeId === mid)) {
        return { ok: false, reason: `Nodo ${mid} tiene cargas nodales` };
      }
    }
    for (const dia of model.diaphragms.values()) {
      if (dia.nodes.includes(mid)) {
        return { ok: false, reason: `Nodo ${mid} pertenece a un diafragma` };
      }
    }
    for (const e of model.elements.values()) {
      if (elemIds.includes(e.id)) continue;
      if (e.n1 === mid || e.n2 === mid) {
        return { ok: false, reason: `Nodo ${mid} conecta otros elementos` };
      }
    }
  }

  // Sin rótulas en extremos intermedios (sería una rótula real dentro del nuevo elemento)
  const endRel = (e, atN1) => {
    const off = atN1 ? 0 : 6;
    return (e.releases || []).slice(off, off + 6);
  };
  for (let k = 0; k < ordered.length; k++) {
    const { elem, flip } = ordered[k];
    const innerStart = k > 0;                      // extremo que toca nodo intermedio
    const innerEnd   = k < ordered.length - 1;
    const relStart = endRel(elem, !flip);          // extremo en el lado "inicio de cadena"
    const relEnd   = endRel(elem, flip);
    if ((innerStart && relStart.some(r => r)) || (innerEnd && relEnd.some(r => r))) {
      return { ok: false, reason: 'Hay liberaciones (rótulas) en nodos intermedios' };
    }
  }

  // Cargas distribuidas: unir solo si todos los tramos tienen la MISMA carga por caso
  const mergedDist = [];   // {lcId, dir, w}
  for (const lc of model.loadCases.values()) {
    const sig = e => {
      const dl = lc.loads.filter(l => l.type === 'dist' && l.elemId === e.id);
      if (dl.length === 0) return 'none';
      if (dl.length > 1) return 'multi';
      return `${dl[0].dir || 'gravity'}|${dl[0].w}`;
    };
    const sigs = elems.map(sig);
    if (!sigs.every(s => s === sigs[0])) {
      return { ok: false, reason: `Cargas distribuidas distintas entre tramos (caso "${lc.name}")` };
    }
    if (sigs[0] !== 'none' && sigs[0] !== 'multi') {
      const [dir, w] = sigs[0].split('|');
      mergedDist.push({ lcId: lc.id, dir, w: +w });
    }
  }

  // Liberaciones del elemento unido: extremos exteriores de la cadena
  const relA = endRel(ordered[0].elem, !ordered[0].flip);
  const last = ordered[ordered.length - 1];
  const relB = endRel(last.elem, last.flip);

  // ── Ejecutar la unión ──
  for (const { elem } of ordered) {
    for (const lc of model.loadCases.values()) {
      lc.loads = lc.loads.filter(l => !(l.type === 'dist' && l.elemId === elem.id));
    }
    model.elements.delete(elem.id);
  }
  for (const mid of midNodes) model.nodes.delete(mid);

  const merged = model.addElement(startId, endId, ref.matId, ref.secId);
  for (let i = 0; i < 6; i++) {
    merged.releases[i]     = relA[i] ? 1 : 0;
    merged.releases[6 + i] = relB[i] ? 1 : 0;
  }
  for (const md of mergedDist) {
    model.loadCases.get(md.lcId)?.loads.push({ type: 'dist', elemId: merged.id, dir: md.dir, w: md.w });
  }

  return { ok: true, elemId: merged.id, removedNodes: midNodes };
}
