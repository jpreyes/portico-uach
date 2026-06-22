// ──────────────────────────────────────────────────────────────────────────────
// tendon.js — PRETENSADO POR TENDONES (cargas equivalentes) · #60
//
// Método de las CARGAS EQUIVALENTES (load balancing, T.Y. Lin): un tendón con
// trazado curvo y fuerza P ejerce sobre el hormigón un sistema de cargas que se
// puede sustituir por cargas de gravedad/nodales ordinarias, que el solver lineal
// ya resuelve. Para un tendón PARABÓLICO con flecha `a` (sagita respecto a la
// cuerda que une las anclas) sobre una luz L:
//
//     w_eq = 8 · P · a / L²      (hacia ARRIBA si el tendón cuelga bajo el eje)
//
// más, en cada ancla, una fuerza axial P (compresión) y, si el ancla es excéntrica
// (e≠0), un momento P·e.  El momento PRIMARIO en cualquier sección es M(x)=P·e(x).
//
// Trazado POLIGONAL: en cada punto de quiebre el tendón aplica una carga puntual
// transversal igual a P·(cambio de pendiente).
//
// Pérdidas: modelo de fricción/ondulación  P(x) = P0·e^(−(μ·θ + k·x))  (θ = cambio
// angular acumulado) y/o una fracción global a tanto alzado.  La carga equivalente
// usa la fuerza EFECTIVA media a lo largo del tendón.
//
// Convención del modelo: la viga se asume aproximadamente horizontal (eje global X)
// con la flexión en el plano X–Z; `e` (excentricidad) es POSITIVA hacia ABAJO (−Z),
// el caso típico de un tendón que cuelga.  Carga `gravity` w>0 = ↓, w<0 = ↑.
// ──────────────────────────────────────────────────────────────────────────────

// Longitud de un elemento.
function elemLen(model, el) {
  const a = model.nodes.get(el.n1), b = model.nodes.get(el.n2);
  return Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z);
}

// Fuerza efectiva del tendón con pérdidas por fricción/ondulación + tanto alzado.
//   tendon.jack  = P0 (kN) en el gato (extremo activo)
//   tendon.P     = fuerza efectiva directa (si se da, se usa tal cual)
//   tendon.friction = { mu, k }  (μ por radián, k por metro)
//   tendon.lumpSum  = fracción de pérdida global adicional (0..1) (retracción,
//                     fluencia, relajación a largo plazo, acortamiento elástico)
//   a = sagita, L = luz  (para estimar el cambio angular total de la parábola)
// Devuelve { P0, Pavg, Pend } (Pend = en el extremo pasivo, la más baja).
export function tendonForce(tendon, L, a) {
  if (tendon.P != null && tendon.jack == null) {
    const P = +tendon.P; return { P0: P, Pavg: P, Pend: P };
  }
  const P0 = +tendon.jack || 0;
  const fr = tendon.friction || {};
  const mu = +fr.mu || 0, k = +fr.k || 0;
  // Cambio angular total de una parábola de sagita a y luz L: pendiente de extremo
  // = 4a/L; de un extremo al otro el giro acumulado ≈ 8a/L.
  const thetaTot = L > 0 ? 8 * Math.abs(a) / L : 0;
  const Pend = P0 * Math.exp(-(mu * thetaTot + k * L));
  // P media a lo largo (integral de la exponencial ≈ promedio extremos para β chico).
  let Pavg = (P0 + Pend) / 2;
  const lump = Math.min(Math.max(+tendon.lumpSum || 0, 0), 1);
  Pavg *= (1 - lump); const PendEff = Pend * (1 - lump);
  return { P0, Pavg, Pend: PendEff };
}

// Eccentricidad del trazado en la fracción s∈[0,1] de la luz (positiva ↓).
//   parábola: e(s) = e1 + (e2−e1)·s + 4·a·s(1−s),  a = em − (e1+e2)/2  (sagita)
export function tendonEcc(tendon, s) {
  const e1 = +tendon.e?.start || 0, e2 = +tendon.e?.end || 0, em = +tendon.e?.mid || 0;
  const a = em - (e1 + e2) / 2;
  return e1 + (e2 - e1) * s + 4 * a * s * (1 - s);
}

/**
 * Cargas equivalentes de un tendón → arreglo de cargas de modelo
 * ({type:'dist'…} y {type:'nodal'…}) listas para añadir a un caso de carga.
 *
 * @param {Model}  model
 * @param {object} tendon
 *    elems   : [elemId…] elementos COLINEALES que forman la viga (en orden).
 *    profile : 'parabola' (def.) | 'polygon'
 *    e       : { start, mid, end }  excentricidades ↓+ (m) — parábola
 *    points  : [{s,e}]              quiebres (s∈[0,1], e ↓+) — poligonal
 *    P | jack/friction/lumpSum      fuerza efectiva o de tesado + pérdidas
 * @returns {{loads: Array, P: number, weq: number, L: number}}
 */
export function tendonEquivalentLoads(model, tendon) {
  const elems = (tendon.elems || []).map(id => model.elements.get(id)).filter(Boolean);
  if (!elems.length) throw new Error('tendón sin elementos válidos');

  // Luz total y posición acumulada (s) de cada nodo a lo largo del tendón.
  const lens = elems.map(e => elemLen(model, e));
  const L = lens.reduce((x, y) => x + y, 0);
  if (!(L > 0)) throw new Error('tendón de luz nula');

  // Nodos extremos y dirección axial del primer/último elemento (unitaria).
  const first = elems[0], last = elems[elems.length - 1];
  const nA = model.nodes.get(first.n1), nB = model.nodes.get(last.n2);
  const axA = unit(nA, model.nodes.get(first.n2));   // del nodo A hacia el interior
  const axB = unit(nB, model.nodes.get(last.n1));    // del nodo B hacia el interior

  const e1 = +tendon.e?.start || 0, e2 = +tendon.e?.end || 0, em = +tendon.e?.mid || 0;
  const a  = em - (e1 + e2) / 2;                       // sagita respecto a la cuerda
  const { Pavg: P } = tendonForce(tendon, L, a);

  const loads = [];
  let weq = 0;

  if ((tendon.profile || 'parabola') === 'parabola') {
    // Carga equivalente uniforme (↑ si a>0): w = 8 P a / L².  En el modelo,
    // gravity w>0 = ↓, así que un empuje hacia arriba es w_model = −8Pa/L².
    weq = 8 * P * a / (L * L);                         // magnitud física (↑+)
    for (const e of elems) loads.push({ type: 'dist', elemId: e.id, dir: 'gravity', w: -weq });
  } else {
    // Poligonal: carga puntual en cada quiebre = P·(Δpendiente). Se reparte al
    // elemento que contiene el punto por funciones de forma (consistente).
    const pts = (tendon.points || []).slice().sort((p, q) => p.s - q.s);
    for (let i = 1; i < pts.length - 1; i++) {
      const sPrev = pts[i - 1].s, sCur = pts[i].s, sNext = pts[i + 1].s;
      const slopeIn  = (pts[i].e - pts[i - 1].e) / ((sCur - sPrev) * L);
      const slopeOut = (pts[i + 1].e - pts[i].e) / ((sNext - sCur) * L);
      const dSlope = slopeOut - slopeIn;              // cambio de pendiente (↓+)
      const Fup = P * dSlope;                          // fuerza ↑ (si el quiebre cuelga)
      addTransversePoint(model, elems, lens, L, sCur * L, -Fup, loads);
    }
  }

  // Anclas: axial P (compresión) + momento P·e en cada extremo (si excéntrico).
  // Fuerza axial hacia el interior del miembro en cada ancla → compresión P.
  loads.push({ type: 'nodal', nodeId: nA.id, F: [P * axA[0], P * axA[1], P * axA[2], 0, 0, 0] });
  loads.push({ type: 'nodal', nodeId: nB.id, F: [P * axB[0], P * axB[1], P * axB[2], 0, 0, 0] });
  // Momento primario de ancla M = P·e (flexión en plano X–Z → momento global My).
  // e ↓+: una compresión P aplicada bajo el eje genera momento que tracciona la
  // fibra inferior (sagging). Con la convención My del modelo el momento de ancla
  // que reproduce M(x)=P·e en el extremo es −P·e en el nodo inicial y +P·e en el final.
  if (Math.abs(e1) > 1e-12) addEndMoment(loads, nA.id, -P * e1);
  if (Math.abs(e2) > 1e-12) addEndMoment(loads, nB.id, +P * e2);

  return { loads, P, weq, L };
}

// Conveniencia: aplica las cargas del tendón a un caso de carga del modelo.
export function applyTendon(model, lcId, tendon) {
  const { loads } = tendonEquivalentLoads(model, tendon);
  const lc = model.loadCases.get(lcId);
  if (!lc) throw new Error('caso de carga inexistente: ' + lcId);
  for (const ld of loads) lc.loads.push(ld);
  return loads;
}

// ── helpers ─────────────────────────────────────────────────────────────────
function unit(from, to) {
  const dx = to.x - from.x, dy = to.y - from.y, dz = to.z - from.z;
  const L = Math.hypot(dx, dy, dz) || 1; return [dx / L, dy / L, dz / L];
}

function addEndMoment(loads, nodeId, My) {
  loads.push({ type: 'nodal', nodeId, F: [0, 0, 0, 0, My, 0] });
}

// Reparte una carga transversal puntual (vertical global Z, signo en `Fz`) ubicada
// a distancia `x` del inicio del tendón, al elemento que la contiene, por funciones
// de forma cúbicas de Hermite (fuerzas y momentos nodales consistentes).
function addTransversePoint(model, elems, lens, L, x, Fz, loads) {
  let acc = 0, idx = 0;
  for (; idx < elems.length; idx++) { if (x <= acc + lens[idx] + 1e-9) break; acc += lens[idx]; }
  idx = Math.min(idx, elems.length - 1);
  const el = elems[idx], Le = lens[idx];
  const xi = Math.min(Math.max((x - acc) / Le, 0), 1);
  const n1 = model.nodes.get(el.n1), n2 = model.nodes.get(el.n2);
  // Funciones de forma de viga (carga ↑ vertical global): reparto consistente.
  const N1 = 1 - 3 * xi * xi + 2 * xi ** 3, N2 = 3 * xi * xi - 2 * xi ** 3;
  const M1 =  Le * (xi - 2 * xi * xi + xi ** 3), M2 = Le * (-xi * xi + xi ** 3);
  // Fz vertical global (Z); momento de flexión asociado → My (plano X–Z).
  loads.push({ type: 'nodal', nodeId: n1.id, F: [0, 0, Fz * N1, 0, -Fz * M1, 0] });
  loads.push({ type: 'nodal', nodeId: n2.id, F: [0, 0, Fz * N2, 0, -Fz * M2, 0] });
}
