// ──────────────────────────────────────────────────────────────────────────────
// links.js — LINKS / COUPLINGS entre nodos (restricciones cinemáticas) · puentes
//
// Permiten ligar dos nodos que ocupan posiciones distintas, sin un elemento entre
// ellos.  El caso típico de puentes: el TABLERO se modela en su eje (más arriba) y
// la VIGA en el suyo; un **link rígido** transmite fuerzas y momentos entre ambos
// respetando el brazo (offset) → la fuerza del tablero llega a la viga como fuerza
// + momento = F·brazo.  También sirve para apoyos excéntricos, vigas con offset de
// extremo (end offsets, 1-010), insertion points (1-011) y ejes de cálculo.
//
// Dos tipos:
//   · rigid = true  → LINK RÍGIDO: el esclavo sigue al maestro como sólido rígido
//        u_s = u_m + θ_m × r,   θ_s = θ_m      (r = posición_esclavo − maestro)
//     Transmite los 6 GDL con el brazo. `dofs` puede acotar cuáles GDL del esclavo
//     se ligan (p.ej. sólo traslaciones → rótula que transmite fuerza, no momento).
//   · rigid = false → COUPLING simple: iguala los GDL seleccionados sin brazo
//        gdl_s = gdl_m    (para los `dofs` marcados)
//
// Implementación por PENALIZACIÓN (idéntica a los diafragmas): para cada ecuación
// de restricción g·u = 0 se suma α·gᵀg a K.  α = max(diag K)·1e5.  AUTÓNOMO.
// ──────────────────────────────────────────────────────────────────────────────

const PENALTY_FACTOR = 1e5;   // mismo factor que los diafragmas (error <0.001%)

function denseWriter(K, nDOF) {
  return { add: (i, j, v) => { K[i * nDOF + j] += v; }, diag: (i) => K[i * nDOF + i] };
}

function _addPenalty(W, alpha, dofs, coeffs) {
  for (let i = 0; i < dofs.length; i++)
    for (let j = 0; j < dofs.length; j++)
      W.add(dofs[i], dofs[j], alpha * coeffs[i] * coeffs[j]);
}

// Ecuaciones de restricción de un link → lista de {dofs:[gi…], coeffs:[…]}.
// gdl(im) = 6·im + {0..5} = [ux,uy,uz, rx,ry,rz].
function _linkEquations(link, master, slave, im, is) {
  const M = k => 6 * im + k, S = k => 6 * is + k;
  const d = link.dofs || { ux: 1, uy: 1, uz: 1, rx: 1, ry: 1, rz: 1 };
  const eqs = [];
  if (link.rigid) {
    const dx = slave.x - master.x, dy = slave.y - master.y, dz = slave.z - master.z;
    // u_s = u_m + θ_m × r   (θ×r)=(ry·dz−rz·dy, rz·dx−rx·dz, rx·dy−ry·dx)
    if (d.ux) eqs.push({ dofs: [S(0), M(0), M(4), M(5)], coeffs: [1, -1, -dz, dy] });
    if (d.uy) eqs.push({ dofs: [S(1), M(1), M(5), M(3)], coeffs: [1, -1, -dx, dz] });
    if (d.uz) eqs.push({ dofs: [S(2), M(2), M(3), M(4)], coeffs: [1, -1, -dy, dx] });
    if (d.rx) eqs.push({ dofs: [S(3), M(3)], coeffs: [1, -1] });
    if (d.ry) eqs.push({ dofs: [S(4), M(4)], coeffs: [1, -1] });
    if (d.rz) eqs.push({ dofs: [S(5), M(5)], coeffs: [1, -1] });
  } else {
    // coupling simple: gdl_s = gdl_m para los GDL marcados
    ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'].forEach((k, idx) => { if (d[k]) eqs.push({ dofs: [S(idx), M(idx)], coeffs: [1, -1] }); });
  }
  return eqs;
}

// ── API densa (usada por assembleK) ─────────────────────────────────────────────
export function applyLinkConstraints(K, model, nodeIndex, nDOF) {
  applyLinkConstraintsW(denseWriter(K, nDOF), model, nodeIndex, nDOF);
}

// Variante con writer (denso o disperso), para el camino en banda (sparse.js).
export function applyLinkConstraintsW(W, model, nodeIndex, nDOF) {
  if (!model.links || model.links.size === 0) return;
  let maxKii = 0;
  for (let i = 0; i < nDOF; i++) { const v = W.diag(i); if (v > maxKii) maxKii = v; }
  const alpha = maxKii > 0 ? maxKii * PENALTY_FACTOR : 1e12;

  for (const link of model.links.values()) {
    const master = model.nodes.get(link.master), slave = model.nodes.get(link.slave);
    if (!master || !slave || link.master === link.slave) continue;
    const im = nodeIndex.get(link.master), is = nodeIndex.get(link.slave);
    if (im == null || is == null) continue;
    for (const eq of _linkEquations(link, master, slave, im, is)) _addPenalty(W, alpha, eq.dofs, eq.coeffs);
  }
}
