// ──────────────────────────────────────────────────────────────────────────────
// polygon_props.js — Propiedades de una sección POLIGONAL arbitraria (#70).
//
// Para el «Section Designer»: dado un contorno (polígono) y, opcionalmente, huecos
// (polígonos interiores), calcula por el teorema de Green TODAS las propiedades de
// sección: área, centroide (y,z), momentos de inercia Iz, Iy y el PRODUCTO Iyz,
// ejes principales (I₁, I₂, θ), módulos elásticos Sz/Sy, módulos PLÁSTICOS Zz/Zy
// (eje neutro de áreas iguales por bisección + recorte de semiplano), perímetro y
// envolvente. Soporta secciones cóncavas y con huecos.
//
// Convención de ejes del proyecto: z = eje fuerte (horizontal), Iz = ∫y²dA; y = eje
// débil (vertical), Iy = ∫x²dA. Coordenadas del polígono: x = horizontal, y = vertical.
// Unidades: las del input (m → m², m⁴).
// ──────────────────────────────────────────────────────────────────────────────

// Integrales de un lazo sobre el ORIGEN (CCW → positivas). Normaliza a área positiva.
function loopIntegrals(loop) {
  let A = 0, Qx = 0, Qy = 0, Ixx = 0, Iyy = 0, Ixy = 0, per = 0;
  const n = loop.length;
  for (let i = 0; i < n; i++) {
    const [x0, y0] = loop[i], [x1, y1] = loop[(i + 1) % n];
    const cr = x0 * y1 - x1 * y0;
    A += cr;
    Qx += (x0 + x1) * cr;                          // 6·∫x dA
    Qy += (y0 + y1) * cr;                          // 6·∫y dA
    Iyy += (x0 * x0 + x0 * x1 + x1 * x1) * cr;     // 12·∫x² dA
    Ixx += (y0 * y0 + y0 * y1 + y1 * y1) * cr;     // 12·∫y² dA
    Ixy += (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * cr;   // 24·∫xy dA
    per += Math.hypot(x1 - x0, y1 - y0);
  }
  let r = { A: A / 2, Qx: Qx / 6, Qy: Qy / 6, Iyy: Iyy / 12, Ixx: Ixx / 12, Ixy: Ixy / 24, per };
  if (r.A < 0) { r.A = -r.A; r.Qx = -r.Qx; r.Qy = -r.Qy; r.Iyy = -r.Iyy; r.Ixx = -r.Ixx; r.Ixy = -r.Ixy; }
  return r;
}

// Recorta un lazo al semiplano coord(axis) ≥ c (Sutherland–Hodgman). axis: 'x'|'y'.
function clipGE(loop, axis, c) {
  const val = p => axis === 'y' ? p[1] : p[0];
  const out = []; const n = loop.length;
  for (let i = 0; i < n; i++) {
    const a = loop[i], b = loop[(i + 1) % n];
    const va = val(a), vb = val(b), inA = va >= c, inB = vb >= c;
    if (inA) out.push(a);
    if (inA !== inB) { const t = (c - va) / (vb - va); out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]); }
  }
  return out.length >= 3 ? out : null;
}

// Combina contorno (+) y huecos (−) para una integral; clipAxis/clipC opcionales.
function combine(outline, holes, field, clipAxis, clipC) {
  const loops = [{ loop: outline, s: 1 }, ...holes.map(h => ({ loop: h, s: -1 }))];
  let v = 0;
  for (const { loop, s } of loops) {
    const L = (clipAxis ? clipGE(loop, clipAxis, clipC) : loop);
    if (!L) continue;
    v += s * loopIntegrals(L)[field];
  }
  return v;
}

/**
 * @param {object} o  { outline:[[x,y]…], holes?:[[[x,y]…]…] }
 * @returns propiedades de sección (ver cabecera).
 */
export function polygonProps({ outline, holes = [] }) {
  if (!outline || outline.length < 3) throw new Error('contorno poligonal inválido (≥3 vértices)');
  const loops = [{ loop: outline, s: 1 }, ...holes.map(h => ({ loop: h, s: -1 }))];
  let A = 0, Qx = 0, Qy = 0, Ixx0 = 0, Iyy0 = 0, Ixy0 = 0, per = 0;
  for (const { loop, s } of loops) {
    const g = loopIntegrals(loop);
    A += s * g.A; Qx += s * g.Qx; Qy += s * g.Qy; Ixx0 += s * g.Ixx; Iyy0 += s * g.Iyy; Ixy0 += s * g.Ixy;
    per += g.per;
  }
  if (!(A > 1e-12)) throw new Error('área de la sección nula o negativa (¿orden de vértices?)');
  const cx = Qx / A, cy = Qy / A;
  // Centroidales (proyecto: Iz=∫y²dA, Iy=∫x²dA, Iyz=∫(y)(x)dA).
  const Iz = Ixx0 - A * cy * cy;
  const Iy = Iyy0 - A * cx * cx;
  const Iyz = Ixy0 - A * cx * cy;
  // Ejes principales.
  const avg = (Iz + Iy) / 2, dif = (Iz - Iy) / 2;
  const R = Math.hypot(dif, Iyz);
  const I1 = avg + R, I2 = avg - R;
  const theta = 0.5 * Math.atan2(-Iyz, dif);     // ángulo del eje principal mayor
  // Envolvente.
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const p of outline) { xmin = Math.min(xmin, p[0]); xmax = Math.max(xmax, p[0]); ymin = Math.min(ymin, p[1]); ymax = Math.max(ymax, p[1]); }
  const Sz = Iz / Math.max(ymax - cy, cy - ymin);
  const Sy = Iy / Math.max(xmax - cx, cx - xmin);
  // Eje neutro plástico (áreas iguales) por bisección y módulo plástico por recorte.
  const areaGE = (axis, c) => combine(outline, holes, 'A', axis, c);
  const bisect = (axis, lo, hi) => { for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; (areaGE(axis, m) > A / 2 ? lo = m : hi = m); } return (lo + hi) / 2; };
  const yp = bisect('y', ymin, ymax), xp = bisect('x', xmin, xmax);
  // Zz = 2·Qy_above − A·cy ;  Zy = 2·Qx_right − A·cx   (Qy_above=∫_{y≥yp} y dA)
  const QyAbove = combine(outline, holes, 'Qy', 'y', yp);
  const QxRight = combine(outline, holes, 'Qx', 'x', xp);
  const Zz = Math.abs(2 * QyAbove - A * cy);
  const Zy = Math.abs(2 * QxRight - A * cx);
  return {
    shape: 'polygon', A, cx, cy, Iz, Iy, Iyz, I1, I2, theta,
    Sz, Sy, Zz, Zy, perimeter: per,
    h: ymax - ymin, b: xmax - xmin, xmin, xmax, ymin, ymax,
  };
}
