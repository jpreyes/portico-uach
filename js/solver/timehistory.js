// ──────────────────────────────────────────────────────────────────────────────
// timehistory.js — Análisis dinámico en el TIEMPO (time-history) LINEAL por
// SUPERPOSICIÓN MODAL (G12 · #48a).
//
// Excitación: acelerograma UNIFORME en la base a_g(t) (mismo movimiento en todos
// los apoyos). Para cada modo i (ωᵢ, factor de participación Γᵢ = φᵢᵀM·r / φᵢᵀMφᵢ)
// la coordenada modal qᵢ(t) resuelve la SDOF
//
//     q̈ᵢ + 2ζωᵢ q̇ᵢ + ωᵢ² qᵢ = −Γᵢ · a_g(t)
//
// integrada por la DUHAMEL evaluada con la recurrencia EXACTA de Nigam–Jennings
// (Chopra, «Dynamics of Structures», Tabla 5.4.1): es exacta cuando a_g(t) es
// lineal por tramos (la interpolación natural de un registro muestreado a Δt
// constante) e incondicionalmente estable. La respuesta física se reconstruye
// por superposición:  u(t) = Σ φᵢ qᵢ(t).
//
// Núcleo AUTÓNOMO (sin dependencias) para verificarlo en Node contra soluciones
// analíticas SDOF (escalón → DLF=2; armónico → función de transferencia;
// vibración libre → decremento logarítmico). Ver test_timehistory.mjs.
// ──────────────────────────────────────────────────────────────────────────────

// ── Coeficientes de la recurrencia exacta de Nigam–Jennings ───────────────────
// Para  ü + 2ζω u̇ + ω² u = p(t)  con p lineal por tramos y paso Δt constante:
//   u_{k+1} = A·u_k + B·u̇_k + C·p_k + D·p_{k+1}
//   u̇_{k+1} = A'·u_k + B'·u̇_k + C'·p_k + D'·p_{k+1}
// Válido para 0 ≤ ζ < 1. (Se asume masa modal unitaria: p está en unidades de
// aceleración, como −Γ·a_g.)
export function njCoeffs(omega, zeta, dt) {
  const w1 = Math.sqrt(Math.max(1 - zeta * zeta, 1e-300));   // √(1−ζ²)
  const wd = omega * w1;                                       // ω amortiguada
  const e = Math.exp(-zeta * omega * dt);
  const s = Math.sin(wd * dt), c = Math.cos(wd * dt);
  const zr = zeta / w1;                                        // ζ/√(1−ζ²)
  const w2 = omega * omega;
  const zod = 2 * zeta / (omega * dt);                         // 2ζ/(ωΔt)

  const A = e * (zr * s + c);
  const B = e * (s / wd);
  const C = (1 / w2) * (zod + e * (((1 - 2 * zeta * zeta) / (wd * dt) - zr) * s - (1 + zod) * c));
  const D = (1 / w2) * (1 - zod + e * (((2 * zeta * zeta - 1) / (wd * dt)) * s + zod * c));

  const Ap = -e * (omega / w1) * s;                            // −e·(ω/√(1−ζ²))·s
  const Bp = e * (c - zr * s);
  const Cp = (1 / w2) * (-1 / dt + e * ((omega / w1 + zr / dt) * s + (1 / dt) * c));
  const Dp = (1 / (w2 * dt)) * (1 - e * (zr * s + c));

  return { A, B, C, D, Ap, Bp, Cp, Dp };
}

// ── Respuesta SDOF a una carga p[k] (lineal por tramos) ───────────────────────
// Resuelve ü + 2ζω u̇ + ω² u = p(t) y devuelve u (desplazamiento) y v (velocidad)
// en los mismos instantes que p. u0/v0 = condiciones iniciales (def. reposo).
export function sdofResponse(omega, zeta, dt, p, u0 = 0, v0 = 0) {
  const n = p.length;
  const u = new Float64Array(n), v = new Float64Array(n);
  const { A, B, C, D, Ap, Bp, Cp, Dp } = njCoeffs(omega, zeta, dt);
  u[0] = u0; v[0] = v0;
  for (let k = 0; k < n - 1; k++) {
    u[k + 1] = A * u[k] + B * v[k] + C * p[k] + D * p[k + 1];
    v[k + 1] = Ap * u[k] + Bp * v[k] + Cp * p[k] + Dp * p[k + 1];
  }
  return { u, v };
}

// ── Espectro de respuesta a partir de un acelerograma ─────────────────────────
// Para cada período T de `periods`, integra la SDOF de frecuencia ω=2π/T con
// amortiguamiento ζ bajo a_g(t) y devuelve Sd, Sv (pseudo), Sa (pseudo) máximas.
// (Auxiliar de verificación y de uso futuro.)
export function responseSpectrum(ag, dt, periods, zeta = 0.05) {
  return periods.map(T => {
    const w = 2 * Math.PI / T;
    const p = new Float64Array(ag.length);
    for (let k = 0; k < ag.length; k++) p[k] = -ag[k];   // ü + … = −a_g
    const { u } = sdofResponse(w, zeta, dt, p);
    let Sd = 0; for (let k = 0; k < u.length; k++) Sd = Math.max(Sd, Math.abs(u[k]));
    return { T, w, Sd, Sv: w * Sd, Sa: w * w * Sd };
  });
}

// ── Time-history modal completo ───────────────────────────────────────────────
/**
 * @param {object} o
 *   modes  [{ omega, gamma, phi }]   ω (rad/s), Γ participación en la dirección de
 *                                    excitación, phi = forma modal (Float64Array nDOF)
 *   ag     Float64Array              acelerograma basal (m/s²) muestreado a Δt
 *   dt     number                    paso de tiempo del registro (s)
 *   zeta   number | number[]         amortiguamiento (escalar o por modo, def. 0.05)
 * @returns {object}
 *   t       Float64Array(nSteps)     instantes
 *   q       Float64Array[nModes]     coordenada modal qᵢ(t)
 *   nSteps, nModes
 *   nodalDOF(dof)        → Float64Array(nSteps)   historia de un GDL global
 *   uAt(step)           → Float64Array(nDOF)      desplazamientos en un instante
 *   peakModal           Float64Array(nModes)      |qᵢ| máximo por modo
 */
export function modalTimeHistory(o) {
  const modes = o.modes, ag = o.ag, dt = o.dt;
  const nSteps = ag.length, nModes = modes.length;
  const zArr = Array.isArray(o.zeta) ? o.zeta : modes.map(() => (o.zeta ?? 0.05));

  const t = new Float64Array(nSteps);
  for (let k = 0; k < nSteps; k++) t[k] = k * dt;

  const q = [];
  const peakModal = new Float64Array(nModes);
  for (let i = 0; i < nModes; i++) {
    const G = modes[i].gamma;
    const p = new Float64Array(nSteps);
    for (let k = 0; k < nSteps; k++) p[k] = -G * ag[k];     // −Γ·a_g
    const { u } = sdofResponse(modes[i].omega, zArr[i], dt, p);
    q.push(u);
    let pk = 0; for (let k = 0; k < nSteps; k++) pk = Math.max(pk, Math.abs(u[k]));
    peakModal[i] = pk;
  }

  const nDOF = modes.length ? modes[0].phi.length : 0;
  return {
    t, q, nSteps, nModes, peakModal,
    // Historia de un GDL global por superposición:  u_dof(t) = Σ φᵢ[dof]·qᵢ(t)
    nodalDOF(dof) {
      const h = new Float64Array(nSteps);
      for (let i = 0; i < nModes; i++) {
        const c = modes[i].phi[dof], qi = q[i];
        if (c === 0) continue;
        for (let k = 0; k < nSteps; k++) h[k] += c * qi[k];
      }
      return h;
    },
    // Vector de desplazamientos completo en el paso `step`.
    uAt(step) {
      const u = new Float64Array(nDOF);
      for (let i = 0; i < nModes; i++) {
        const qi = q[i][step], phi = modes[i].phi;
        if (qi === 0) continue;
        for (let d = 0; d < nDOF; d++) u[d] += phi[d] * qi;
      }
      return u;
    }
  };
}
