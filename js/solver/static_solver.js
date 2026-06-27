// ──────────────────────────────────────────────────────────────────────────────
// StaticSolver — direct stiffness method for linear static analysis
// Solver:  K_ff · u_f = F_f  (Gaussian elimination via numeric.js)
// ──────────────────────────────────────────────────────────────────────────────
import { buildNodeIndex, assembleK, assembleF, getNodeDOFs } from './assembler.js?v=211';
import { Results } from './postprocess.js?v=211';

// Curva de suelo no lineal (#4): fuerza restauradora y tangente por interpolación lineal
// de la tabla [[d,F],…] (monótona en d). Fuera del rango → fuerza última (saturación).
function soilForce(curve, d) {
  const n = curve.length;
  if (n === 0) return [0, 0];
  if (d <= curve[0][0]) return [curve[0][1], 0];
  if (d >= curve[n - 1][0]) return [curve[n - 1][1], 0];
  for (let i = 1; i < n; i++) {
    if (d <= curve[i][0]) {
      const d0 = curve[i - 1][0], f0 = curve[i - 1][1], d1 = curve[i][0], f1 = curve[i][1];
      const kt = (f1 - f0) / (d1 - d0);
      return [f0 + kt * (d - d0), kt];
    }
  }
  return [curve[n - 1][1], 0];
}

export class StaticSolver {
  /**
   * @param {Model}       model
   * @param {number|null} lcId        Load case ID (null → pure self-weight)
   * @param {boolean}     selfWeight
   */
  solve(model, lcId = null, selfWeight = false) {
    // ── Build index & global matrices ─────────────────────────────────────
    const nodeIndex = buildNodeIndex(model);
    const { K, nDOF } = assembleK(model, nodeIndex);
    const F = assembleF(model, nodeIndex, lcId, selfWeight);

    // ── Classify DOFs ──────────────────────────────────────────────────────
    const freeDOF  = [];
    const fixedDOF = [];

    // Modelo 2D (pórtico plano X–Z): los GDL fuera del plano (uy, rx, rz) se
    // restringen automáticamente en todos los nodos — el usuario solo maneja
    // los GDL del plano: ux, uz y el giro ry.
    const is2D = model.mode === '2D';

    // Desplazamiento prescrito (#54): valor impuesto por GDL restringido.
    // up[gi] = desplazamiento conocido del GDL soporte (0 = apoyo normal).
    const up = new Float64Array(nDOF);
    let hasPresc = false;
    const dofNames = ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'];

    for (const node of model.nodes.values()) {
      const d    = getNodeDOFs(nodeIndex, node.id);
      const r    = node.restraints;
      const pd   = node.prescDisp;
      const rArr = [
        r.ux,
        is2D ? 1 : r.uy,
        r.uz,
        is2D ? 1 : r.rx,
        r.ry,
        is2D ? 1 : r.rz,
      ];
      d.forEach((gi, li) => {
        const pv = pd ? (+pd[dofNames[li]] || 0) : 0;   // valor prescrito de este GDL
        if (rArr[li] || pv !== 0) {
          fixedDOF.push(gi);
          if (pv !== 0) { up[gi] = pv; hasPresc = true; }
        } else {
          freeDOF.push(gi);
        }
      });
    }

    if (freeDOF.length === 0) {
      throw new Error('El modelo no tiene grados de libertad libres (¿todos los nodos están empotrados?)');
    }

    // ── Extract K_ff and F_f ──────────────────────────────────────────────
    const nF  = freeDOF.length;
    const Kff = Array.from({ length: nF }, (_, i) =>
      Array.from({ length: nF }, (_, j) => K[freeDOF[i] * nDOF + freeDOF[j]])
    );
    // F_f efectivo = F_f − K_fp·u_p   (traslada el desplazamiento prescrito al RHS)
    const Ff = freeDOF.map((di, i) => {
      let f = F[di];
      if (hasPresc) for (const dj of fixedDOF) { if (up[dj]) f -= K[di * nDOF + dj] * up[dj]; }
      return f;
    });

    // ── Solve ─────────────────────────────────────────────────────────────
    const num = window.numeric;
    if (!num) throw new Error('numeric.js no está disponible');

    // Resortes UNILATERALES (#3): se resuelven por CONJUNTO ACTIVO. assembleK los metió
    // como bilaterales; aquí se desactivan (restando su k de la diagonal de Kff) los que
    // estén «descomprimidos» y se re-resuelve hasta que el conjunto activo se estabiliza
    // (uplift/gap). 'c' = solo-compresión (activo si u≤0), 't' = solo-tracción (u≥0).
    const freeIdx = new Map(freeDOF.map((g, i) => [g, i]));
    const uni = [];
    for (const node of model.nodes.values()) {
      const su = node.springUni, sp = node.springs;
      if (!su || !sp) continue;
      const ks = [sp.kux, sp.kuy, sp.kuz, sp.krx, sp.kry, sp.krz];
      const md = [su.ux, su.uy, su.uz, su.rx, su.ry, su.rz];
      const d  = getNodeDOFs(nodeIndex, node.id);
      for (let i = 0; i < 6; i++)
        if ((md[i] === 'c' || md[i] === 't') && ks[i] > 0 && freeIdx.has(d[i]))
          uni.push({ fi: freeIdx.get(d[i]), k: ks[i], mode: md[i] });
    }

    const solveLinear = () => {
      let x;
      try { x = num.solve(Kff, Ff); } catch (e) { throw new Error(`Solver falló: ${e.message}. Verifique que el modelo es estable.`); }
      if (!x || x.some(v => !Number.isFinite(v)))
        throw new Error(
          'Estructura INESTABLE (matriz singular): existe un mecanismo. ' +
          'Revise apoyos y liberaciones — p.ej. liberar el mismo giro en ambos ' +
          'extremos de elementos contiguos permite rotación libre' +
          (uni.length ? ', o un apoyo unilateral que se despegó dejó la estructura sin sustento.' : '.')
        );
      return x;
    };

    let uf;
    const inactive = new Set();   // GDL globales con resorte unilateral desactivado
    if (!uni.length) {
      uf = solveLinear();
    } else {
      const active = uni.map(() => true);   // arranca con todos activos (Kff tal cual)
      for (let pass = 0; pass <= uni.length + 2; pass++) {
        for (let s = 0; s < uni.length; s++) if (!active[s]) Kff[uni[s].fi][uni[s].fi] -= uni[s].k;
        try { uf = solveLinear(); }
        finally { for (let s = 0; s < uni.length; s++) if (!active[s]) Kff[uni[s].fi][uni[s].fi] += uni[s].k; }
        let changed = false;
        for (let s = 0; s < uni.length; s++) {
          const ud = uf[uni[s].fi];
          const want = uni[s].mode === 'c' ? (ud <= 1e-12) : (ud >= -1e-12);   // ¿debe estar activo?
          if (want !== active[s]) { active[s] = want; changed = true; }
        }
        if (!changed) break;
      }
      uni.forEach((s, k) => { if (!active[k]) inactive.add(freeDOF[s.fi]); });
    }

    // ── Resortes de suelo NO LINEALES (#4): Newton sobre la curva p-y/t-z/q-z ──
    // No están en assembleK (son no lineales): aquí se itera  Kt·Δu = Fext − K·u − R(u)
    // con R y la tangente kt del segmento actual de la curva. Saturan a la fuerza última.
    const soil = [];
    for (const node of model.nodes.values()) {
      const ss = node.soilSpring; if (!ss) continue;
      const d = getNodeDOFs(nodeIndex, node.id);
      ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'].forEach((nm, i) => {
        if (Array.isArray(ss[nm]) && ss[nm].length >= 2 && freeIdx.has(d[i])) soil.push({ fi: freeIdx.get(d[i]), gi: d[i], curve: ss[nm] });
      });
    }
    if (soil.length) {
      const fref = Math.sqrt(Ff.reduce((a, b) => a + b * b, 0)) || 1;
      for (let it = 0; it < 60; it++) {
        const Kt = Kff.map(row => row.slice());          // tangente (no muta Kff)
        const r = new Float64Array(nF);
        for (let i = 0; i < nF; i++) { let s = 0; const Ki = Kff[i]; for (let j = 0; j < nF; j++) s += Ki[j] * uf[j]; r[i] = Ff[i] - s; }
        for (const sp of soil) { const [F, kt] = soilForce(sp.curve, uf[sp.fi]); r[sp.fi] -= F; Kt[sp.fi][sp.fi] += kt; }
        let rn = 0; for (let i = 0; i < nF; i++) rn += r[i] * r[i];
        if (Math.sqrt(rn) / fref < 1e-9) break;
        let du; try { du = num.solve(Kt, Array.from(r)); } catch { break; }
        if (!du || du.some(v => !Number.isFinite(v))) break;
        for (let i = 0; i < nF; i++) uf[i] += du[i];
      }
    }

    // ── Assemble full displacement vector ──────────────────────────────────
    const u = new Float64Array(nDOF);
    freeDOF.forEach((d, i) => { u[d] = uf[i]; });
    if (hasPresc) for (const d of fixedDOF) if (up[d]) u[d] = up[d];   // GDL prescritos (#54)

    // ── Compute reactions ──────────────────────────────────────────────────
    const reactions = new Float64Array(nDOF);
    for (let i = 0; i < nDOF; i++) {
      let Ku_i = 0;
      for (let j = 0; j < nDOF; j++) Ku_i += K[i * nDOF + j] * u[j];
      reactions[i] = Ku_i - F[i];
    }

    // Reacciones de apoyos elásticos: SOLO en GDL libres (donde el resorte es el
    // apoyo). Ahí el balance Ku−F vale 0 y la reacción real es −k·u. Si el GDL
    // también está rígidamente restringido, el resorte no actúa (u=0): se conserva
    // la reacción rígida Ku−F y NO se sobrescribe (de lo contrario daría 0).
    const freeSet = new Set(freeDOF);
    for (const node of model.nodes.values()) {
      const sp = node.springs;
      if (!sp) continue;
      const ks = [sp.kux, sp.kuy, sp.kuz, sp.krx, sp.kry, sp.krz];
      if (!ks.some(k => k > 0)) continue;
      const d = getNodeDOFs(nodeIndex, node.id);
      for (let i = 0; i < 6; i++) {
        if (ks[i] > 0 && freeSet.has(d[i])) reactions[d[i]] = inactive.has(d[i]) ? 0 : -ks[i] * u[d[i]];   // resorte unilateral despegado → 0 (#3)
      }
    }
    // Reacción de los resortes de suelo no lineales (#4): −R(u) en su GDL libre.
    for (const sp of soil) if (freeSet.has(sp.gi)) reactions[sp.gi] = -soilForce(sp.curve, u[sp.gi])[0];

    return new Results(model, nodeIndex, u, reactions, F, lcId, selfWeight);
  }
}


// ── Default load case manager helper ─────────────────────────────────────────
export function ensureDefaultLC(model) {
  // CM (carga muerta) incluye el peso propio por defecto
  if (model.loadCases.size === 0) model.addLoadCase('CM', true);
  return model.loadCases.keys().next().value;
}

