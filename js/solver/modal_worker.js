// ──────────────────────────────────────────────────────────────────────────────
// modal_worker.js — Classic Web Worker for Stodola inverse power iteration
//
// Runs the heavy modal solver off the main thread so the browser stays
// responsive during analysis.  Loaded as a CLASSIC worker (no {type:'module'})
// so that numeric.js can be imported via importScripts().
//
// Protocol:
//   Main → Worker: { Kff_flat: Float64Array, Mff_flat: Float64Array, nF, nModes }
//   Worker → Main: { modes: [{omega2, vec}] }  OR  { error: string }
//
// P4-15 Float64Array optimization:
//   · Kff_flat / Mff_flat are already Float64Arrays; used directly in _mv_f64
//     (flat row-major access → better cache locality, JIT-vectorizable).
//   · Only the Array-of-Arrays copies are built for numeric.LU / numeric.LUsolve.
//   · All iteration vectors (x, y, Mx, Kx, My) are Float64Arrays.
// ──────────────────────────────────────────────────────────────────────────────

importScripts('/lib/numeric.js');

// ── Flat-matrix × Float64Array vector  (P4-15: typed-array hot path) ─────────
function _mv_f64(A_flat, x, n) {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0;
    const off = i * n;
    for (let j = 0; j < n; j++) s += A_flat[off + j] * x[j];
    y[i] = s;
  }
  return y;
}

function _dot_f64(a, b, n) {
  let s = 0;
  for (let i = 0; i < n; i++) s += a[i] * b[i];
  return s;
}

// M-normalise in-place; returns same typed array
function _mNorm_f64(x, Mf, n) {
  const Mx   = _mv_f64(Mf, x, n);
  const norm = Math.sqrt(Math.max(_dot_f64(x, Mx, n), 0));
  if (norm < 1e-30) return x;
  for (let i = 0; i < n; i++) x[i] /= norm;
  return x;
}

// M-orthogonalise x (typed array) against found modes in-place
function _mOrtho_f64(x, found, Mf, n) {
  for (const { vec: phi } of found) {
    const Mphi = _mv_f64(Mf, phi, n);
    const c    = _dot_f64(x, Mphi, n);   // xᵀ M φ
    for (let i = 0; i < n; i++) x[i] -= c * phi[i];
  }
}

// ── Stodola inverse power iteration with M-orthogonal deflation ───────────────
function _stodola(Kf, KLU, Mf, Kff_aoa, nF, nModes, num) {
  const found = [];

  for (let modeNum = 0; modeNum < nModes; modeNum++) {
    let bestOmega2 = Infinity, bestVec = null;

    for (let attempt = 0; attempt < 6; attempt++) {
      const phase = (modeNum + 1 + attempt * 7) * 0.7 + attempt * 0.41;

      // Build start vector as Float64Array
      let x = new Float64Array(nF);
      for (let i = 0; i < nF; i++) {
        x[i] = Math.sin(phase * (i + 1)) + Math.cos((attempt + 1) * (i + 0.5) * 1.1) * 0.5 + 0.1;
      }
      _mOrtho_f64(x, found, Mf, nF);

      const Mx0 = _mv_f64(Mf, x, nF);
      const n0  = Math.sqrt(Math.max(_dot_f64(x, Mx0, nF), 0));
      if (n0 < 1e-10) continue;

      _mNorm_f64(x, Mf, nF);

      let omega2 = 0, converged = false;

      for (let iter = 0; iter < 150; iter++) {
        const Mx = _mv_f64(Mf, x, nF);

        // numeric.LUsolve expects regular Array; convert once per iter
        const Mx_arr = Array.from(Mx);
        const y_arr  = num.LUsolve(KLU, Mx_arr);

        // Back to Float64Array for fast orthogonalisation
        const y = new Float64Array(y_arr);
        _mOrtho_f64(y, found, Mf, nF);

        const My = _mv_f64(Mf, y, nF);
        const yn = Math.sqrt(Math.max(_dot_f64(y, My, nF), 0));
        if (yn < 1e-30) break;

        const xNew = y.slice();         // Float64Array.slice() — typed copy
        _mNorm_f64(xNew, Mf, nF);

        const Kx = _mv_f64(Kf, xNew, nF);
        const w2 = _dot_f64(xNew, Kx, nF);

        if (!isFinite(w2) || w2 < 0) break;

        const relChange = Math.abs(w2 - omega2) / Math.max(w2, 1e-10);
        omega2 = w2;
        x      = xNew;

        if (relChange < 1e-7 && iter >= 4)  { converged = true; break; }
        if (relChange < 1e-4 && iter >= 20) { converged = true; break; }
      }

      if (converged && isFinite(omega2) && omega2 >= 0 && omega2 < 1e12) {
        if (omega2 < bestOmega2) {
          bestOmega2 = omega2;
          bestVec    = Array.from(x);   // store as plain Array for postMessage
        }
      }
    }

    if (!bestVec) break;
    found.push({ omega2: bestOmega2, vec: new Float64Array(bestVec) });
  }

  // Convert found vecs back to plain Arrays for structured-clone transfer
  return found.map(m => ({ omega2: m.omega2, vec: Array.from(m.vec) }));
}

// ── Message handler ────────────────────────────────────────────────────────────
self.onmessage = function (e) {
  const { Kff_flat, Mff_flat, nF, nModes } = e.data;
  try {
    // Build Array-of-Arrays for numeric.LU (required by numeric.js API)
    const Kff_aoa = [];
    const Mff_aoa = [];
    for (let i = 0; i < nF; i++) {
      Kff_aoa.push(Array.from(Kff_flat.subarray(i * nF, (i + 1) * nF)));
      Mff_aoa.push(Array.from(Mff_flat.subarray(i * nF, (i + 1) * nF)));
    }

    // Regularise M diagonal
    let maxMd = 0;
    for (let i = 0; i < nF; i++) maxMd = Math.max(maxMd, Math.abs(Mff_aoa[i][i]));
    if (maxMd < 1e-30) {
      self.postMessage({ error: 'Matriz de masas nula. Asigne densidad ρ a los materiales o masa a los diafragmas.' });
      return;
    }
    const eps = maxMd * 1e-8;
    for (let i = 0; i < nF; i++) {
      if (Math.abs(Mff_aoa[i][i]) < eps) {
        Mff_aoa[i][i] = eps;
        Mff_flat[i * nF + i] = eps;   // keep flat copy in sync
      }
    }

    // Factor K once (Array-of-Arrays required by numeric.js)
    let KLU;
    try { KLU = numeric.LU(Kff_aoa); } catch (e2) {
      self.postMessage({ error: 'Factorización K falló: ' + (e2.message || e2) + '. Verifique estabilidad del modelo.' });
      return;
    }

    // Stodola iteration — hot path uses Kff_flat / Mff_flat directly
    const modes = _stodola(Kff_flat, KLU, Mff_flat, Kff_aoa, nF, nModes, numeric);
    if (!modes.length) {
      self.postMessage({ error: 'Sin modos estructurales. Verifique masa (ρ en material o diafragmas) y apoyos.' });
      return;
    }

    self.postMessage({ modes });
  } catch (err) {
    self.postMessage({ error: (err && err.message) ? err.message : String(err) });
  }
};
