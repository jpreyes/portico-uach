// ──────────────────────────────────────────────────────────────────────────────
// subspace.js — NÚCLEO COMPARTIDO de la iteración de subespacio (Bathe).
//
// Reúne las primitivas que el análisis MODAL (K·φ = ω²·M·φ) y el de PANDEO lineal
// ((K + λ·Kg)·φ = 0 → K·φ = λ·(−Kg)·φ) usan en común:
//   · `smallGenEig` — autovalores generalizados de un problema PEQUEÑO q×q
//     A·v = λ·B·v con B SPD (reducción de Cholesky + Jacobi). Es el corazón de la
//     fase Rayleigh-Ritz de la iteración de subespacio.
//   · `mvBand` / `dot` — producto matriz·vector dentro de la banda y producto punto.
//
// Pensado para correr en un Web Worker (módulo) y también para verificarse en Node.
// La iteración de subespacio MODAL vive en `modal_worker.js`; la de PANDEO en
// `buckling.js`. Ambas comparten estas piezas para que el método sea idéntico.
// ──────────────────────────────────────────────────────────────────────────────

// Producto A·x usando la extensión por filas (banda variable) → O(n·b).
export function mvBand(A, x, n, lo, hi) {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0; const off = i * n, a = lo[i], b = hi[i];
    for (let j = a; j <= b; j++) s += A[off + j] * x[j];
    y[i] = s;
  }
  return y;
}

export function dot(a, b, n) { let s = 0; for (let i = 0; i < n; i++) s += a[i] * b[i]; return s; }

// ── Autovalores generalizados pequeños (q×q):  A·v = λ·B·v, con B SPD ─────────
// Reducción de Cholesky B=L·Lᵀ → problema estándar C·y=λ·y con C=L⁻¹A L⁻ᵀ, luego
// Jacobi clásico.  Devuelve { vals (ascendente), vecs (columnas, B-ortonormales) }.
// (Antes era `_smallGenEig` privado en modal_worker.js; extraído sin cambios para
//  reutilizarlo en el motor de pandeo.)
export function smallGenEig(A, B, n) {
  // Cholesky de B (SPD tras regularización)
  const L = []; for (let i = 0; i < n; i++) L.push(new Float64Array(n));
  for (let j = 0; j < n; j++) {
    let s = B[j][j]; for (let k = 0; k < j; k++) s -= L[j][k] * L[j][k];
    L[j][j] = Math.sqrt(Math.max(s, 1e-300));
    for (let i = j + 1; i < n; i++) {
      let t = B[i][j]; for (let k = 0; k < j; k++) t -= L[i][k] * L[j][k];
      L[i][j] = t / L[j][j];
    }
  }
  const fwd = b => { const y = new Float64Array(n); for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; } return y; };
  const bwd = b => { const y = new Float64Array(n); for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * y[k]; y[i] = s / L[i][i]; } return y; };
  // C = L⁻¹ A L⁻ᵀ : por columnas, C = fwd(A·(L⁻ᵀ e_j)) — se simetriza al final
  const C = []; for (let i = 0; i < n; i++) C.push(new Float64Array(n));
  for (let j = 0; j < n; j++) {
    const ej = new Float64Array(n); ej[j] = 1;
    const w = bwd(ej);                                  // w = L⁻ᵀ e_j
    const Aw = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = 0; for (let k = 0; k < n; k++) s += A[i][k] * w[k]; Aw[i] = s; }
    const col = fwd(Aw);                                // C[:,j]
    for (let i = 0; i < n; i++) C[i][j] = col[i];
  }
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const m = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = m; }
  // Jacobi simétrico estándar → autovalores en diag(C), autovectores en V
  const V = []; for (let i = 0; i < n; i++) { V.push(new Float64Array(n)); V[i][i] = 1; }
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += C[p][q] * C[p][q];
    if (off < 1e-22) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(C[p][q]) < 1e-18) continue;
      const th = (C[q][q] - C[p][p]) / (2 * C[p][q]);
      const tt = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(tt * tt + 1), s = tt * c;
      for (let i = 0; i < n; i++) { const cip = C[i][p], ciq = C[i][q]; C[i][p] = c * cip - s * ciq; C[i][q] = s * cip + c * ciq; }
      for (let i = 0; i < n; i++) { const cpi = C[p][i], cqi = C[q][i]; C[p][i] = c * cpi - s * cqi; C[q][i] = s * cpi + c * cqi; }
      for (let i = 0; i < n; i++) { const vip = V[i][p], viq = V[i][q]; V[i][p] = c * vip - s * viq; V[i][q] = s * vip + c * viq; }
    }
  }
  // autovectores del problema generalizado: v = L⁻ᵀ y  (y = columnas de V)
  const pairs = [];
  for (let j = 0; j < n; j++) {
    const y = new Float64Array(n); for (let i = 0; i < n; i++) y[i] = V[i][j];
    pairs.push({ lam: C[j][j], v: bwd(y) });
  }
  pairs.sort((a, b) => a.lam - b.lam);
  return { vals: pairs.map(p => p.lam), vecs: pairs.map(p => p.v) };
}
