// runners.mjs — corre los solvers de Pórtico HEADLESS (Node), reusando el código
// real de la app. numeric.js se carga como global (shim de window) una sola vez.
import { ModalSolver } from '../../js/solver/modal_solver.js';

let _num = false;
export async function ensureNumeric() {
  if (_num) return;
  globalThis.window = globalThis;
  await import('../../lib/numeric.js');           // define `numeric` global
  globalThis.window.numeric = globalThis.numeric;
  _num = true;
}

// Análisis modal — devuelve el ModalResults real (period[], freq[], getModeShape…).
export async function runModal(model, nModes = 6) {
  await ensureNumeric();
  return new ModalSolver().solve(model, nModes);
}

// Despacho por tipo de análisis (se irá ampliando: static, buckling, …).
export async function runAnalysis(model, spec) {
  switch (spec.analysis) {
    case 'modal': return { type: 'modal', res: await runModal(model, spec.nModes || 6) };
    default: throw new Error('Análisis no soportado en el harness: ' + spec.analysis);
  }
}
