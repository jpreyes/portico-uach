// ──────────────────────────────────────────────────────────────────────────────
// static_worker.js — resuelve K·u = F para TODOS los casos estáticos fuera del
// hilo principal (la UI no se congela). Module worker (importa linsolve.js ESM).
//
//   Main → Worker: { Kflat: Float64Array(nDOF²), nDOF, freeDOF: Int32Array, Flist: [Float64Array(nDOF)] }
//   Worker → Main: { progress, done, total }   (avance)
//                  { ok:true, uList, reactionsList, bandwidth }   (éxito)
//                  { ok:false, error? }   (no SPD / inestable → el main usa respaldo)
//
// Estrategia: extrae K_ff, factoriza UNA vez (Cholesky en banda con RCM) y
// resuelve cada lado derecho. Reacciones = K·u − F (en el worker, no bloquea).
// ──────────────────────────────────────────────────────────────────────────────
import { factorSolveMany } from './linsolve.js?v=73';

self.onmessage = (e) => {
  const { Kflat, nDOF, freeDOF, Flist } = e.data;
  try {
    const nF = freeDOF.length;
    if (nF === 0) { self.postMessage({ ok: false, error: 'sin GDL libres' }); return; }

    // Extraer K_ff (libre–libre)
    const Kff = new Float64Array(nF * nF);
    for (let i = 0; i < nF; i++) {
      const rowK = freeDOF[i] * nDOF, rowF = i * nF;
      for (let j = 0; j < nF; j++) Kff[rowF + j] = Kflat[rowK + freeDOF[j]];
    }
    // Reducir cada F a los GDL libres
    const FfList = Flist.map(F => { const ff = new Float64Array(nF); for (let i = 0; i < nF; i++) ff[i] = F[freeDOF[i]]; return ff; });

    self.postMessage({ progress: 'factorizando', done: 0, total: Flist.length });
    const res = factorSolveMany(Kff, nF, FfList);
    if (!res.ok) { self.postMessage({ ok: false }); return; }   // no SPD → respaldo en el main

    const uList = [], reactionsList = [];
    for (let c = 0; c < Flist.length; c++) {
      const uf = res.uList[c];
      const u = new Float64Array(nDOF);
      for (let i = 0; i < nF; i++) u[freeDOF[i]] = uf[i];
      // reacciones = K·u − F
      const F = Flist[c];
      const reac = new Float64Array(nDOF);
      for (let i = 0; i < nDOF; i++) {
        let s = 0; const off = i * nDOF;
        for (let j = 0; j < nDOF; j++) s += Kflat[off + j] * u[j];
        reac[i] = s - F[i];
      }
      uList.push(u); reactionsList.push(reac);
      self.postMessage({ progress: 'resolviendo', done: c + 1, total: Flist.length });
    }
    self.postMessage({ ok: true, uList, reactionsList, bandwidth: res.bandwidth });
  } catch (err) {
    self.postMessage({ ok: false, error: String((err && err.message) || err) });
  }
};
