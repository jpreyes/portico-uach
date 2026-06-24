// ──────────────────────────────────────────────────────────────────────────────
// io/ifc/ifcWarnings.js — recolector de ADVERTENCIAS de importación IFC (#76, G19)
//
// Centraliza los avisos no bloqueantes (material genérico, sección aproximada, curva
// segmentada, etc.).  La filosofía de G19 es ROBUSTA: nunca se aborta la importación por
// un dato que falta o no se entiende; en su lugar se deja una advertencia y se sigue con
// un valor genérico razonable.  Las advertencias se muestran por elemento (en la tabla del
// diálogo) y agregadas (resumen global).  AUTÓNOMO (Node + navegador).
// ──────────────────────────────────────────────────────────────────────────────

/** Pequeño acumulador de advertencias con conteo de repetidas (para el resumen global). */
export class Warnings {
  constructor() { this.list = []; this._counts = new Map(); }

  /** Agrega una advertencia (texto libre). Devuelve `this` para encadenar. */
  add(msg) {
    if (!msg) return this;
    this.list.push(msg);
    this._counts.set(msg, (this._counts.get(msg) || 0) + 1);
    return this;
  }

  get length() { return this.list.length; }
  get empty()  { return this.list.length === 0; }

  /** Texto plano de todas las advertencias (para `_alert`). */
  toText() { return this.list.join('\n'); }

  /** Resumen agrupado «(×N) mensaje», ordenado por frecuencia descendente. */
  summary() {
    return [...this._counts.entries()]
      .sort((a, b) => b[1] - a[1])
      .map(([m, c]) => (c > 1 ? `(×${c}) ${m}` : m));
  }
}
