// ──────────────────────────────────────────────────────────────────────────────
// cargas.js — magnitudes de cargas normativas a partir de la ficha + reglas.json
//   · cargaNieveNCh431  → pf, ps (kN/m²)
//   · cargaVientoNCh432 → q (N/m²) y presiones por zona
//   · espectroNCh433    → curva elástica Sa(T) [g] + saFactor (g/R*) para PÓRTICO
// Módulo ES puro (Node + navegador). Sin DOM. Determinista y auditable.
// ──────────────────────────────────────────────────────────────────────────────

const G_GRAV = 9.80665; // m/s² (también kN por tonelada-fuerza)

// ── helpers de tablas escalón / bandas ────────────────────────────────────────

/** Índice de banda "lo-hi" (ej. "38-42") que contiene x; -1 si ninguna. */
function indiceBandaRango(bandas, x) {
  for (let i = 0; i < bandas.length; i++) {
    const [lo, hi] = String(bandas[i]).split('-').map(parseFloat);
    if (x >= lo && (i === bandas.length - 1 ? x <= hi : x < hi)) return i;
  }
  return -1;
}

/** Clave de banda de altitud ("0-300"…">4000") que contiene alt. */
function bandaAltitud(porAltitud, alt) {
  for (const k of Object.keys(porAltitud)) {
    if (k.startsWith('>')) { if (alt >= parseFloat(k.slice(1))) return k; }
    else { const [lo, hi] = k.split('-').map(parseFloat); if (alt >= lo && alt < hi) return k; }
  }
  return null;
}

/** Lookup en tabla escalón [[limSupExcl|null, val]] (último null = ∞). */
function escalon(tabla, x) {
  for (const [ub, val] of tabla) if (ub === null || x < ub) return val;
  return tabla[tabla.length - 1][1];
}

// ── NIEVE (NCh431.Of2010) ─────────────────────────────────────────────────────

/**
 * @returns {object} { pg, Ce, Ct, I, Cs, pf, ps, _notas }  (kN/m²); pf/ps null si pg null.
 */
export function cargaNieveNCh431(ficha, reglas) {
  const n = reglas.cargas.nieve_NCh431;
  const ub = ficha.ubicacion || {};
  const lat = ub.latitud_sur_deg, alt = ub.altitud_msnm ?? 0;
  const notas = [];

  const t1 = n.pg_carga_basica_terreno_Tabla1;
  const li = indiceBandaRango(t1.latitudes_sur, lat);
  const ab = bandaAltitud(t1.por_altitud_msnm, alt);
  if (li < 0 || ab == null) return { pg: null, pf: null, ps: null, _notas: ['latitud/altitud fuera de la Tabla 1'] };
  const pg = t1.por_altitud_msnm[ab][li];
  if (pg == null) notas.push(`Tabla 1 sin dato para lat ${lat}°, altitud ${alt} m (banda ${t1.latitudes_sur[li]} / ${ab})`);

  // Factores: defaults conservadores; la ficha puede afinarlos.
  const expo = ub.exposicion || 'C';
  const proteccion = ficha.cargas?.proteccion_nieve || 'Parcialmente expuesto';
  const Ce = n.Ce_factor_exposicion_Tabla4[expo]?.[proteccion] ?? 1.0;
  const Ct = n.Ct_factor_termico_Tabla2['Todas las estructuras (salvo las siguientes)'] ?? 1.0;
  const cat = ficha.sismo?.categoria || ficha.categoria_ocupacion || 'II';
  const I = n.I_factor_importancia_Tabla3[cat] ?? 1.0;

  // Pendiente del techo → Cs (forma ASCE; VALIDAR Figura 1).
  const pend = ficha.geometria?.pendiente_techo_deg ?? 0;
  const csQuiebre = Ct >= 1.2 ? 45 : Ct >= 1.1 ? 37.5 : 30;
  const Cs = pend <= csQuiebre ? 1.0 : pend >= 70 ? 0 : (70 - pend) / (70 - csQuiebre);

  const pf = pg == null ? null : +(0.7 * Ce * Ct * I * pg).toFixed(4);
  const ps = pf == null ? null : +(Cs * pf).toFixed(4);
  return { pg, Ce, Ct, I, Cs: +Cs.toFixed(3), pf, ps, _notas: notas, _formula: 'pf=0.7·Ce·Ct·I·pg ; ps=Cs·pf' };
}

// ── VIENTO (NCh432.Of2010) ────────────────────────────────────────────────────

/** Velocidad básica V (m/s): por estación (match de ciudad) o por banda de latitud. */
function velocidadViento(w, ub) {
  const ciudad = (ub.ciudad || '').trim().toLowerCase();
  if (ciudad) {
    for (const [k, v] of Object.entries(w.V_basica_m_s.por_estacion))
      if (k.toLowerCase() === ciudad) return { V: v, fuente: `estación ${k}` };
  }
  // bandas de latitud: límites superiores alineados al array (robusto a edición de V)
  const lat = ub.latitud_sur_deg;
  const ubExcl = [27, 35, 42, 50, Infinity];
  const arr = w.V_basica_m_s.por_latitud;
  for (let i = 0; i < arr.length && i < ubExcl.length; i++)
    if (lat < ubExcl[i]) return { V: arr[i].V, fuente: `latitud ${arr[i].rango}` };
  return { V: arr[arr.length - 1].V, fuente: `latitud ${arr[arr.length - 1].rango}` };
}

/**
 * @param {number} h_techo  altura de techo (m)
 * @returns {object} { V, Kz, Kzt, Kd, I, q_Nm2, GCpi, presiones:{zona:p_Nm2}, _notas }
 */
export function cargaVientoNCh432(ficha, reglas, h_techo) {
  const w = reglas.cargas.viento_NCh432;
  const ub = ficha.ubicacion || {};
  const expo = ub.exposicion || 'C';
  const e = w.exposicion[expo];
  const { V, fuente } = velocidadViento(w, ub);

  const z = Math.max(h_techo, 4.6);
  const Kz = 2.01 * Math.pow(z / e.Zg_m, 2 / e.alfa);
  const Kzt = 1.0; // sin datos topográficos en la ficha → terreno plano
  const Kd = w.Kd_direccionalidad['Edificio SPRFV'] ?? 0.85;
  const cat = ficha.sismo?.categoria || ficha.categoria_ocupacion || 'II';
  const I = w.factor_importancia_I[cat] ?? 1.0;
  const q = 0.613 * Kz * Kzt * Kd * V * V * I; // N/m²

  const cierre = ficha.cargas?.cierre_viento || 'Cerrado';
  const GCpi = w.GCpi_presion_interna[cierre] ?? 0.18;

  // tramo de pendiente para GCpf
  const pend = ficha.geometria?.pendiente_techo_deg ?? 0;
  const tr = pend < 5 ? 0 : pend < 20 ? 1 : pend < 45 ? 2 : 3;
  const presiones = {};
  for (const [zona, vals] of Object.entries(w.GCpf_coef_presion_externa.zonas)) {
    const g = vals[tr];
    const a = q * (g - GCpi), b = q * (g + GCpi);
    presiones[zona] = +(Math.abs(a) > Math.abs(b) ? a : b).toFixed(2);
  }
  return {
    V, fuente_V: fuente, Kz: +Kz.toFixed(4), Kzt, Kd, I,
    q_Nm2: +q.toFixed(2), GCpi, cierre, tramo_pendiente: tr, presiones,
    _formula: 'q=0.613·Kz·Kzt·Kd·V²·I [N/m²] ; p=q·(GCpf−GCpi)',
    _notas: ['Kzt=1 (terreno plano: la ficha no trae topografía)'],
  };
}

// ── SISMO (NCh433/DS61) ───────────────────────────────────────────────────────

/** R* = 1 + T* / (0.10·To + T* / Ro), con T* del análisis modal. */
export function Rstar(Tstar, To, Ro = 11.0) {
  return 1 + Tstar / (0.10 * To + Tstar / Ro);
}

/**
 * Espectro ELÁSTICO de diseño NCh433: Sa(T) = S·Ao·I·α(T) en unidades g.
 * La reducción por R* (que requiere T* del modal) se aplica como saFactor.
 * @returns {object} { curva:[{T,Sa}], texto, params, saFactor_nota, Rstar_formula }
 */
export function espectroNCh433(ficha, reglas, { Tmax = 3.0, dT = 0.02 } = {}) {
  const s = reglas.cargas.sismica_NCh433;
  const p = ficha.sismo || {};
  const def = s.parametros_defecto || {};
  const suelo = s.tabla_suelos[p.suelo];
  if (!suelo) throw new Error(`Suelo NCh433 no válido: "${p.suelo}"`);
  const Ao = s.tabla_zona_Ao_g[String(p.zona)];
  if (Ao == null) throw new Error(`Zona sísmica no válida: ${p.zona}`);
  const I = s.tabla_categoria_I[p.categoria || 'II'] ?? 1.0;
  const { S, To, Tp, n: nn, p: pp } = suelo;

  const alpha = (T) => (1 + 4.5 * Math.pow(T / To, pp)) / (1 + Math.pow(T / To, 3));
  const curva = [];
  for (let T = 0; T <= Tmax + 1e-9; T += dT) {
    const Tr = +T.toFixed(4);
    curva.push({ T: Tr, Sa: +(S * Ao * I * alpha(Tr)).toFixed(6) });
  }
  const texto = curva.map((q) => `${q.T}\t${q.Sa}`).join('\n');
  const Ro = def.Ro ?? 11.0, R = p.R ?? def.R ?? 7.0;
  return {
    curva, texto,
    params: { S, Ao, I, To, Tp, n: nn, p: pp, Ro, R },
    Rstar_formula: 'R* = 1 + T*/(0.10·To + T*/Ro), con T* del modal',
    saFactor_nota: `Pegar la curva en F7; saFactor = ${G_GRAV.toFixed(5)}/R* (convierte g→m/s² y aplica R*). R* tras el modal con Rstar(T*, To=${To}, Ro=${Ro}).`,
  };
}
