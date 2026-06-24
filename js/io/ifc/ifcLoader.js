// ──────────────────────────────────────────────────────────────────────────────
// io/ifc/ifcLoader.js — lector de archivos IFC (STEP / ISO-10303-21) · #75, G19
//
// Los .ifc son TEXTO (IFC-SPF, ISO-10303-21), no binario.  Para barras estructurales
// (eje de IfcBeam/IfcColumn/IfcMember) basta parsear el texto: no hace falta el motor de
// geometría completo (web-ifc/WASM).  Este parser es JS PURO, sin dependencias, así que
// corre offline en el navegador y en Node (patrón `test_*.mjs`), y encaja en `io/`
// (texto → modelo), igual que los demás adaptadores.
//
// Resultado: un `IfcModel` con un Map  id → { id, type, args }  donde cada `arg` es:
//   { ref:N }        referencia a otra instancia (#N)
//   "texto"          string (ya des-escapado de '' y \X2\…\X0\)
//   123.4            número
//   "BEAM"           enumeración (.BEAM. → 'BEAM')   ── string, se distingue por contexto
//   null             $ (sin valor) o * (derivado)
//   [ … ]            lista/agregado (recursivo)
//   { type, value }  valor tipado (p.ej. IFCBOOLEAN(.T.))
//
// Sólo se interpreta la sección DATA; del HEADER se conserva el esquema (IFC2X3/IFC4).
// AUTÓNOMO (Node + navegador).
// ──────────────────────────────────────────────────────────────────────────────

/** Modelo IFC parseado: instancias indexadas por id + utilidades de des-referencia. */
export class IfcModel {
  constructor() {
    this.entities = new Map();   // id → { id, type, args }
    this.byType = new Map();     // 'IFCBEAM' → [entity, …]
    this.schema = 'IFC4';
    this.header = {};
  }

  /** Instancia por id (acepta number o {ref:N}). `null` si no existe. */
  get(ref) {
    const id = (ref && typeof ref === 'object' && 'ref' in ref) ? ref.ref : ref;
    return this.entities.get(id) || null;
  }

  /** Todas las instancias de un tipo (mayúsculas, p.ej. 'IFCBEAM'). */
  ofType(type) { return this.byType.get(type) || []; }

  /** ¿`arg` es una referencia #N? */
  isRef(a) { return a && typeof a === 'object' && 'ref' in a; }
}

// ── Decodificación de strings IFC (ISO-10303-21 + extensión \X2\ de IFC) ──────────
function decodeIfcString(s) {
  // '' → '  (escape de comilla simple del propio STEP)
  let out = s.replace(/''/g, "'");
  if (out.indexOf('\\') < 0) return out;
  // \S\c  → carácter c con bit alto (Latin-1 +128); poco usado, aproximamos a c.
  out = out.replace(/\\S\\(.)/g, (_, c) => String.fromCharCode(c.charCodeAt(0) + 128));
  // \X\HH → un byte hex
  out = out.replace(/\\X\\([0-9A-Fa-f]{2})/g, (_, h) => String.fromCharCode(parseInt(h, 16)));
  // \X2\HHHH…\X0\ → secuencia de unidades UTF-16 (acentos, etc.)
  out = out.replace(/\\X2\\([0-9A-Fa-f]+)\\X0\\/g, (_, hex) => {
    let r = '';
    for (let i = 0; i + 4 <= hex.length; i += 4) r += String.fromCharCode(parseInt(hex.substr(i, 4), 16));
    return r;
  });
  return out;
}

// ── Parser de la lista de argumentos de una instancia ────────────────────────────
// Recorre `s` desde el índice `pos.i`, asumiendo que `s[pos.i]` es '(' y devuelve el
// array de valores hasta el ')' que cierra (recursivo para listas anidadas).
function parseList(s, pos) {
  const arr = [];
  pos.i++; // saltar '('
  skipWs(s, pos);
  if (s[pos.i] === ')') { pos.i++; return arr; }
  for (;;) {
    arr.push(parseValue(s, pos));
    skipWs(s, pos);
    const c = s[pos.i];
    if (c === ',') { pos.i++; skipWs(s, pos); continue; }
    if (c === ')') { pos.i++; break; }
    // tolerante: si algo raro, cortar para no colgar
    if (c === undefined) break;
    pos.i++;
  }
  return arr;
}

function skipWs(s, pos) {
  for (;;) {
    const c = s[pos.i];
    if (c === ' ' || c === '\t' || c === '\r' || c === '\n') { pos.i++; continue; }
    // comentario /* … */
    if (c === '/' && s[pos.i + 1] === '*') {
      const end = s.indexOf('*/', pos.i + 2);
      pos.i = end < 0 ? s.length : end + 2;
      continue;
    }
    break;
  }
}

function parseValue(s, pos) {
  skipWs(s, pos);
  const c = s[pos.i];
  if (c === '#') {                                   // referencia #N
    let j = pos.i + 1; while (j < s.length && s[j] >= '0' && s[j] <= '9') j++;
    const id = parseInt(s.slice(pos.i + 1, j), 10); pos.i = j;
    return { ref: id };
  }
  if (c === "'") {                                   // string '…'  (con '' internos)
    let j = pos.i + 1, buf = '';
    for (;;) {
      if (j >= s.length) break;
      if (s[j] === "'") {
        if (s[j + 1] === "'") { buf += "''"; j += 2; continue; }  // comilla escapada
        break;
      }
      buf += s[j]; j++;
    }
    pos.i = j + 1;
    return decodeIfcString(buf);
  }
  if (c === '(') return parseList(s, pos);           // lista anidada
  if (c === '$' || c === '*') { pos.i++; return null; } // sin valor / derivado
  if (c === '.') {                                   // enumeración .NAME.
    const end = s.indexOf('.', pos.i + 1);
    const name = s.slice(pos.i + 1, end < 0 ? s.length : end);
    pos.i = end < 0 ? s.length : end + 1;
    return name;
  }
  if (c === '-' || c === '+' || c === '.' || (c >= '0' && c <= '9')) { // número
    let j = pos.i; while (j < s.length && /[0-9+\-.eE]/.test(s[j])) j++;
    const num = parseFloat(s.slice(pos.i, j)); pos.i = j;
    return Number.isFinite(num) ? num : 0;
  }
  // identificador: o bien valor tipado TIPO(args), o palabra suelta
  let j = pos.i; while (j < s.length && /[A-Za-z0-9_]/.test(s[j])) j++;
  const ident = s.slice(pos.i, j); pos.i = j;
  skipWs(s, pos);
  if (s[pos.i] === '(') { const value = parseList(s, pos); return { type: ident.toUpperCase(), value }; }
  return ident;                                      // p.ej. enteros sin punto ya cubiertos arriba
}

// ── Recorre el texto y separa las sentencias `#id=TIPO(...);` ─────────────────────
// Escanea respetando strings y comentarios para encontrar el ';' de cierre de cada
// sentencia (los ';' dentro de '…' o /*…*/ no cuentan).
function* statements(body) {
  let i = 0; const n = body.length;
  while (i < n) {
    // saltar espacios / comentarios entre sentencias
    while (i < n) {
      const c = body[i];
      if (c === ' ' || c === '\t' || c === '\r' || c === '\n') { i++; continue; }
      if (c === '/' && body[i + 1] === '*') { const e = body.indexOf('*/', i + 2); i = e < 0 ? n : e + 2; continue; }
      break;
    }
    if (i >= n) break;
    const start = i;
    let inStr = false;
    while (i < n) {
      const c = body[i];
      if (inStr) {
        if (c === "'") { if (body[i + 1] === "'") { i += 2; continue; } inStr = false; }
        i++; continue;
      }
      if (c === "'") { inStr = true; i++; continue; }
      if (c === '/' && body[i + 1] === '*') { const e = body.indexOf('*/', i + 2); i = e < 0 ? n : e + 2; continue; }
      if (c === ';') break;
      i++;
    }
    const stmt = body.slice(start, i).trim();
    i++; // saltar ';'
    if (stmt) yield stmt;
  }
}

/**
 * Parsea texto IFC-SPF a un `IfcModel`.
 * @param {string} text  contenido del .ifc
 * @returns {IfcModel}
 */
export function parseIFC(text) {
  if (typeof text !== 'string' || !text.length) throw new Error('IFC: archivo vacío');
  if (text.indexOf('ISO-10303-21') < 0 && text.indexOf('IFC') < 0)
    throw new Error('IFC: no parece un archivo IFC (falta cabecera ISO-10303-21)');

  const model = new IfcModel();

  // Esquema desde FILE_SCHEMA(('IFC4')) del HEADER.
  const sch = text.match(/FILE_SCHEMA\s*\(\s*\(\s*'([^']+)'/i);
  if (sch) model.schema = sch[1].toUpperCase();

  // Sólo la sección DATA contiene instancias.
  const dStart = text.indexOf('DATA;');
  const dEnd = text.indexOf('ENDSEC', dStart >= 0 ? dStart : 0);
  const body = text.slice(dStart >= 0 ? dStart + 5 : 0, dEnd >= 0 ? dEnd : text.length);

  for (const stmt of statements(body)) {
    // #id = TIPO ( … )
    const eq = stmt.indexOf('=');
    if (stmt[0] !== '#' || eq < 0) continue;
    const id = parseInt(stmt.slice(1, eq), 10);
    if (!Number.isFinite(id)) continue;
    const rest = stmt.slice(eq + 1).trimStart();
    // nombre del tipo hasta el '('
    const paren = rest.indexOf('(');
    if (paren < 0) continue;
    const type = rest.slice(0, paren).trim().toUpperCase();
    const pos = { i: paren };
    let args;
    try { args = parseList(rest, pos); }
    catch { args = []; }
    const ent = { id, type, args };
    model.entities.set(id, ent);
    if (!model.byType.has(type)) model.byType.set(type, []);
    model.byType.get(type).push(ent);
  }

  if (model.entities.size === 0) throw new Error('IFC: no se encontraron instancias en la sección DATA');
  return model;
}

// ── Unidades: factor de longitud a METROS desde IfcUnitAssignment ─────────────────
const SI_PREFIX = { EXA: 1e18, PETA: 1e15, TERA: 1e12, GIGA: 1e9, MEGA: 1e6, KILO: 1e3, HECTO: 1e2, DECA: 1e1, DECI: 1e-1, CENTI: 1e-2, MILLI: 1e-3, MICRO: 1e-6, NANO: 1e-9 };

/**
 * Factor para convertir longitudes del archivo a METROS y nombre de la unidad.
 * Busca el IFCSIUNIT de tipo LENGTHUNIT (o IFCCONVERSIONBASEDUNIT para pies/pulgadas).
 * @returns {{ factor:number, name:string }}
 */
export function lengthUnit(model) {
  for (const u of model.ofType('IFCSIUNIT')) {
    // IfcSIUnit(Dimensions, UnitType, Prefix, Name)
    const unitType = u.args[1];           // .LENGTHUNIT.
    if (unitType !== 'LENGTHUNIT') continue;
    const prefix = u.args[2];             // .MILLI. | null
    const f = prefix ? (SI_PREFIX[prefix] || 1) : 1;
    return { factor: f, name: (prefix ? prefix.toLowerCase() : '') + 'metre' };
  }
  for (const u of model.ofType('IFCCONVERSIONBASEDUNIT')) {
    // IfcConversionBasedUnit(Dimensions, UnitType, Name, ConversionFactor)
    if (u.args[1] !== 'LENGTHUNIT') continue;
    const name = (u.args[2] || '').toString().toLowerCase();
    const mr = model.get(u.args[3]);      // IfcMeasureWithUnit(ValueComponent, UnitComponent)
    let f = 1;
    if (mr && typeof mr.args[0] === 'object' && mr.args[0].value) f = +mr.args[0].value[0] || 1;
    else if (/foot|feet/.test(name)) f = 0.3048;
    else if (/inch/.test(name)) f = 0.0254;
    return { factor: f, name };
  }
  return { factor: 1, name: 'metre (asumido)' };
}
