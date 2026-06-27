// docs/refresh_ui_texts.mjs — mantiene docs/textos_ui.csv al día SIN perder la curación.
//   1) Re-localiza cada fila existente en su archivo (refresca el número de línea).
//   2) Anexa los textos de UI NUEVOS que aún no estén en el CSV.
// Conserva la 4.ª columna `texto_propuesto` (tu revisión). Re-ejecutable.
//
//   node docs/refresh_ui_texts.mjs
import { readFileSync, writeFileSync } from 'node:fs';

const CSV = 'docs/textos_ui.csv';
// Archivos «tipo HTML»: se extrae el texto visible entre etiquetas (poco ruido).
const HTML_LIKE = ['index.html', 'js/ui/ifcImportDialog.js'];
// Archivos de código: sólo mensajes + title= + placeholder= (evita ruido de plantillas).
const CODE = ['js/app.js', 'js/ui/properties.js', 'js/ui/viewport.js', 'js/ui/menu.js'];

const ENT = { '&amp;': '&', '&lt;': '<', '&gt;': '>', '&quot;': '"', '&#39;': "'", '&nbsp;': ' ', '&times;': '×', '&middot;': '·' };
const decode = s => s.replace(/&[a-z#0-9]+;/gi, m => ENT[m] ?? m);
const norm = s => decode(s).replace(/\$\{[^}]*\}/g, '…').replace(/\s+/g, ' ').trim();
const hasLetter = s => /[A-Za-zÀ-ÿ]/.test(s);
const meaningful = s => s && hasLetter(s) && s !== '…' && !/^[…\s.,:;|/×·–—\-+]+$/.test(s)
  && !/^[a-zà-ÿ]{1,2}$/.test(s);   // descarta tokens sueltos de 1–2 letras minúsculas (ruido), conserva XY/Hz/3D…

const fileCache = {};
const linesOf = f => (fileCache[f] ??= (() => { try { return readFileSync(f, 'utf8').split(/\r?\n/); } catch { return null; } })());

// localiza la 1.ª línea del archivo que contiene el texto (usa el fragmento literal más
// largo, partiendo por «…» que representa interpolaciones)
function locate(file, text) {
  const src = linesOf(file); if (!src) return null;
  let frags = text.split('…').map(s => s.trim()).filter(s => s.length >= 4 && hasLetter(s));
  if (!frags.length) frags = [text.trim()];
  frags.sort((a, b) => b.length - a.length);
  for (let i = 0; i < src.length; i++) for (const f of frags) if (src[i].includes(f)) return i + 1;
  return null;
}

// ── parse CSV (todos los campos entre comillas, "" escapa comilla) ──────────────
function parseCSV(txt) {
  const out = []; const lines = txt.split(/\r?\n/);
  for (let i = 1; i < lines.length; i++) {
    const ln = lines[i]; if (!ln.trim()) continue;
    const f = []; let j = 0;
    while (j < ln.length) {
      if (ln[j] !== '"') { j++; continue; }
      j++; let v = '';
      while (j < ln.length) { if (ln[j] === '"') { if (ln[j + 1] === '"') { v += '"'; j += 2; continue; } j++; break; } v += ln[j++]; }
      f.push(v);
      if (ln[j] === ',') j++;
    }
    if (f.length >= 3) out.push({ ubic: f[0], cat: f[1], text: f[2], prop: f[3] || '' });
  }
  return out;
}

// texto VISIBLE combinado de una línea: se quitan las etiquetas (los hijos inline como
// <kbd> se funden con un espacio, igual que el texto que ve el usuario). Una fila por línea.
function htmlText(line) {
  if (!line.includes('<')) return '';
  // sólo líneas que tengan algún cierre/contenido entre etiquetas (evita aperturas sueltas)
  if (!/>[^<>]*[A-Za-zÀ-ÿ][^<>]*</.test(line)) return '';
  const t = line.replace(/<[^>]*>/g, ' ');
  // descarta artefactos de etiquetas multilínea (atributos CSS que quedaron sin tag de apertura)
  if (/["']>|:var\(|style=|white-space:|cursor:|;\s*(font|color|padding|margin|background|flex|width|height|border)|px[;"]|\$\{/.test(t)) return '';
  return t;
}
// mensajes JS
const MSG = [/\btoast\(\s*(['"`])([^'"`]+?)\1/g, /\b_confirm\(\s*(['"`])([^'"`]+?)\1/g,
  /\b_alert\(\s*(['"`])([^'"`]+?)\1/g, /\b_promptModal\(\s*(['"`])([^'"`]+?)\1/g,
  /\b_prompt\(\s*(['"`])([^'"`]+?)\1/g, /modal-title'\)\.textContent\s*=\s*(['"`])([^'"`]+?)\1/g];

// ── 1) refrescar líneas de las filas existentes ─────────────────────────────────
let rows = parseCSV(readFileSync(CSV, 'utf8').replace(/^﻿/, ''));   // ignora BOM si lo hay
// poda filas-ruido (texto no significativo) SÓLO si no tienen propuesta escrita
rows = rows.filter(r => r.prop.trim() !== '' || meaningful(r.text));
let refreshed = 0;
for (const r of rows) {
  const file = r.ubic.replace(/\s*\(L\d+\)\s*$/, '');
  const ln = locate(file, r.text);
  if (ln) { const nu = `${file} (L${ln})`; if (nu !== r.ubic) refreshed++; r.ubic = nu; }
}

// índice de lo ya presente (por categoría + texto normalizado)
const have = new Set(rows.map(r => `${r.cat}|${norm(r.text)}`));
const add = (file, ln, cat, raw) => {
  const text = norm(raw); if (!meaningful(text)) return;
  const k = `${cat}|${text}`; if (have.has(k)) return; have.add(k);
  rows.push({ ubic: `${file} (L${ln})`, cat, text, prop: '', _new: true });
};

// ── 2) anexar textos NUEVOS ──────────────────────────────────────────────────────
for (const file of HTML_LIKE) {
  const src = linesOf(file); if (!src) continue;
  let inSvg = false;
  src.forEach((line, i) => {
    if (/<svg\b/i.test(line)) inSvg = true; const wasSvg = inSvg; if (/<\/svg>/i.test(line)) inSvg = false; if (wasSvg) return;
    for (const m of line.matchAll(/\btitle="([^"]+)"/g)) add(file, i + 1, 'tooltip', m[1]);
    for (const m of line.matchAll(/\bplaceholder="([^"]+)"/g)) add(file, i + 1, 'placeholder', m[1]);
    { const t = htmlText(line); if (t) add(file, i + 1, "etiqueta/HTML", t); }
  });
}
for (const file of CODE) {
  const src = linesOf(file); if (!src) continue;
  src.forEach((line, i) => {
    if (line.trim().startsWith('//')) return;
    for (const m of line.matchAll(/\btitle="([^"]+)"/g)) add(file, i + 1, 'tooltip', m[1]);
    for (const m of line.matchAll(/\bplaceholder="([^"]+)"/g)) add(file, i + 1, 'placeholder', m[1]);
    for (const re of MSG) for (const m of line.matchAll(re)) add(file, i + 1, 'mensaje/JS', m[2]);
  });
}
// etiquetas de plantilla NUEVAS en properties.js (no son title/placeholder/mensaje)
for (const t of ['Irx (ton·m²)', 'Iry (ton·m²)', 'Irz (ton·m²)', 'Unilateral (a los GDL con resorte)',
  'Bilateral (normal)', 'Solo-compresión (uplift)', 'Solo-tracción', 'k incl. (kN/m)', '↗ Aplicar resorte inclinado']) {
  const ln = locate('js/ui/properties.js', t); if (ln) add('js/ui/properties.js', ln, 'etiqueta/HTML', t);
}

// ── orden: por archivo (orden de aparición) y luego por línea ────────────────────
const order = ['index.html', 'js/app.js', 'js/ui/properties.js', 'js/ui/viewport.js', 'js/ui/menu.js', 'js/ui/ifcImportDialog.js'];
const fileOf = r => r.ubic.replace(/\s*\(L\d+\)\s*$/, '');
const lineOf = r => { const m = r.ubic.match(/\(L(\d+)\)/); return m ? +m[1] : 0; };
rows.sort((a, b) => { const fa = order.indexOf(fileOf(a)), fb = order.indexOf(fileOf(b)); return (fa - fb) || (lineOf(a) - lineOf(b)); });

const q = s => '"' + String(s).replace(/"/g, '""') + '"';
const csv = ['ubicacion,categoria,texto_actual,texto_propuesto',
  ...rows.map(r => [q(r.ubic), q(r.cat), q(r.text), q(r.prop)].join(','))].join('\r\n') + '\r\n';
// UTF-8 CON BOM (﻿) + saltos CRLF → Excel en Windows lo abre con los tildes correctos
// (sin BOM lo lee como ANSI/Latin-1 y los acentos/símbolos salen ilegibles).
writeFileSync(CSV, '﻿' + csv);

const nuevos = rows.filter(r => r._new).length;
console.log(`textos_ui.csv: ${rows.length} filas · líneas refrescadas: ${refreshed} · filas nuevas: ${nuevos}`);
const byFile = {}; for (const r of rows.filter(r => r._new)) { const f = fileOf(r); byFile[f] = (byFile[f] || 0) + 1; }
console.log('Nuevas por archivo:', byFile);
