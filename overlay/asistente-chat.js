// ──────────────────────────────────────────────────────────────────────────────
// asistente-chat.js — capa overlay UACh: el CHAT del asistente (opennodex).
//
// Reemplaza el asistente de un turno por el agente conversacional de opennodex.
// El agente real (herramientas + motor nodex) corre como servicio HTTP; el Worker
// de PÓRTICO hace de proxy en /api/assistant/agent, así la credencial vive del lado
// servidor y el navegador habla con su propio origen.
//
// PUENTE .ndx↔portico (todo del lado cliente, con la API pública del core):
//   • antes de cada turno  → exportModel(app.model,'ndx')  se manda como `src`.
//   • al volver updatedSrc  → importModel(updatedSrc,'ndx') reemplaza app.model y
//                             repinta el visor (mismo flujo que importar un archivo).
//
// LIMITACIÓN (a propósito, hasta que portico-core actualice su importador .ndx):
// el puente SÓLO transfiere GEOMETRÍA (nodos y barras). El importador .ndx de
// portico-core entiende su propio subconjunto (fix/line…) y NO la gramática que
// escribe opennodex (support/load…), así que los APOYOS y las CARGAS del deck del
// agente NO se importan al modelo de PÓRTICO: hay que reespecificarlos en PÓRTICO.
//
// FALLBACK: si el backend no responde, se ofrece el generador NL→spec del core
// (app.assistantFromText), que es serverless y no necesita el agente.
//
// NO forkea core: usa window.app y las utilidades públicas del submodule.
// ──────────────────────────────────────────────────────────────────────────────
import { exportModel, importModel } from '../vendor/portico-core/js/io/index.js?v=2';
import { esc } from '../vendor/portico-core/js/utils/escape.js?v=2';

const LS_MODEL   = 'portico_asis_model';   // último modelo elegido
const LS_MODE    = 'portico_asis_mode';    // último modo elegido
const LS_CFG     = 'portico_asis_cfg';     // { endpoint, token } del proveedor

// Config de proveedor. Vacío ⇒ servicio gratis por defecto vía el Worker
// (/api/assistant/*), con la clave del lado servidor. Si el usuario define un
// `endpoint`, el chat habla DIRECTO con ese backend opennodex (que expone
// /api/agent y /api/models) usando su `token` como Authorization.
function loadCfg() {
  try { const c = JSON.parse(localStorage.getItem(LS_CFG) || '{}'); return { endpoint: c.endpoint || '', token: c.token || '' }; }
  catch { return { endpoint: '', token: '' }; }
}
let cfg = loadCfg();
const trimBase = (u) => String(u || '').trim().replace(/\/+$/, '');
const agentUrl  = () => (cfg.endpoint ? `${trimBase(cfg.endpoint)}/api/agent`  : '/api/assistant/agent');
const agentStreamUrl = () => (cfg.endpoint ? `${trimBase(cfg.endpoint)}/api/agent/stream` : '/api/assistant/agent/stream');
const modelsUrl = () => (cfg.endpoint ? `${trimBase(cfg.endpoint)}/api/models` : '/api/assistant/models');
const authHeaders = () => (cfg.endpoint && cfg.token ? { Authorization: `Bearer ${cfg.token}` } : {});

// El generador NL→spec (fallback) del core resuelve su endpoint desde este
// localStorage. Como el Worker sirve la PWA en el MISMO origen, dejamos un default
// razonable para que el fallback funcione sin configuración manual.
if (!localStorage.getItem('portico_n8n_endpoint')) {
  try { localStorage.setItem('portico_n8n_endpoint', location.origin + '/api/assistant'); } catch { /* modo privado */ }
}

let conversation = [];   // historial hilado que devuelve el agente (memoria entre turnos)
let sending = false;

// ── Markdown mínimo (la respuesta del agente es GFM): fences, inline code,
//    negritas/itálicas, encabezados, listas y saltos de línea. Todo escapado. ──
function renderMarkdown(src) {
  const fences = [];
  const OPEN = "\uF8F0", CLOSE = "\uF8F1";   // centinelas (uso privado) fuera del texto normal
  let s = String(src || "").replace(/```(\w+)?\n?([\s\S]*?)```/g, (_m, lang, code) => {
    fences.push(`<pre class="asis-code"><code>${esc(code.replace(/\n$/, ""))}</code></pre>`);
    return `${OPEN}${fences.length - 1}${CLOSE}`;
  });
  s = esc(s)
    .replace(/`([^`]+)`/g, "<code>$1</code>")
    .replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>")
    .replace(/(^|[^*])\*([^*\n]+)\*/g, "$1<em>$2</em>")
    .replace(/^######\s?(.*)$/gm, "<h6>$1</h6>")
    .replace(/^#{1,5}\s?(.*)$/gm, "<h4>$1</h4>")
    .replace(/^\s*[-*]\s+(.*)$/gm, "<li>$1</li>");
  s = s.replace(/(<li>[\s\S]*?<\/li>)/g, (m) => `<ul>${m}</ul>`).replace(/<\/ul>\s*<ul>/g, "");
  s = s.replace(/\n{2,}/g, "<br><br>").replace(/\n/g, "<br>");
  return s.replace(new RegExp(OPEN + "(\\d+)" + CLOSE, "g"), (_m, i) => fences[+i]);
}

const $ = (id) => document.getElementById(id);

function bubble(role, html, cls = '') {
  const log = $('asis-chat-log'); if (!log) return null;
  const div = document.createElement('div');
  div.className = `asis-msg asis-${role} ${cls}`.trim();
  div.innerHTML = html;
  log.appendChild(div);
  log.scrollTop = log.scrollHeight;
  return div;
}

// Semilla para "modelo nuevo": el backend SÓLO devuelve updatedSrc cuando recibe un
// `src` no vacío (escribe ese deck y lo relee). Con un deck vacío no hay archivo que
// releer, así que mandamos siempre esta semilla y el agente construye sobre ella.
const NEW_DECK_SEED = '// PÓRTICO: modelo nuevo — describe la estructura y la construiré.\n';

// Exporta el modelo actual a .ndx (semilla si no hay nada construido todavía).
function currentDeck() {
  const app = window.app;
  if (!app || !app.model || !app.model.nodes || app.model.nodes.size === 0) return NEW_DECK_SEED;
  try { return exportModel(app.model, 'ndx').text || NEW_DECK_SEED; }
  catch (e) { console.warn('[asistente] no se pudo exportar el modelo a .ndx:', e); return NEW_DECK_SEED; }
}

// Si el agente no escribió el deck (updatedSrc null) pero propuso uno en su texto
// dentro de un bloque ```ndx, lo recuperamos para ofrecer cargarlo (como el UI de opennodex).
function extractDeck(text) {
  const blocks = [...String(text || '').matchAll(/```(?:ndx)?\s*\n([\s\S]*?)```/g)].map((x) => x[1]);
  return blocks.find((b) => /(^|\n)\s*(model|node|beam|column|section|material|solve)\b/.test(b)) || null;
}

// Ofrece (con un botón, porque reemplaza el modelo) cargar un deck propuesto.
function offerLoadDeck(src) {
  const div = bubble('tool', 'El agente propuso un deck (no lo aplicó). ¿Cargarlo en PÓRTICO?');
  if (!div) return;
  const btn = document.createElement('button');
  btn.className = 'btn-add';
  btn.style.cssText = 'margin-top:6px;background:rgba(14,127,192,0.12);color:var(--accent);border-color:var(--accent)';
  btn.textContent = '⤓ Cargar deck propuesto';
  btn.onclick = () => {
    btn.disabled = true;
    try { const { info, warnings } = applyDeck(src); bubble('tool', `✎ ${esc(info)}`); if (warnings.length) bubble('tool', `⚠ ${esc(warnings.slice(0, 8).join(' · '))}`, 'warn'); noteGeometryOnly(src); }
    catch (e) { bubble('tool', `⚠ no se pudo importar el deck: ${esc(String(e.message || e))}`, 'warn'); }
  };
  div.appendChild(document.createElement('br'));
  div.appendChild(btn);
}

// Reimporta el .ndx que devolvió el agente y repinta el visor (mismo flujo que el
// importador de archivos del core). Devuelve un aviso legible o null.
function applyDeck(src) {
  const app = window.app;
  const { model, warnings } = importModel(src, 'ndx');
  if (typeof app.snapshot === 'function') app.snapshot();   // integrar en undo/redo
  app.model = model;
  app.viewport?.renderModel?.(model);
  app.panel?.showNothing?.();
  app.panel?.refresh?.(model);
  app.markDirty?.();
  app._updateStats?.();
  app.viewport?.zoomExtents?.();
  const n = model.nodes?.size ?? 0, e = model.elements?.size ?? 0;
  return { info: `modelo actualizado: ${n} nodos, ${e} barras`, warnings: warnings || [] };
}

// El deck del agente (gramática nodex) trae apoyos/cargas que el importador de
// portico-core NO reconoce → se pierden. Detectarlos para avisar de forma explícita.
function droppedBoundaryLoads(src) {
  const s = String(src || '');
  return /(^|\n)\s*support\b/i.test(s) || /(^|\n)\s*load\b[^\n]*\bon\b/i.test(s);
}

// Aviso explícito y persistente: el puente sólo transfiere geometría.
function noteGeometryOnly(src) {
  const extra = droppedBoundaryLoads(src)
    ? ' El deck traía <b>apoyos y/o cargas</b> que no se importaron.'
    : '';
  bubble('tool', `ℹ El puente sólo transfiere la <b>geometría</b> (nodos y barras).${extra} Reespecifica apoyos y cargas en PÓRTICO.`, 'warn');
}

// Ofrece el generador NL→spec del core cuando el agente no está disponible.
function offerFallback(rawPrompt) {
  const div = bubble('tool', 'El agente no está disponible. Puedes generar un modelo con el <b>modo básico</b> (sin IA conversacional):');
  if (!div) return;
  const btn = document.createElement('button');
  btn.className = 'btn-add';
  btn.style.cssText = 'margin-top:6px;background:rgba(14,127,192,0.12);color:var(--accent);border-color:var(--accent)';
  btn.textContent = '✦ Generar con el modo básico';
  btn.onclick = () => { btn.disabled = true; window.app?.assistantFromText?.(rawPrompt); };
  div.appendChild(document.createElement('br'));
  div.appendChild(btn);
}

// Aplica el resultado final del agente (texto + deck) a la UI/modelo.
function renderResult(data, mode, prevSrc) {
  conversation = Array.isArray(data.history) ? data.history : conversation;   // hilar memoria
  const footer = `<div class="asis-footer">[${esc(data.mode || mode)}${data.runId ? ` · run ${esc(String(data.runId))}` : ''}${data.model ? ` · ${esc(String(data.model))}` : ''}]</div>`;
  // Si veníamos transmitiendo texto en vivo, reusar esa burbuja; si no, crear una.
  if (data._streamEl) { data._streamEl.classList.add('md'); data._streamEl.innerHTML = renderMarkdown(data.text || data._streamed || '(sin respuesta)') + footer; }
  else bubble('agent', renderMarkdown(data.text || '(sin respuesta)') + footer, 'md');

  if (typeof data.updatedSrc === 'string' && data.updatedSrc && data.updatedSrc !== prevSrc) {
    try {
      const { info, warnings } = applyDeck(data.updatedSrc);
      bubble('tool', `✎ ${esc(info)}`);
      if (warnings.length) bubble('tool', `⚠ ${esc(warnings.slice(0, 8).join(' · '))}`, 'warn');
      noteGeometryOnly(data.updatedSrc);
    } catch (e) {
      bubble('tool', `⚠ el agente propuso un deck que PÓRTICO no pudo importar (${esc(String(e.message || e))}). El .ndx quedó en la respuesta.`, 'warn');
    }
  } else {
    const proposed = extractDeck(data.text);
    if (proposed) offerLoadDeck(proposed);
  }
}

async function send() {
  if (sending) return;
  const ta = $('asis-chat-input');
  const raw = ta?.value.trim();
  if (!raw) return;
  const mode  = $('asis-chat-mode')?.value || 'apply';
  const model = $('asis-chat-model')?.value || '';
  ta.value = '';
  sending = true;
  const sendBtn = $('asis-chat-send'); if (sendBtn) sendBtn.disabled = true;

  bubble('user', `<span class="asis-tag">[${esc(mode)}]</span> ${esc(raw)}`);
  const stat = $('asis-chat-stat');
  const log = $('asis-chat-log');
  const t0 = Date.now();
  const frames = '⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'; let fi = 0;
  let label = 'trabajando…', liveTok = 0;
  if (stat) stat.classList.add('busy');
  const tick = stat ? setInterval(() => {
    stat.textContent = `${frames[fi = (fi + 1) % frames.length]} ${label} · ${((Date.now() - t0) / 1000).toFixed(1)}s${liveTok ? ` · ~${liveTok} tok` : ''}`;
  }, 100) : null;

  const body = { prompt: raw, src: currentDeck(), mode, history: conversation };
  if (model) body.model = model;

  try {
    // Streaming (NDJSON): muestra rondas, herramientas y texto EN VIVO — clave en
    // modo aplicar, que puede tardar >1 min. Si el stream falla, cae al no-stream.
    let result = null, streamEl = null, streamed = '';
    let res;
    try {
      res = await fetch(agentStreamUrl(), { method: 'POST', headers: { 'Content-Type': 'application/json', ...authHeaders() }, body: JSON.stringify(body) });
    } catch (netErr) { res = null; }

    if (res && res.ok && res.body) {
      const reader = res.body.getReader();
      const dec = new TextDecoder();
      let buf = '';
      for (;;) {
        const { value, done } = await reader.read();
        if (done) break;
        buf += dec.decode(value, { stream: true });
        let nl;
        while ((nl = buf.indexOf('\n')) >= 0) {
          const line = buf.slice(0, nl).trim(); buf = buf.slice(nl + 1);
          if (!line) continue;
          let obj; try { obj = JSON.parse(line); } catch { continue; }
          if (obj.ev) {
            const e = obj.ev;
            if (e.type === 'round') label = `trabajando… (ronda ${e.round}/${e.of})`;
            else if (e.type === 'usage') liveTok = (e.usage.input || 0) + (e.usage.output || 0);
            else if (e.type === 'text_delta') { if (!streamEl) streamEl = bubble('agent', ''); streamed += e.text; streamEl.textContent = streamed; if (log) log.scrollTop = log.scrollHeight; }
            else if (e.type === 'tool') { label = `→ ${e.name}…`; bubble('tool', `→ ${esc(e.name)}`); }
            else if (e.type === 'tool_error') bubble('tool', `✗ ${esc(e.name)}: ${esc(String(e.text || ''))}`, 'warn');
          } else if (obj.error) { throw new Error(obj.error); }
          else if (obj.done) result = obj.done;
        }
      }
      if (!result) throw new Error('el stream terminó sin resultado');
    } else {
      // Respaldo sin streaming
      const r = await fetch(agentUrl(), { method: 'POST', headers: { 'Content-Type': 'application/json', ...authHeaders() }, body: JSON.stringify(body) });
      const data = await r.json().catch(() => ({}));
      if (!r.ok || data.error) throw new Error(data.error || `HTTP ${r.status}`);
      result = data;
    }

    if (tick) clearInterval(tick);
    const tok = result.usage ? (result.usage.input || 0) + (result.usage.output || 0) : liveTok;
    if (stat) { stat.classList.remove('busy'); stat.textContent = `listo · ${((Date.now() - t0) / 1000).toFixed(1)}s${tok ? ` · ${tok} tok` : ''}`; }
    result._streamEl = streamEl; result._streamed = streamed;
    renderResult(result, mode, body.src);
  } catch (e) {
    if (tick) clearInterval(tick);
    if (stat) { stat.classList.remove('busy'); stat.textContent = 'error'; }
    bubble('agent', `<span class="asis-err">⚠ ${esc(String(e.message || e))}</span>`, 'md');
    offerFallback(raw);
  } finally {
    sending = false;
    if (sendBtn) sendBtn.disabled = false;
    ta?.focus();
  }
}

async function loadModels() {
  const sel = $('asis-chat-model'); if (!sel) return;
  try {
    const r = await fetch(modelsUrl(), { headers: authHeaders() });
    const data = await r.json();
    const models = Array.isArray(data.models) ? data.models : [];
    if (!models.length) { sel.innerHTML = '<option value="">(modelo por defecto)</option>'; return; }
    const last = localStorage.getItem(LS_MODEL) || (data.fixed || '');
    sel.innerHTML = models.map((m) => `<option value="${esc(m)}"${m === last ? ' selected' : ''}>${esc(m)}</option>`).join('');
    if (data.fixed) sel.disabled = true;   // el host fijó un único modelo
  } catch {
    sel.innerHTML = '<option value="">(modelo por defecto)</option>';
  }
}

function wire() {
  const sendBtn = $('asis-chat-send'), ta = $('asis-chat-input');
  if (!sendBtn || !ta) return false;    // el panel aún no está en el DOM
  sendBtn.addEventListener('click', send);
  ta.addEventListener('keydown', (e) => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); send(); } });
  $('asis-chat-new')?.addEventListener('click', () => {
    conversation = [];
    const log = $('asis-chat-log'); if (log) log.innerHTML = '';
    const stat = $('asis-chat-stat'); if (stat) stat.textContent = '';
    bubble('tool', 'Conversación nueva.');
  });
  const modeSel = $('asis-chat-mode');
  if (modeSel) {
    modeSel.value = localStorage.getItem(LS_MODE) || 'apply';
    modeSel.addEventListener('change', () => localStorage.setItem(LS_MODE, modeSel.value));
  }
  $('asis-chat-model')?.addEventListener('change', (e) => localStorage.setItem(LS_MODEL, e.target.value));

  // ── Proveedor: endpoint + token (o servicio gratis por defecto) ──
  const epIn = $('asis-cfg-endpoint'), tkIn = $('asis-cfg-token');
  if (epIn) epIn.value = cfg.endpoint;
  if (tkIn) tkIn.value = cfg.token;
  $('asis-cfg-save')?.addEventListener('click', () => {
    cfg = { endpoint: trimBase(epIn?.value), token: (tkIn?.value || '').trim() };
    localStorage.setItem(LS_CFG, JSON.stringify(cfg));
    const sel = $('asis-chat-model'); if (sel) sel.disabled = false;
    bubble('tool', cfg.endpoint ? `Proveedor: ${esc(cfg.endpoint)}${cfg.token ? ' (con token)' : ''}.` : 'Proveedor: servicio gratis por defecto.');
    loadModels();
  });
  $('asis-cfg-clear')?.addEventListener('click', () => {
    cfg = { endpoint: '', token: '' };
    localStorage.removeItem(LS_CFG);
    if (epIn) epIn.value = ''; if (tkIn) tkIn.value = '';
    const sel = $('asis-chat-model'); if (sel) sel.disabled = false;
    bubble('tool', 'Proveedor: servicio gratis por defecto.');
    loadModels();
  });

  loadModels();
  return true;
}

// El panel vive en el index.html del overlay; window.app lo crea el core en
// DOMContentLoaded. Esperamos a ambos (mismo patrón que abrirAsistente del index).
function boot() {
  if (wire()) return;
  let tries = 0;
  const iv = setInterval(() => { if (wire() || ++tries > 40) clearInterval(iv); }, 100);
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
else boot();
