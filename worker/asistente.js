// ──────────────────────────────────────────────────────────────────────────────
// Cloudflare Worker (capa UACh) — sirve la PWA (assets) y expone la API del asistente.
//   POST /api/assistant            body { message } → LLM → spec  (RAG few-shot)
//                                  body { spec }    → devuelve el spec tal cual
//   POST /api/assistant/modificar  body { message, model?, selection? } → { ops }
//   POST /api/assistant/feedback   body { id, comentario? } → marca 'incorrecto'
//   GET  /api/assistant/log?token= → registro KV (revisión semanal)
//   resto de rutas                 → assets estáticos (PÓRTICO: overlay + vendor/portico-core)
//
// El LLM SÓLO traduce lenguaje natural → spec estructurado (spec.schema.json del
// core); el modelo lo construye el CLIENTE de forma determinista (generator.js del
// core). Por eso el worker NO importa el generador ni genera nada server-side.
//
// La API key vive como SECRETO del Worker (env.OPENAI_API_KEY / env.OPENROUTER_API_KEY),
// nunca en el código ni en el navegador:  npx wrangler secret put OPENROUTER_API_KEY
//
// Compatibilidad de rutas: el core pide al endpoint que el usuario configura
// (p.ej. .../api/assistant) y le añade /modificar, /feedback. Aceptamos además el
// segmento histórico /api/asistente (español) como alias.
// ──────────────────────────────────────────────────────────────────────────────

// Modelos gratis de OpenRouter en cascada: si uno está rate-limited (429) o caído
// (5xx), se prueba el siguiente. env.OPENROUTER_MODEL fuerza uno solo.
const MODELS_FREE = [
  'deepseek/deepseek-chat-v3-0324:free',
  'meta-llama/llama-3.3-70b-instruct:free',
  'google/gemini-2.0-flash-exp:free',
  'qwen/qwen-2.5-72b-instruct:free',
  'meta-llama/llama-3.1-8b-instruct:free',
];

// SYSTEM prompt — produce el "spec" estructurado del core (assistant/spec.schema.json,
// campos en INGLÉS). El LLM rellena campos acotados; NO genera geometría nodo a nodo
// ni inventa valores de ingeniería. Responde SOLO con el JSON del spec.
const SYSTEM = `You convert a natural-language description of a structure into a SPEC JSON for PORTICO. Reply with ONLY the spec JSON, no prose, no markdown. Top-level fields: project, mode("2D"|"3D"), typology, location{city,latitude_deg,altitude_masl,exposure("B"|"C"|"D")}, geometry{levels:[{height_m,use_class?,live_load_kN_m2?,dead_extra_kN_m2?}],spans_x?,spans_y?,base_plan?{Lx_m,Ly_m},top_plan?,roof_slope_deg?}, sections{material,beams,columns}, base_support("fixed"|"pinned"), rigid_diaphragm, loads{dead_extra_kN_m2,use_class,wind_enclosure,snow,wind,seismic}, seismic{zone(1|2|3),soil("A".."E"),category("I".."IV"),R}. Required: mode and geometry. Omit anything not mentioned; do not invent engineering values.
spans_x/spans_y: list of spans [3,3,3,4] OR uniform {count,span_m} (e.g. "4 bays of 3 m in X" -> spans_x:{count:4,span_m:3}). Each level may have its own height and use (e.g. level 1 classrooms, level 3 light storage) -> levels[k].use_class.
sections.beams/columns: a steel profile by name ('IPE300','HEB200') OR a rectangular concrete section {b_cm,h_cm} (e.g. RC beam 20x40cm -> {b_cm:20,h_cm:40}; column 30x30 -> {b_cm:30,h_cm:30}). material: steel (S235/S275/S355/A630-420H) or concrete grade G (NCh170:2016) G20/G25/G30/G40 (cylinder fc in MPa -> G{fc}, e.g. fc=30 -> "G30"). If sections are concrete, material must be a Gxx. base_support: "fixed" (empotrado) or "pinned" (rotulado).
TYPOLOGY (field "typology"): default "frame" (steel/concrete columns+beams).
- "timber_walls": light stud framing (timber OR steel/Metalcon). Do NOT use sections.beams/columns; instead add stud_walls{nominal_size,spacing_m,perimeter,diagonals?,diagonal_bays?,interior:[{level,dir("X"|"Y"),pos_m,openings:[{type("door"|"window"),width_m,opening_height_m?,center_m?}]}]} and floors{nominal_size,spacing_m,dir("X"|"Y")}. Stud sizes IN INCHES as string: "2x4","2x6","2x8","2x10". Plan in geometry.base_plan{Lx_m,Ly_m}. material:"Pino Radiata" (timber) or "S275"/"acero" (steel framing). Roof: flat by default; for gabled Warren-truss roof add roof{type:"truss",slope_pct,spacing_m,chord_size,diagonal_size}.
- "truss": standalone gabled/Warren-Pratt-Howe truss → truss{span_m,slope_pct OR ridge_height_m,n_panels,spacing_m,chord_size,diagonal_size,truss_type("warren"|"pratt"|"howe")}. mode "2D" by default; if the user explicitly asks 3D, set mode "3D".
- "warehouse": industrial shed → warehouse{span_m,length_m,column_height_m,frame_spacing_m,slope_pct,truss_type,column_section?,chord_size?,diagonal_size?}.
- "torre": transmission tower / space lattice → tower{height_m,base_m,top_m,panels,bracing("X"),pinned(true=truss|false=rigid),chord_profile?,diagonal_profile?,crossarms:[{z_m,length_m,vertical_load_kN?,transverse_load_kN?}]}.
- "bridge": bridge → bridge{length_m,width_m,pier_height_m,pier_span_m OR n_piers,type}. type: "deck" (2 side girders+transverse), "central_beam" (ONE central beam+transverse+piers), "arch" (deck-arch w/ spandrel posts, Salginatobel), "tied_arch" (bowstring, vertical hangers), "network" (bowstring, crossed inclined hangers), "cable_stayed" (central pylon + fan stays), "suspension" (parabolic cable + towers + hangers). For arch/cable use rise_m, n_hangers, pylon_height_m/tower_height_m, arch_section/pylon_section. Concrete sections always as {b_cm,h_cm}.
- "primitives": FREE structure for ANY case without a template (custom masts, grillages, N-girder bridges). Use ONLY when no template fits or the user gives explicit coordinates/counts. default_material + elements[] + supports[]. element:{type("bar"|"beam"|"column"|"repeated_beams"),from:[x,y,z],to:[x,y,z] (m),section:{b_cm,h_cm}|"IPE300"|"2x4",material?,n?,load_kN_m?}. repeated_beams repeats the base bar along step_dir("X"|"Y"|"Z") every step, n_repeat (or to_coord). support:{en:[[x,y,z],...] or z:0 (whole level),type("fixed"|"pinned"|"roller")}. Bars sharing a point are joined (same node).
All lattices (truss/roof/bridge/warehouse) accept truss_type warren|pratt|howe.
EXAMPLE timber ("2-story house 3m, 8x6 plan, 2x4 studs @40cm, 2x8 floor joists @60cm, residential, radiata pine, central partition with an 80cm door"):
{"mode":"3D","typology":"timber_walls","sections":{"material":"Pino Radiata"},"geometry":{"base_plan":{"Lx_m":8,"Ly_m":6},"levels":[{"height_m":3},{"height_m":3}]},"stud_walls":{"nominal_size":"2x4","spacing_m":0.4,"perimeter":true,"interior":[{"level":1,"dir":"Y","pos_m":4.0,"openings":[{"type":"door","width_m":0.8,"opening_height_m":2.0,"center_m":3.0}]}]},"floors":{"nominal_size":"2x8","spacing_m":0.6,"dir":"X"},"loads":{"use_class":"Habitacionales/Viviendas"}}
EXAMPLE bridge ("100m central-beam bridge, RC fc=50, 50x200cm central beam, 30x80cm transverse @2m, 10m wide, 100x100cm piers @20m, 500 kN/m on the transverse beams"):
{"mode":"3D","typology":"bridge","sections":{"material":"G50"},"bridge":{"type":"central_beam","length_m":100,"width_m":10,"pier_height_m":5,"pier_span_m":20,"central_beam_section":{"b_cm":50,"h_cm":200},"transverse_section":{"b_cm":30,"h_cm":80},"pier_section":{"b_cm":100,"h_cm":100},"transverse_spacing_m":2,"transverse_load_kN_m":500}}`;

// Prompt para MODIFICAR un modelo ya construido: orden NL → lista de OPERACIONES
// (las ejecuta el cliente de forma determinista, ver js/model/model_ops.js del core).
const MOD_SYSTEM = `You convert a modification ORDER on an ALREADY-BUILT PORTICO structural model into a JSON list of OPERATIONS. Reply with ONLY JSON {"ops":[...]}, no prose, no markdown. Each operation has "op". Valid types:
- {"op":"add_load","target","case?","dir?","w","w2?"}: distributed load (kN/m). target: "selection" (what the user selected), "all_beams" (all beams/horizontals), "columns", "all", or a list of element ids. dir: "gravity" (default, downward -Z), "globalX","globalY","localY","localZ". w = intensity; w2 = intensity at the j-end (TRAPEZOIDAL); omit w2 if uniform. case = load-case name (default the live-load case L).
- {"op":"add_story","height","copies?"}: append story(ies) ON TOP replicating the top level a height (m). copies = number of stories (default 1).
- {"op":"add_bay","dir","span","copies?"}: append lateral bay(s) extending the plan. dir: "x" or "y". span = bay length (m). copies default 1.
- {"op":"set_modifiers","target","mods":{"A?","Iy?","Iz?","J?"}}: stiffness factors (cracked section, etc.) on target elements.
- {"op":"set_mass","target","mass":{"mx?","my?","mz?"}}: nodal mass (ton) on NODES. node target: "selection","all" or a list of node ids.
Use the model summary (levels_z, axes_x, axes_y, cases, sections, bbox, units) and the selection to choose coherent target/values. If the order mentions "selection"/"selected" use target "selection". Do not invent operations outside this list; if something can't be expressed, omit it. Return {"ops":[]} if nothing applies.
Examples:
"add 20 kN/m live load to all beams" -> {"ops":[{"op":"add_load","target":"all_beams","w":20}]}
"append a 3 m story on top" -> {"ops":[{"op":"add_story","height":3}]}
"add a 5 m bay to the right in X" -> {"ops":[{"op":"add_bay","dir":"x","span":5}]}
"triangular load 0 to 10 kN/m on the selection" -> {"ops":[{"op":"add_load","target":"selection","w":0,"w2":10}]}
"apply Iz 0.5 stiffness modifier to the selection" -> {"ops":[{"op":"set_modifiers","target":"selection","mods":{"Iz":0.5}}]}
"put 2 ton horizontal mass on all nodes" -> {"ops":[{"op":"set_mass","target":"all","mass":{"mx":2,"my":2}}]}`;

const json = (obj, status = 200, extra = {}) =>
  new Response(JSON.stringify(obj), { status, headers: { 'Content-Type': 'application/json; charset=utf-8', ...extra } });

// Normaliza /api/asistente* (alias histórico) → /api/assistant* (canónico del core).
const normPath = (p) => p.replace(/^\/api\/asistente\b/, '/api/assistant');

// ── RAG: corpus de ejemplos + recuperación léxica (few-shot dinámico) ──────────
let _CORPUS = null;   // cache por instancia del Worker
async function cargarCorpus(env, base) {
  if (_CORPUS) return _CORPUS;
  try {
    const r = await env.ASSETS.fetch(new Request(new URL('/assistant/examples.json', base)));
    const data = await r.json();
    _CORPUS = Array.isArray(data.ejemplos) ? data.ejemplos : (Array.isArray(data) ? data : []);
  } catch { _CORPUS = []; }
  return _CORPUS;
}
const STOP_RAG = new Set(['de', 'del', 'la', 'el', 'los', 'las', 'y', 'o', 'en', 'con', 'un', 'una', 'para', 'por', 'que', 'es', 'al', 'm', 'cm', 'mm', 'cada', 'tipo', 'the', 'a', 'an', 'of', 'with', 'and', 'in', 'on', 'to', 'for']);
const tokRAG = (s) => String(s || '').toLowerCase().normalize('NFD').replace(/[^a-z0-9 ]/g, ' ').split(/\s+/).filter((t) => t.length > 1 && !STOP_RAG.has(t));
// El contenido del par few-shot puede venir bajo 'spec' (nuevo) o 'ficha' (compat).
const ejSpec = (ej) => ej.spec ?? ej.ficha ?? {};
// Recupera los k ejemplos más parecidos (solape de tokens) al mensaje del usuario.
function recuperarEjemplos(mensaje, corpus, k = 3) {
  const q = new Set(tokRAG(mensaje));
  if (!q.size || !corpus.length) return { ejemplos: [], score: 0 };
  const rank = corpus.map((ej) => {
    const t = tokRAG(`${ej.desc} ${JSON.stringify(ejSpec(ej).typology || ejSpec(ej).tipologia || '')}`);
    let s = 0; for (const w of t) if (q.has(w)) s++;
    return { ej, score: s / Math.sqrt(t.length || 1) };
  }).sort((a, b) => b.score - a.score);
  const top = rank.filter((r) => r.score > 0).slice(0, k);
  return { ejemplos: top.map((r) => r.ej), score: top.length ? top[0].score : 0 };
}

// Proveedores del LLM EN ORDEN DE PRIORIDAD. Se intentan en cascada: primero OpenAI
// (si hay OPENAI_API_KEY), y si TODOS sus modelos fallan, se pasa a OpenRouter
// (modelos gratis). Ambos usan el formato chat/completions de OpenAI.
function proveedoresLLM(env) {
  const lista = [];
  if (env.OPENAI_API_KEY) lista.push({
    nombre: 'OpenAI',
    url: 'https://api.openai.com/v1/chat/completions',
    key: env.OPENAI_API_KEY,
    modelos: [env.OPENAI_MODEL || 'gpt-4o-mini'],
    extraHeaders: {},
  });
  if (env.OPENROUTER_API_KEY) lista.push({
    nombre: 'OpenRouter',
    url: 'https://openrouter.ai/api/v1/chat/completions',
    key: env.OPENROUTER_API_KEY,
    modelos: env.OPENROUTER_MODEL ? [env.OPENROUTER_MODEL] : MODELS_FREE,
    extraHeaders: { 'X-Title': 'PORTICO-UACh Asistente' },
  });
  return lista;
}

// Llamada genérica chat/completions con un system prompt arbitrario (JSON mode)
// + mensajes few-shot opcionales (RAG).
async function chatJSON(prov, modelo, system, userContent, fewshotMsgs = []) {
  return fetch(prov.url, {
    method: 'POST',
    headers: { Authorization: `Bearer ${prov.key}`, 'Content-Type': 'application/json', ...prov.extraHeaders },
    body: JSON.stringify({
      model: modelo,
      temperature: 0,
      response_format: { type: 'json_object' },
      messages: [{ role: 'system', content: system }, ...fewshotMsgs, { role: 'user', content: userContent }],
    }),
  });
}

// Extrae el primer objeto JSON del contenido devuelto por el LLM.
function parseJSONObj(data) {
  let raw = String(data.choices?.[0]?.message?.content ?? '').trim()
    .replace(/^```(json)?/i, '').replace(/```$/, '').trim();
  const i = raw.indexOf('{'), j = raw.lastIndexOf('}');
  if (i < 0 || j < 0) return null;
  try { return JSON.parse(raw.slice(i, j + 1)); } catch { return null; }
}

// Mensaje del usuario → spec (cascada de proveedores/modelos, con few-shot RAG).
async function specDesdeLLM(mensaje, env, fewshot = []) {
  const provs = proveedoresLLM(env);
  if (!provs.length) throw new Error('Falta el secreto OPENAI_API_KEY u OPENROUTER_API_KEY en el Worker.');
  const fewMsgs = [];
  for (const ej of fewshot) fewMsgs.push({ role: 'user', content: ej.desc }, { role: 'assistant', content: JSON.stringify(ejSpec(ej)) });
  const intentos = [];
  for (const prov of provs) {
    for (const modelo of prov.modelos) {
      const r = await chatJSON(prov, modelo, SYSTEM, mensaje, fewMsgs);
      if (r.ok) {
        const data = await r.json();
        if (data.error) { intentos.push(`${prov.nombre}/${modelo}: ${data.error.message || JSON.stringify(data.error)}`); continue; }
        const spec = parseJSONObj(data);
        if (!spec) { intentos.push(`${prov.nombre}/${modelo}: no devolvió JSON`); continue; }
        return { spec, llm: { proveedor: prov.nombre, modelo: data.model || modelo, intentos } };
      }
      intentos.push(`${prov.nombre}/${modelo}: HTTP ${r.status} ${(await r.text()).slice(0, 160)}`);
      if (r.status === 401 || r.status === 403) break;   // credencial/política: no insistir con este proveedor
    }
  }
  throw new Error(`Ningún modelo disponible. Intentos: ${intentos.join(' | ')}`);
}

// Orden de modificación → { ops:[...] } (cascada, como el spec).
async function opsDesdeLLM(payload, env) {
  const provs = proveedoresLLM(env);
  if (!provs.length) throw new Error('Falta el secreto OPENAI_API_KEY u OPENROUTER_API_KEY en el Worker.');
  const userContent = JSON.stringify(payload);   // { message, model, selection }
  const intentos = [];
  for (const prov of provs) {
    for (const modelo of prov.modelos) {
      const r = await chatJSON(prov, modelo, MOD_SYSTEM, userContent);
      if (r.ok) {
        const data = await r.json();
        if (data.error) { intentos.push(`${prov.nombre}/${modelo}: ${data.error.message || JSON.stringify(data.error)}`); continue; }
        const parsed = parseJSONObj(data);
        if (!parsed) { intentos.push(`${prov.nombre}/${modelo}: no devolvió JSON`); continue; }
        const ops = Array.isArray(parsed) ? parsed : (parsed.ops || parsed.operaciones || []);
        return { ops, llm: { proveedor: prov.nombre, modelo: data.model || modelo, intentos } };
      }
      intentos.push(`${prov.nombre}/${modelo}: HTTP ${r.status} ${(await r.text()).slice(0, 160)}`);
      if (r.status === 401 || r.status === 403) break;
    }
  }
  throw new Error(`Ningún modelo disponible. Intentos: ${intentos.join(' | ')}`);
}

// ── Registro de consultas en KV (revisión semanal) ────────────────────────────
// estado: 'ok' (generó), 'error' (falló LLM), 'incorrecto' (feedback del usuario:
// no era lo solicitado). 'novedoso' es ortogonal (score RAG bajo).
// Devuelve la clave del registro (para que la app pueda enviar feedback luego).
async function registrarConsulta(env, { mensaje, spec = null, rag = null, llm = null, estado = 'ok', error = null }) {
  if (!env.ASIS_LOG) return null;
  const ts = Date.now();
  const key = `q:${ts}-${Math.random().toString(36).slice(2, 6)}`;
  const registro = {
    id: key, ts, fecha: new Date(ts).toISOString(), estado,
    mensaje, spec: spec || null, tipologia: spec?.typology || spec?.tipologia || null,
    score: rag?.score ?? null, novedoso: !!rag?.novedoso,
    modelo: llm?.modelo || null, error: error || null, comentario: null,
  };
  try {
    await env.ASIS_LOG.put(key, JSON.stringify(registro), {
      expirationTtl: 60 * 60 * 24 * 180,
      metadata: { estado, novedoso: registro.novedoso, tipologia: registro.tipologia, score: registro.score },
    });
    return key;
  } catch { return null; }
}

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const path = normPath(url.pathname);   // acepta /api/asistente* como alias

    // ── Lectura del registro de consultas (revisión semanal) ──
    // GET /api/assistant/log?token=TOKEN[&solo_novedosos=1][&estado=ok|error|incorrecto][&limite=50]
    if (path === '/api/assistant/log') {
      if (!env.ASIS_LOG) return json({ error: 'KV ASIS_LOG no está configurado en el Worker.' }, 400);
      const token = url.searchParams.get('token') || request.headers.get('x-asis-token');
      if (!env.ASIS_LOG_TOKEN || token !== env.ASIS_LOG_TOKEN) return json({ error: 'Token inválido (defina el secreto ASIS_LOG_TOKEN y páselo en ?token=).' }, 401);
      const soloNov = url.searchParams.get('solo_novedosos') === '1';
      const estadoF = url.searchParams.get('estado');
      const limite = Math.min(500, Math.max(1, parseInt(url.searchParams.get('limite') || '100', 10)));
      const lista = await env.ASIS_LOG.list({ prefix: 'q:', limit: 1000 });
      const conteo = { ok: 0, error: 0, incorrecto: 0, novedoso: 0 };
      for (const k of lista.keys) {
        const e = (k.metadata && k.metadata.estado) || 'ok';
        if (conteo[e] != null) conteo[e]++;
        if (k.metadata && k.metadata.novedoso) conteo.novedoso++;
      }
      let keys = lista.keys.sort((a, b) => b.name.localeCompare(a.name));
      if (soloNov) keys = keys.filter((k) => k.metadata && k.metadata.novedoso);
      if (estadoF) keys = keys.filter((k) => ((k.metadata && k.metadata.estado) || 'ok') === estadoF);
      keys = keys.slice(0, limite);
      const items = await Promise.all(keys.map(async (k) => { try { return JSON.parse(await env.ASIS_LOG.get(k.name)); } catch { return null; } }));
      const reg = items.filter(Boolean);
      // Candidatos a corpus: SOLO novedosos que generaron bien. Formato {desc, spec}
      // listo para revisar y pegar en assistant/examples.json.
      const corpus_sugerido = reg.filter((r) => r.novedoso && r.estado === 'ok').map((r) => ({ desc: r.mensaje, spec: r.spec }));
      const revisar = reg.filter((r) => r.estado === 'error' || r.estado === 'incorrecto')
        .map((r) => ({ id: r.id, estado: r.estado, mensaje: r.mensaje, error: r.error || null, comentario: r.comentario || null, tipologia: r.tipologia, spec: r.spec }));
      return json({ total: reg.length, conteo, registros: reg, corpus_sugerido, revisar });
    }

    // ── Generar spec desde texto (o devolver el spec recibido) ──
    if (path === '/api/assistant') {
      if (request.method !== 'POST') return json({ error: 'Use POST' }, 405);
      let body;
      try { body = await request.json(); } catch { return json({ error: 'JSON inválido en la solicitud' }, 400); }
      const message = body.message ?? body.mensaje ?? null;          // inglés (core) o español (compat)
      let spec = body.spec ?? body.ficha ?? null;
      const desdeMensaje = !spec && !!message;   // solo registramos lo que vino del LLM
      let llm = null, rag = null;
      try {
        if (desdeMensaje) {
          // RAG: recuperar ejemplos parecidos del corpus e inyectarlos como few-shot
          const corpus = await cargarCorpus(env, request.url);
          const rec = recuperarEjemplos(message, corpus, 3);
          rag = { usados: rec.ejemplos.map((e) => String(e.desc).slice(0, 60)), score: +rec.score.toFixed(3), novedoso: rec.score < 0.5 };
          const res = await specDesdeLLM(message, env, rec.ejemplos);
          spec = res.spec; llm = res.llm;
        }
        if (!spec) return json({ error: 'Envíe { message } o { spec }' }, 400);
        // Registro OK (candidato a corpus). Devuelve la clave para feedback posterior.
        const logId = desdeMensaje ? await registrarConsulta(env, { mensaje: message, spec, rag, llm, estado: 'ok' }) : null;
        const hdr = llm ? { 'X-Asistente-Proveedor': llm.proveedor, 'X-Asistente-Modelo': String(llm.modelo) } : {};
        // El CLIENTE construye el modelo desde spec (generator.js del core).
        return json({ spec, _llm: llm, _rag: rag, _logId: logId }, 200, hdr);
      } catch (e) {
        const msg = String(e.message || e);
        if (desdeMensaje) { try { await registrarConsulta(env, { mensaje: message, spec, rag, llm, estado: 'error', error: msg }); } catch { /* no bloquear */ } }
        return json({ error: msg }, 500);
      }
    }

    // ── MODIFICAR el modelo ya construido: orden NL → operaciones ──
    // POST /api/assistant/modificar  body { message, model?, selection? } → { ops }
    if (path === '/api/assistant/modificar') {
      if (request.method !== 'POST') return json({ error: 'Use POST' }, 405);
      let body;
      try { body = await request.json(); } catch { return json({ error: 'JSON inválido en la solicitud' }, 400); }
      const message = body.message ?? body.mensaje ?? null;
      if (!message) return json({ error: 'Envíe { message }' }, 400);
      const model = body.model ?? body.modelo ?? null;
      const selection = body.selection ?? body.seleccion ?? null;
      try {
        const { ops, llm } = await opsDesdeLLM({ message, model, selection }, env);
        const hdr = llm ? { 'X-Asistente-Proveedor': llm.proveedor, 'X-Asistente-Modelo': String(llm.modelo) } : {};
        return json({ ops, _llm: llm }, 200, hdr);
      } catch (e) {
        return json({ error: String(e.message || e) }, 500);
      }
    }

    // ── Feedback del usuario: marcar una consulta como 'incorrecto' ──
    // POST /api/assistant/feedback  body { id, comentario? }
    if (path === '/api/assistant/feedback') {
      if (request.method !== 'POST') return json({ error: 'Use POST' }, 405);
      if (!env.ASIS_LOG) return json({ error: 'Registro KV no configurado' }, 400);
      let body;
      try { body = await request.json(); } catch { return json({ error: 'JSON inválido' }, 400); }
      const id = String(body.id || '');
      if (!id.startsWith('q:')) return json({ error: 'id inválido' }, 400);
      const raw = await env.ASIS_LOG.get(id);
      if (!raw) return json({ error: 'registro no encontrado' }, 404);
      let reg; try { reg = JSON.parse(raw); } catch { return json({ error: 'registro corrupto' }, 500); }
      reg.estado = 'incorrecto';
      reg.comentario = body.comentario ? String(body.comentario).slice(0, 500) : null;
      reg.feedback_ts = Date.now();
      try {
        await env.ASIS_LOG.put(id, JSON.stringify(reg), {
          expirationTtl: 60 * 60 * 24 * 180,
          metadata: { estado: 'incorrecto', novedoso: !!reg.novedoso, tipologia: reg.tipologia, score: reg.score },
        });
      } catch (e) { return json({ error: String(e.message || e) }, 500); }
      return json({ ok: true, id, estado: reg.estado });
    }

    // Resto: servir la PWA (assets estáticos: overlay + vendor/portico-core)
    return env.ASSETS.fetch(request);
  },
};
