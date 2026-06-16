// ──────────────────────────────────────────────────────────────────────────────
// Cloudflare Worker — sirve la PWA (assets) y expone la API del asistente.
//   POST /api/asistente   body { mensaje }  → LLM (OpenRouter) → ficha → modelo
//                         body { ficha }    → genera directo (sin LLM)
//   resto de rutas        → assets estáticos (la app PÓRTICO)
//
// La API key de OpenRouter vive como SECRETO del Worker (env.OPENROUTER_API_KEY),
// nunca en el código ni en el navegador:
//   npx wrangler secret put OPENROUTER_API_KEY
// ──────────────────────────────────────────────────────────────────────────────
import { generarModelo } from '../asistente/generador.js';

const MODEL_DEFAULT = 'meta-llama/llama-3.3-70b-instruct:free';

const SYSTEM = `Eres un asistente que convierte la descripcion de una estructura en una FICHA JSON para PORTICO. Responde SOLO con el JSON de la ficha, sin texto ni markdown. Campos: proyecto, modo (2D|3D), ubicacion{ciudad,latitud_sur_deg,altitud_msnm,exposicion(B|C|D)}, geometria{niveles:[{altura_m,uso_NCh1537?,sobrecarga_uso_kN_m2?}],vanos_x?,vanos_y?,planta_inferior?{Lx_m,Ly_m},planta_superior?,pendiente_techo_deg?}, secciones{material,vigas,pilares}, apoyo_base(empotrado|rotulado), diafragma_rigido, cargas{muerta_adicional_kN_m2,uso_NCh1537,cierre_viento,nieve,viento,sismo}, sismo{zona(1|2|3),suelo(A..E),categoria(I..IV),R}. vanos_x/vanos_y: lista de luces [3,3,3,4] o uniforme {cantidad,luz_m} (ej. "4 vanos de 3 m en X" -> vanos_x:{cantidad:4,luz_m:3}). Cada nivel puede tener distinta altura y distinto uso (ej. nivel 1 Salas de Clases, nivel 3 Bodegas livianas) -> ponlo en niveles[k].uso_NCh1537. Omite lo no mencionado; no inventes valores de ingenieria. Materiales tipicos: S275, A630-420H; perfiles IPE300, HEB200.`;

function parseCSV(txt) {
  const lines = txt.split(/\r?\n/).filter((l) => l.trim() && !l.startsWith('#'));
  const head = lines[0].split(',').map((s) => s.trim());
  return lines.slice(1).map((l) => {
    const c = l.split(',').map((s) => s.trim());
    return Object.fromEntries(head.map((h, i) => [h, c[i]]));
  });
}

const json = (obj, status = 200) =>
  new Response(JSON.stringify(obj), { status, headers: { 'Content-Type': 'application/json; charset=utf-8' } });

async function cargarBibliotecas(env, base) {
  const get = (p) => env.ASSETS.fetch(new Request(new URL(p, base)));
  const [reglas, pTxt, mTxt, sTxt] = await Promise.all([
    get('/asistente/reglas.json').then((r) => r.json()),
    get('/asistente/perfiles.csv').then((r) => r.text()),
    get('/asistente/materiales.csv').then((r) => r.text()),
    get('/asistente/sobrecargas_NCh1537.csv').then((r) => r.text()),
  ]);
  return { reglas, perfiles: parseCSV(pTxt), materiales: parseCSV(mTxt), sobrecargas: parseCSV(sTxt) };
}

async function fichaDesdeLLM(mensaje, env) {
  if (!env.OPENROUTER_API_KEY) throw new Error('Falta el secreto OPENROUTER_API_KEY en el Worker.');
  const r = await fetch('https://openrouter.ai/api/v1/chat/completions', {
    method: 'POST',
    headers: {
      Authorization: `Bearer ${env.OPENROUTER_API_KEY}`,
      'Content-Type': 'application/json',
      'X-Title': 'PORTICO Asistente',
    },
    body: JSON.stringify({
      model: env.OPENROUTER_MODEL || MODEL_DEFAULT,
      temperature: 0,
      response_format: { type: 'json_object' },
      messages: [{ role: 'system', content: SYSTEM }, { role: 'user', content: mensaje }],
    }),
  });
  if (!r.ok) throw new Error(`OpenRouter HTTP ${r.status}: ${(await r.text()).slice(0, 300)}`);
  const data = await r.json();
  let raw = String(data.choices?.[0]?.message?.content ?? '').trim()
    .replace(/^```(json)?/i, '').replace(/```$/, '').trim();
  const i = raw.indexOf('{'), j = raw.lastIndexOf('}');
  if (i < 0 || j < 0) throw new Error('El LLM no devolvió JSON: ' + raw.slice(0, 200));
  return JSON.parse(raw.slice(i, j + 1));
}

export default {
  async fetch(request, env) {
    const url = new URL(request.url);

    if (url.pathname === '/api/asistente') {
      if (request.method !== 'POST') return json({ error: 'Use POST' }, 405);
      try {
        const body = await request.json();
        const ficha = body.ficha ?? (body.mensaje ? await fichaDesdeLLM(body.mensaje, env) : null);
        if (!ficha) return json({ error: 'Envíe { mensaje } o { ficha }' }, 400);
        const libs = await cargarBibliotecas(env, request.url);
        const modelo = generarModelo(ficha, libs);
        return json({ ficha, resumen: modelo._generado?.resumen, modelo });
      } catch (e) {
        return json({ error: String(e.message || e) }, 500);
      }
    }

    // Resto: servir la PWA (assets estáticos)
    return env.ASSETS.fetch(request);
  },
};
