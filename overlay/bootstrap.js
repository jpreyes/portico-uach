// ──────────────────────────────────────────────────────────────────────────────
// bootstrap.js — capa overlay UACh sobre portico-core.
//
// NO forkea app.js: importa los seams de `js/ext/extensions.js` del core (mismo
// singleton, vía la misma URL del submodule) y REGISTRA las contribuciones UACh.
// Se carga como módulo ANTES de app.js en index.html; como los módulos son
// diferidos y se ejecutan en orden, el registro ocurre antes de que App se cree
// en DOMContentLoaded y las pinte (App._initExtensions → #ext-badges).
// ──────────────────────────────────────────────────────────────────────────────
import { extensions } from '../vendor/portico-core/js/ext/extensions.js?v=2';
import { EN } from '../vendor/portico-core/js/i18n/dict.en.js?v=2';

// i18n del overlay: el motor traduce por «source string» (español = fuente). Los
// textos UACh (bloque legal de Ayuda → Acerca de) no están en el diccionario del
// core, así que inyectamos sus traducciones EN mutando el MISMO objeto EN (mismo
// módulo del submodule). Se hace antes de App._initI18n; _lookup lee EN en vivo, así
// que también surte efecto al cambiar de idioma. Clave = string ES tal como aparece
// en el DOM (cada frase es un único nodo de texto en index.html).
Object.assign(EN, {
  'Código fuente y licencia': 'Source code and license',
  'PÓRTICO-UACh es la capa académica de la Universidad Austral de Chile (Instituto de Obras Civiles), construida sobre el motor open source portico-core. Ambos se distribuyen bajo licencia GNU AGPL-3.0: usted puede usarlo, estudiarlo, modificarlo y redistribuirlo, y tiene derecho a obtener el código fuente correspondiente, incluso al usar la aplicación a través de la red (AGPL §13).':
    'PÓRTICO-UACh is the academic layer of the Universidad Austral de Chile (Instituto de Obras Civiles), built on the open source engine portico-core. Both are distributed under the GNU AGPL-3.0 license: you may use, study, modify and redistribute it, and you have the right to obtain the corresponding source code, even when using the application over a network (AGPL §13).',
  'Motor (portico-core):': 'Engine (portico-core):',
  'Capa académica (este sitio):': 'Academic layer (this site):',
  'Licencia:': 'License:',
  'Bibliotecas de terceros:': 'Third-party libraries:',
  '© Instituto de Obras Civiles, Universidad Austral de Chile. PÓRTICO-UACh es 100% académico: todas las funciones están disponibles, sin licencias de pago. El asistente IA, cuando está habilitado, usa un servicio LLM externo; el resto corre íntegramente en su navegador.':
    '© Instituto de Obras Civiles, Universidad Austral de Chile. PÓRTICO-UACh is 100% academic: all features are available, with no paid licenses. The AI assistant, when enabled, uses an external LLM service; everything else runs entirely in your browser.',
});

// Insignia académica (Universidad Austral de Chile). PÓRTICO-UACh es 100%
// académico: sin modo profesional ni tokens; todas las funciones disponibles.
extensions.registerBadge({
  id: 'uach-academic',
  html:
    '<span class="badge-academic" data-i18n-skip ' +
    'title="Versión académica — Instituto de Obras Civiles, Universidad Austral de Chile" ' +
    'style="display:inline-block;background:var(--teal,#0ea5a4);color:#fff;' +
    'font-size:11px;font-weight:600;letter-spacing:.04em;border-radius:6px;' +
    'padding:3px 8px;margin-right:8px;vertical-align:middle;white-space:nowrap">' +
    'ACADÉMICO</span>',
});

export default extensions;
