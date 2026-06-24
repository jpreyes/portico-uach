// ──────────────────────────────────────────────────────────────────────────────
// io/index.js — punto de entrada del módulo de INTEROPERABILIDAD (#74, G18)
//
// Importa los adaptadores (que se auto-registran por efecto colateral) y re-exporta la
// API del registro.  Para agregar un motor nuevo: crear `formats/<motor>.js` que llame
// a `registerFormat({ id, name, ext, write, read })` y añadirlo a la lista de imports.
// ──────────────────────────────────────────────────────────────────────────────
export { registerFormat, getFormat, listFormats, exportModel, importModel } from './registry.js?v=188';
export { modelToNeutral, neutralToModel } from './neutral.js?v=188';

// Adaptadores de formato (auto-registro):
import './formats/vector.js?v=188';
import './formats/abaqus.js?v=188';
import './formats/sap2000.js?v=188';
// Próximos: ETABS (.e2k), OpenSees (.tcl), SOFISTIK (.dat)
