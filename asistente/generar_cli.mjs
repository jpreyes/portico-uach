#!/usr/bin/env node
// ──────────────────────────────────────────────────────────────────────────────
// generar_cli.mjs — wrapper de línea de comandos para el generador determinista.
// Lee una FICHA (JSON) por stdin y emite el modelo .s3d (JSON) por stdout.
// Pensado para el nodo "Execute Command" de n8n:  node asistente/generar_cli.mjs
// Errores → stderr + exit 1 (n8n los captura).
//
// Uso:
//   echo '{"modo":"3D",...}' | node asistente/generar_cli.mjs
//   node asistente/generar_cli.mjs ruta/ficha.json
// ──────────────────────────────────────────────────────────────────────────────
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { generarModelo } from './generador.js';

const dir = path.dirname(fileURLToPath(import.meta.url));
const read = (p) => fs.readFileSync(path.join(dir, p), 'utf8');

function parseCSV(txt) {
  const lines = txt.split(/\r?\n/).filter((l) => l.trim() && !l.startsWith('#'));
  const head = lines[0].split(',').map((s) => s.trim());
  return lines.slice(1).map((l) => {
    const c = l.split(',').map((s) => s.trim());
    return Object.fromEntries(head.map((h, i) => [h, c[i]]));
  });
}

async function leerEntrada() {
  const arg = process.argv[2];
  if (arg) return fs.readFileSync(arg, 'utf8');
  // stdin
  const chunks = [];
  for await (const c of process.stdin) chunks.push(c);
  return Buffer.concat(chunks).toString('utf8');
}

try {
  const fichaText = (await leerEntrada()).trim();
  if (!fichaText) throw new Error('Ficha vacía: pásela por stdin o como argumento.');
  const ficha = JSON.parse(fichaText);

  const libs = {
    reglas: JSON.parse(read('reglas.json')),
    perfiles: parseCSV(read('perfiles.csv')),
    materiales: parseCSV(read('materiales.csv')),
    sobrecargas: parseCSV(read('sobrecargas_NCh1537.csv')),
  };

  const modelo = generarModelo(ficha, libs);
  process.stdout.write(JSON.stringify(modelo, null, 2));
} catch (e) {
  process.stderr.write(`generar_cli: ${e.message}\n`);
  process.exit(1);
}
