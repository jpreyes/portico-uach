// ── scripts/create-icons.js ────────────────────────────────────────────────
// Generates assets/icon.png (1024×1024), assets/icon.ico (Windows),
// and assets/icon.icns (macOS) from an SVG source.
//
// Run once:  node scripts/create-icons.js
// Then commit the generated files — CI uses them directly.
// ──────────────────────────────────────────────────────────────────────────
const sharp     = require('sharp');
const png2icons = require('png2icons');
const fs        = require('fs');
const path      = require('path');

// ── Icon design: 3D portal frame with mode shape ──────────────────────────
// Dark-blue background, blue frame (3D perspective), yellow deflected curve,
// red pinned nodes, red fixed-support symbols at base.
const SVG = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
  <!-- background -->
  <rect width="512" height="512" rx="88" fill="#0D47A1"/>

  <!-- back frame (perspective offset) -->
  <g fill="none" stroke="#42A5F5" stroke-width="5" opacity="0.55">
    <line x1="168" y1="80"  x2="358" y2="80"/>
    <line x1="168" y1="80"  x2="168" y2="240"/>
    <line x1="358" y1="80"  x2="358" y2="240"/>
  </g>

  <!-- perspective connectors -->
  <line x1="96"  y1="170" x2="168" y2="80"  stroke="#64B5F6" stroke-width="5"/>
  <line x1="416" y1="170" x2="358" y2="80"  stroke="#64B5F6" stroke-width="5"/>

  <!-- front frame -->
  <g fill="none" stroke="#90CAF9" stroke-width="10" stroke-linecap="round" stroke-linejoin="round">
    <line x1="96"  y1="170" x2="96"  y2="390"/>
    <line x1="416" y1="170" x2="416" y2="390"/>
    <line x1="96"  y1="170" x2="416" y2="170"/>
  </g>

  <!-- mode shape (1st mode parabola) -->
  <path d="M 96 170 Q 256 95 416 170"
        stroke="#FFD54F" stroke-width="11"
        fill="none" stroke-linecap="round"/>

  <!-- nodes at frame corners -->
  <circle cx="96"  cy="170" r="16" fill="#EF5350"/>
  <circle cx="416" cy="170" r="16" fill="#EF5350"/>
  <!-- anti-node at peak -->
  <circle cx="256" cy="116" r="11" fill="#FFFFFF"/>

  <!-- fixed supports at column bases -->
  <!-- left support -->
  <rect  x="68"  y="391" width="56" height="9"  fill="#EF5350" rx="2"/>
  <line x1="72"  y1="400" x2="63"  y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="84"  y1="400" x2="75"  y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="96"  y1="400" x2="87"  y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="108" y1="400" x2="99"  y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="120" y1="400" x2="111" y2="413" stroke="#EF5350" stroke-width="4"/>
  <!-- right support -->
  <rect  x="388" y="391" width="56" height="9"  fill="#EF5350" rx="2"/>
  <line x1="392" y1="400" x2="383" y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="404" y1="400" x2="395" y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="416" y1="400" x2="407" y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="428" y1="400" x2="419" y2="413" stroke="#EF5350" stroke-width="4"/>
  <line x1="440" y1="400" x2="431" y2="413" stroke="#EF5350" stroke-width="4"/>
</svg>`;

// ── Main ──────────────────────────────────────────────────────────────────
async function main() {
  const assetsDir = path.join(__dirname, '..', 'assets');
  fs.mkdirSync(assetsDir, { recursive: true });

  // 1. SVG → 1024×1024 PNG
  process.stdout.write('Generating PNG...  ');
  const pngBuffer = await sharp(Buffer.from(SVG))
    .resize(1024, 1024)
    .png()
    .toBuffer();
  fs.writeFileSync(path.join(assetsDir, 'icon.png'), pngBuffer);
  console.log('✓  assets/icon.png');

  // 2. PNG → ICO (multi-size: 16, 24, 32, 48, 64, 128, 256)
  process.stdout.write('Generating ICO...  ');
  const ico = png2icons.createICO(pngBuffer, png2icons.BILINEAR, 0, true, false);
  if (!ico) throw new Error('ICO creation failed');
  fs.writeFileSync(path.join(assetsDir, 'icon.ico'), ico);
  console.log('✓  assets/icon.ico');

  // 3. PNG → ICNS (macOS)
  process.stdout.write('Generating ICNS... ');
  const icns = png2icons.createICNS(pngBuffer, png2icons.BILINEAR, 0);
  if (!icns) throw new Error('ICNS creation failed');
  fs.writeFileSync(path.join(assetsDir, 'icon.icns'), icns);
  console.log('✓  assets/icon.icns');

  console.log('\nDone. Commit assets/ and these will be used by electron-builder.');
}

main().catch(err => { console.error(err); process.exit(1); });
