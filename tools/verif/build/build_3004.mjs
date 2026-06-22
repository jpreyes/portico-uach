// Construye el .s3d del caso 3-004 (cilindro de pared gruesa, plane-strain/stress, #58).
// Timoshenko 1956 / MacNeal-Harder 1985 / CSI Example 3-004. Cuarto de cilindro
// (simetría alineada a ejes: borde θ=0 restringe uz, borde θ=90° restringe ux),
// presión interna P=1 ksi. Malla radial de 5 bandas (los radios del original).
import { Model } from '../../../js/model/model.js';
import { Serializer } from '../../../js/model/serializer.js';
import fs from 'fs';

// radii del original (5 bandas), nθ segmentos circunferenciales en 90°.
export function build3004(planeStrain, nu = 0.3, nTheta = 9) {
  const R = [3, 3.5, 4.2, 5.2, 6.75, 9], nR = R.length - 1;
  const P = 1, t = 1, E = 1000;
  const m = new Model();
  m.mode = '3D'; m.units = 'kip-in';
  m.materials.clear(); m.sections.clear();
  const mat = m.addMaterial({ name: 'Mat', E, G: E / (2 * (1 + nu)), nu, alpha: 0, rho: 0 });

  // Nodos: anillo i (radio R[i]) × ángulo j (θ_j en [0,90°]). Sección en plano X-Z.
  const idAt = new Map();
  const key = (i, j) => `${i},${j}`;
  for (let i = 0; i <= nR; i++) for (let j = 0; j <= nTheta; j++) {
    const th = (Math.PI / 2) * j / nTheta;
    const x = R[i] * Math.cos(th), z = R[i] * Math.sin(th);
    const node = m.addNode(x, 0, z, { uy: 1, rx: 1, ry: 1, rz: 1 });   // membrana X-Z
    idAt.set(key(i, j), node.id);
  }
  // Simetría: borde θ=0 (z=0) restringe uz; borde θ=90° (x=0) restringe ux.
  for (let i = 0; i <= nR; i++) {
    m.updateNode(idAt.get(key(i, 0)),      { restraints: { uz: 1 } });
    m.updateNode(idAt.get(key(i, nTheta)), { restraints: { ux: 1 } });
  }
  // Áreas QUAD
  for (let i = 0; i < nR; i++) for (let j = 0; j < nTheta; j++)
    m.addArea([idAt.get(key(i, j)), idAt.get(key(i + 1, j)), idAt.get(key(i + 1, j + 1)), idAt.get(key(i, j + 1))],
      mat.id, { thickness: t, behavior: 'membrane', planeStrain });

  // Presión interna P en la cara interna (anillo i=0): fuerza radial saliente por nodo
  // = P·t·(arco tributario). Dirección radial (cosθ,0,sinθ).
  const lc = m.addLoadCase('Presion', false);
  const r1 = R[0], dth = (Math.PI / 2) / nTheta;
  for (let j = 0; j <= nTheta; j++) {
    const trib = (j === 0 || j === nTheta) ? dth / 2 : dth;   // medio en los extremos
    const F = P * t * r1 * trib;
    const th = dth * j;
    m.addLoad(lc.id, { type: 'nodal', nodeId: idAt.get(key(0, j)), F: [F * Math.cos(th), 0, F * Math.sin(th), 0, 0, 0] });
  }
  const innerTheta0 = idAt.get(key(0, 0));   // nodo cara interna en θ=0 → radial = ux
  return { m, lc, innerTheta0 };
}

if (process.argv[2] === 'write') {
  const { m } = build3004(true, 0.3, 9);
  fs.writeFileSync('examples/verif_3-004_plane_strain_cilindro.s3d', new Serializer().toJSON(m));
  console.log('escrito plane-strain ν=0.3; nodos', m.nodes.size, 'areas', m.areas.size);
}
