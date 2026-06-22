// Construye el .s3d del caso 3-002 (viga voladizo modelada con plane-stress, #58).
// MacNeal-Harder 1985 / CSI Example 3-002. Malla nx×nz de QUADs membrana.
//   3 casos de carga: (1) extensión axial, (2) corte+flexión en plano, (3) momento.
import { Model } from '../../../js/model/model.js';
import { Serializer } from '../../../js/model/serializer.js';
import fs from 'fs';

export function build3002(nx = 6, nz = 1) {
  const L = 6, H = 0.2, t = 0.1;
  const m = new Model();
  m.mode = '3D'; m.units = 'lb-in';
  m.materials.clear(); m.sections.clear();
  const mat = m.addMaterial({ name: 'Plane', E: 10e6, G: 3.846154e6, nu: 0.3, alpha: 0, rho: 0 });

  // Nodos: rejilla (nx+1)×(nz+1) en el plano X-Z. Sólo Ux,Uz activos (resto fijo).
  const id = (i, k) => k * (nx + 1) + i + 1;
  for (let k = 0; k <= nz; k++) for (let i = 0; i <= nx; i++) {
    const x = L * i / nx, z = H * k / nz;
    m.addNode(x, 0, z, { uy: 1, rx: 1, ry: 1, rz: 1 });   // membrana X-Z: sólo ux,uz
  }
  // Empotramiento (x=0): nodo inferior izq Ux,Uz; nodo superior izq sólo Ux (sin Poisson)
  m.updateNode(id(0, 0), { restraints: { ux: 1, uz: 1 } });
  for (let k = 1; k <= nz; k++) m.updateNode(id(0, k), { restraints: { ux: 1 } });

  // Áreas QUAD membrana plane-stress
  for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++)
    m.addArea([id(i, k), id(i + 1, k), id(i + 1, k + 1), id(i, k + 1)],
      mat.id, { thickness: t, behavior: 'membrane', planeStrain: false });

  // Nodos de la punta (x=L) por fila
  const tip = []; for (let k = 0; k <= nz; k++) tip.push(id(nx, k));
  const nTip = tip.length;
  const topLeft = id(0, nz);   // jt 8 del original

  // LC1 — extensión axial: Fx total +1 repartido en la punta
  const lc1 = m.addLoadCase('Axial', false);
  for (const j of tip) m.addLoad(lc1.id, { type: 'nodal', nodeId: j, F: [1 / nTip, 0, 0, 0, 0, 0] });
  // LC2 — corte+flexión: Fz total +1 en la punta; reacción −1 en el nodo sup-izq (sin Poisson)
  const lc2 = m.addLoadCase('Corte', false);
  for (const j of tip) m.addLoad(lc2.id, { type: 'nodal', nodeId: j, F: [0, 0, 1 / nTip, 0, 0, 0] });
  m.addLoad(lc2.id, { type: 'nodal', nodeId: topLeft, F: [0, 0, -1, 0, 0, 0] });
  // LC3 — momento M=1 en plano: par de fuerzas Fx en la punta (±M/H en sup/inf)
  const lc3 = m.addLoadCase('Momento', false);
  const Fcouple = 1 / H;   // M=1 → fuerzas ±1/H separadas H
  m.addLoad(lc3.id, { type: 'nodal', nodeId: id(nx, 0), F: [-Fcouple, 0, 0, 0, 0, 0] });
  m.addLoad(lc3.id, { type: 'nodal', nodeId: id(nx, nz), F: [+Fcouple, 0, 0, 0, 0, 0] });

  return { m, tip, lc1, lc2, lc3 };
}

if (import.meta.url === `file://${process.argv[1].replace(/\\/g, '/')}`) {
  const { m } = build3002(6, 1);
  fs.writeFileSync('examples/verif_3-002_plane_stress.s3d', new Serializer().toJSON(m));
  console.log('escrito 6x1; nodos', m.nodes.size, 'areas', m.areas.size);
}
