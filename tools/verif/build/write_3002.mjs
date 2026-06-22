import { build3002 } from './build_3002.mjs';
import { Serializer } from '../../../js/model/serializer.js';
import fs from 'fs';
const { m } = build3002(6, 1);
fs.writeFileSync('examples/verif_3-002_plane_stress.s3d', new Serializer().toJSON(m));
console.log('escrito 6x1; nodos', m.nodes.size, 'areas', m.areas.size);
