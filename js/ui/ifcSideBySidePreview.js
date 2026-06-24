// ──────────────────────────────────────────────────────────────────────────────
// ui/ifcSideBySidePreview.js — vista comparada IFC ↔ PÓRTICO · #77, G19
//
// Dos mini-visores Three.js lado a lado con CÁMARAS SINCRONIZADAS:
//   • izquierda  → lo que trae el IFC (ejes de todos los elementos con geometría,
//                  coloreados por estado: importable / sin geometría / no soportado).
//   • derecha    → lo que CREARÁ PÓRTICO (sólo los elementos seleccionados e importables,
//                  con sus NODOS tras el snap).  Se actualiza al marcar/desmarcar filas.
// Orbitar/zoom en un lado mueve el otro.  Convención de ejes igual que el viewport
// principal: model(x,y,z) → three(x, z, y).  Sin estado global; se destruye con dispose().
// ──────────────────────────────────────────────────────────────────────────────
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const COL = {
  bg: 0x0a0e1a,
  ok: 0x4caf50,        // importable
  nogeom: 0xef5350,    // sin geometría
  unsup: 0x607d8b,     // no soportado
  node: 0xffc107,      // nodo PÓRTICO
  portico: 0x4fc3f7,   // barra PÓRTICO
};
const m2t = (x, y, z) => new THREE.Vector3(x, z, y);

class SubView {
  constructor(host) {
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setClearColor(COL.bg, 1);
    this.renderer.domElement.style.cssText = 'width:100%;height:100%;display:block';
    host.appendChild(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.scene.add(new THREE.AmbientLight(0xffffff, 1));
    this.camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1e5);
    this.camera.position.set(10, 8, 12);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = false;
    this.group = new THREE.Group();
    this.scene.add(this.group);
    this.host = host;
  }
  clear() { this.group.clear(); }
  resize() {
    const w = this.host.clientWidth || 1, h = this.host.clientHeight || 1;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h; this.camera.updateProjectionMatrix();
  }
  render() { this.renderer.render(this.scene, this.camera); }
  dispose() { this.controls.dispose(); this.renderer.dispose(); this.renderer.forceContextLoss?.(); this.host.removeChild(this.renderer.domElement); }
}

export class IfcPreview {
  constructor(leftHost, rightHost) {
    this.left = new SubView(leftHost);
    this.right = new SubView(rightHost);
    this.items = [];
    this._target = new THREE.Vector3();
    this._syncing = false;

    // sincronización bidireccional de cámara
    const sync = (from, to) => {
      if (this._syncing) return; this._syncing = true;
      to.camera.position.copy(from.camera.position);
      to.controls.target.copy(from.controls.target);
      to.camera.updateProjectionMatrix(); to.controls.update();
      this._syncing = false;
    };
    this.left.controls.addEventListener('change', () => sync(this.left, this.right));
    this.right.controls.addEventListener('change', () => sync(this.right, this.left));

    this._raf = null;
    const loop = () => { this.left.render(); this.right.render(); this._raf = requestAnimationFrame(loop); };
    this._raf = requestAnimationFrame(loop);
    this._onResize = () => this._resize();
    window.addEventListener('resize', this._onResize);
  }

  _resize() { this.left.resize(); this.right.resize(); }

  // construye el lado IZQUIERDO (IFC) y encuadra ambas cámaras
  setData(items) {
    this.items = items;
    this.left.clear();
    const box = new THREE.Box3();
    const lines = { ok: [], nogeom: [], unsup: [] };
    for (const it of items) {
      const bucket = it.status === 'ok' ? lines.ok : (it.status === 'no-geom' ? lines.nogeom : lines.unsup);
      for (const [a, b] of (it.segments || [])) {
        const va = m2t(a[0], a[1], a[2]), vb = m2t(b[0], b[1], b[2]);
        bucket.push(va, vb); box.expandByPoint(va); box.expandByPoint(vb);
      }
    }
    for (const [k, pts] of Object.entries(lines)) if (pts.length) this.left.group.add(this._seg(pts, COL[k === 'ok' ? 'ok' : k === 'nogeom' ? 'nogeom' : 'unsup'], k === 'unsup' ? 0.35 : 1));

    // encuadre
    if (!box.isEmpty()) {
      box.getCenter(this._target);
      const r = Math.max(box.getSize(new THREE.Vector3()).length() * 0.6, 1);
      const cam = this._target.clone().add(new THREE.Vector3(r, r * 0.8, r));
      for (const v of [this.left, this.right]) { v.camera.position.copy(cam); v.controls.target.copy(this._target); v.controls.update(); }
    }
    this._resize();
  }

  // reconstruye el lado DERECHO (PÓRTICO) con la selección actual
  updateSelection(selected) {
    this.right.clear();
    const pts = [], nodeSet = new Map();
    for (const it of this.items) {
      if (it.status !== 'ok') continue;
      if (selected && !selected.has(it.ifcId)) continue;
      for (const [a, b] of (it.segments || [])) {
        const va = m2t(a[0], a[1], a[2]), vb = m2t(b[0], b[1], b[2]);
        pts.push(va, vb);
        nodeSet.set(a.join(','), va); nodeSet.set(b.join(','), vb);
      }
    }
    if (pts.length) this.right.group.add(this._seg(pts, COL.portico, 1));
    if (nodeSet.size) {
      const geo = new THREE.BufferGeometry().setFromPoints([...nodeSet.values()]);
      this.right.group.add(new THREE.Points(geo, new THREE.PointsMaterial({ color: COL.node, size: 6, sizeAttenuation: false })));
    }
  }

  _seg(pts, color, opacity) {
    const geo = new THREE.BufferGeometry().setFromPoints(pts);
    const mat = new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity });
    return new THREE.LineSegments(geo, mat);
  }

  dispose() {
    cancelAnimationFrame(this._raf);
    window.removeEventListener('resize', this._onResize);
    this.left.dispose(); this.right.dispose();
  }
}
