#!/usr/bin/env python3
import os
import argparse
import xml.etree.ElementTree as ET
import numpy as np
import amriso


def leaf_grids(xdmf_path):
    base = os.path.dirname(os.path.abspath(xdmf_path))
    root = ET.parse(xdmf_path).getroot()
    for grid in root.iter("Grid"):
        topo = grid.find("Topology")
        if topo is None or topo.get("Dimensions") is None:
            continue
        ncells = int(topo.get("Dimensions").split()[0])
        geom = os.path.join(base, grid.find("Geometry/DataItem").text.strip())
        tnode = grid.find("Time")
        time = float(tnode.get("Value")) if tnode is not None else 0.0
        attrs = {a.get("Name"): os.path.join(base, a.find("DataItem").text.strip())
                 for a in grid.findall("Attribute")}
        yield ncells, geom, attrs, time


def read_raw(path, count, per_item):
    nbytes = os.path.getsize(path)
    n = count * per_item
    if nbytes == n * 4:
        dt = np.float32
    elif nbytes == n * 8:
        dt = np.float64
    else:
        raise SystemExit("%s: %d bytes != %d * (4 or 8)" % (path, nbytes, n))
    return np.fromfile(path, dtype=dt)


FRAME = '''\
      <Grid>
        <Time
            Value="{time:.16e}"/>
        <Topology
            TopologyType="Polyline"
            NodesPerElement="2"
            Dimensions="{nseg}">
          <DataItem
              Dimensions="{nseg} 2"
              NumberType="Int"
              Format="Binary">
            {name}.seg.raw
          </DataItem>
        </Topology>
        <Geometry
            GeometryType="XY">
          <DataItem
              Dimensions="{nvert} 2"
              Precision="4"
              Format="Binary">
            {name}.xy.raw
          </DataItem>
        </Geometry>
        <Attribute
            Center="Node"
            Name="u">
          <DataItem
              Dimensions="{nvert}"
              Precision="4"
              Format="Binary">
            {name}.attr.raw
          </DataItem>
        </Attribute>
      </Grid>
'''


def dump_raws(prefix, xy, seg, attr):
    xy.astype(np.float32).tofile(prefix + ".xy.raw")
    seg.astype(np.int32).tofile(prefix + ".seg.raw")
    attr.astype(np.float32).tofile(prefix + ".attr.raw")


def dump_master(prefix, frames):
    body = "".join(FRAME.format(time=t, nseg=ns, nvert=nv, name=n)
                   for (n, ns, nv, t) in frames if ns)
    with open(prefix + ".xdmf2", "w") as f:
        f.write('''\
<Xdmf
    Version="2">
  <Domain>
    <Grid
        GridType="Collection"
        CollectionType="Temporal">
''' + body + '''\
    </Grid>
  </Domain>
</Xdmf>
''')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("frames", nargs="+")
    ap.add_argument("--iso", type=float, default=0.0)
    ap.add_argument("--field", default=None)
    ap.add_argument("--out", default="iso")
    a = ap.parse_args()

    grids = [g for f in sorted(a.frames) for g in leaf_grids(f)]
    frames = []
    for i, (ncells, geom, attrs, time) in enumerate(grids):
        name = a.field or next(iter(attrs))
        if name not in attrs:
            raise SystemExit("field '%s' not in %s" % (name, list(attrs)))
        coords = read_raw(geom, ncells, 8).reshape(ncells, 4, 2).astype(np.float32)
        scalar = read_raw(attrs[name], ncells, 1).astype(np.float32)
        xy, seg, attr = amriso.extract2d(coords, scalar, scalar, a.iso)
        prefix = "%s.%04d" % (a.out, i)
        dump_raws(prefix, xy, seg, attr)
        frames.append((os.path.basename(prefix), len(seg), len(xy), time))
        print("%s [%s] iso=%g: %d segments, %d verts"
              % (os.path.basename(geom), name, a.iso, len(seg), len(xy)))

    dump_master(a.out, frames)
    print("-> %s.xdmf2 (%d frames)" % (a.out, len(frames)))


if __name__ == "__main__":
    main()
