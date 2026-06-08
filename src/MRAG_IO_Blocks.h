#pragma once
#include <cstdio>
#include <string>
#include <vector>
#include "MRAGCommon.h"

namespace MRAG {

template <typename TWavelets, typename TBlock> struct IO_Blocks {
  struct Frame {
    double time;
    long nblocks;
    std::string geom, attr;
  };

  static std::string _basename(const std::string &p) {
    size_t s = p.find_last_of("/\\");
    return s == std::string::npos ? p : p.substr(s + 1);
  }

  template <typename GridType>
  static Frame Write(GridType &grid, const std::string &name, double time) {
    std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const long nb = vInfo.size();

    FILE *fxyz = fopen((name + ".xyz.raw").c_str(), "wb");
    FILE *flev = fopen((name + ".level.raw").c_str(), "wb");
    for (long i = 0; i < nb; i++) {
      const BlockInfo &b = vInfo[i];
      const float x0 = b.origin[0], y0 = b.origin[1];
      const float x1 = x0 + TBlock::sizeX * b.h[0];
      const float y1 = y0 + TBlock::sizeY * b.h[1];
      const float q[8] = {x0, y0, x0, y1, x1, y1, x1, y0};
      fwrite(q, sizeof(float), 8, fxyz);
      const float l = (float)b.level;
      fwrite(&l, sizeof(float), 1, flev);
    }
    fclose(fxyz);
    fclose(flev);

    Frame f;
    f.time = time;
    f.nblocks = nb;
    f.geom = _basename(name + ".xyz.raw");
    f.attr = _basename(name + ".level.raw");
    return f;
  }

  static void WriteTemporalMaster(const std::string &name,
                                  const std::vector<Frame> &frames) {
    FILE *x = fopen((name + ".xdmf2").c_str(), "w");
    fprintf(x, "<Xdmf\n    Version=\"2\">\n  <Domain>\n"
               "    <Grid\n        GridType=\"Collection\"\n"
               "        CollectionType=\"Temporal\">\n");
    for (size_t i = 0; i < frames.size(); i++) {
      const Frame &f = frames[i];
      fprintf(x,
              "      <Grid>\n"
              "        <Time\n            Value=\"%.16e\"/>\n"
              "        <Topology\n            TopologyType=\"Quadrilateral\"\n"
              "            Dimensions=\"%ld\"/>\n"
              "        <Geometry\n            GeometryType=\"XY\">\n"
              "          <DataItem\n              Dimensions=\"%ld 2\"\n"
              "              Precision=\"4\"\n              Format=\"Binary\">\n"
              "            %s\n          </DataItem>\n"
              "        </Geometry>\n"
              "        <Attribute\n            Center=\"Cell\"\n"
              "            Name=\"level\">\n"
              "          <DataItem\n              Dimensions=\"%ld\"\n"
              "              Precision=\"4\"\n              Format=\"Binary\">\n"
              "            %s\n          </DataItem>\n"
              "        </Attribute>\n"
              "      </Grid>\n",
              f.time, f.nblocks, 4 * f.nblocks, f.geom.c_str(), f.nblocks,
              f.attr.c_str());
    }
    fprintf(x, "    </Grid>\n  </Domain>\n</Xdmf>\n");
    fclose(x);
  }
};

} // namespace MRAG
