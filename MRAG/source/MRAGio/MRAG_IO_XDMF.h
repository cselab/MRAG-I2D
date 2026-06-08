#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include "MRAGcore/MRAGCommon.h"

namespace MRAG
{
	template<typename TWavelets, typename TBlock>
	class IO_XDMF
	{
		static std::string _basename(const std::string& p)
		{
			const size_t s = p.find_last_of('/');
			return s == std::string::npos ? p : p.substr(s + 1);
		}

	public:
		struct Frame
		{
			double      time;
			long        ncells;
			std::string geom;
			std::string attr;
			std::string channel;
		};

		template<typename GridType>
		Frame Write(GridType& grid, const std::string& fileName,
		            const char* channelName = "phi", double time = 0.0, int step = 0)
		{
			std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
			const int  bx = TBlock::sizeX, by = TBlock::sizeY;
			const long ncells = (long)vInfo.size() * bx * by;

			const std::string xyzPath  = fileName + ".xyz.raw";
			const std::string attrPath = fileName + "." + channelName + ".raw";

			const std::string xdmfPath = fileName + ".xdmf2";
			FILE* xdmf = fopen(xdmfPath.c_str(), "w");
			if (xdmf == NULL) { printf("IO_XDMF: cannot open %s\n", xdmfPath.c_str()); return; }
			fprintf(xdmf,
			        "<Xdmf\n    Version=\"2.0\">\n"
			        "  <Domain>\n"
			        "    <Grid>\n"
			        "      <Time Value=\"%.16e\"/>\n"
			        "      <Information Name=\"Step\" Value=\"%d\"/>\n"
			        "      <Topology\n          Dimensions=\"%ld\"\n          TopologyType=\"Quadrilateral\"/>\n"
			        "      <Geometry\n          GeometryType=\"XY\">\n"
			        "        <DataItem\n            Dimensions=\"%ld 2\"\n            Format=\"Binary\">\n"
			        "          %s\n"
			        "        </DataItem>\n"
			        "      </Geometry>\n"
			        "      <Attribute\n          AttributeType=\"Scalar\"\n          Name=\"%s\"\n          Center=\"Cell\">\n"
			        "        <DataItem\n            Dimensions=\"%ld 1\"\n            Precision=\"%ld\"\n            Format=\"Binary\">\n"
			        "          %s\n"
			        "        </DataItem>\n"
			        "      </Attribute>\n"
			        "    </Grid>\n  </Domain>\n</Xdmf>\n",
			        time, step, ncells, 4 * ncells, _basename(xyzPath).c_str(),
			        channelName, ncells, (long)sizeof(Real), _basename(attrPath).c_str());
			fclose(xdmf);

			const double off = TWavelets::CenteringOffset;
			FILE* fxyz = fopen(xyzPath.c_str(), "wb");
			for (int i = 0; i < (int)vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				const double h = info.h[0];
				for (int iy = 0; iy < by; iy++)
					for (int ix = 0; ix < bx; ix++)
					{
						double x[2];
						info.pos(x, ix, iy);
						const float u0 = x[0] - off * h, v0 = x[1] - off * h;
						const float u1 = u0 + h, v1 = v0 + h;
						const float q[8] = { u0,v0, u0,v1, u1,v1, u1,v0 };
						fwrite(q, sizeof(float), 8, fxyz);
					}
			}
			fclose(fxyz);

			FILE* fattr = fopen(attrPath.c_str(), "wb");
			for (int i = 0; i < (int)vInfo.size(); i++)
			{
				BlockInfo& info  = vInfo[i];
				TBlock&    block = grid.getBlockCollection()[info.blockID];
				for (int iy = 0; iy < by; iy++)
					for (int ix = 0; ix < bx; ix++)
					{
						const Real v = (Real)block(ix, iy);
						fwrite(&v, sizeof(Real), 1, fattr);
					}
			}
			fclose(fattr);

			printf("wrote %s (+ .xyz.raw, .%s.raw), %ld cells\n",
			       xdmfPath.c_str(), channelName, ncells);

			Frame f;
			f.time = time; f.ncells = ncells;
			f.geom = _basename(xyzPath); f.attr = _basename(attrPath); f.channel = channelName;
			return f;
		}

		static void WriteTemporalMaster(const std::string& masterName,
		                                const std::vector<Frame>& frames)
		{
			const std::string path = masterName + ".xdmf2";
			FILE* m = fopen(path.c_str(), "w");
			if (m == NULL) {
				fprintf(stderr, "IO_XDMF: cannot open %s\n", path.c_str());
				return;
			}
			fprintf(m,
			        "<Xdmf\n"
			        "    Version=\"2\">\n"
			        "  <Domain>\n"
			        "    <Grid\n"
			        "        GridType=\"Collection\"\n"
			        "        CollectionType=\"Temporal\">\n");
			for (size_t i = 0; i < frames.size(); i++) {
				const Frame& f = frames[i];
				fprintf(m,
				        "      <Grid>\n"
				        "        <Time\n"
				        "            Value=\"%.16e\"/>\n"
				        "        <Topology\n"
				        "            TopologyType=\"Quadrilateral\"\n"
				        "            Dimensions=\"%ld\"/>\n"
				        "        <Geometry\n"
				        "            GeometryType=\"XY\">\n"
				        "          <DataItem\n"
				        "              Dimensions=\"%ld 2\"\n"
				        "              Format=\"Binary\">\n"
				        "            %s\n"
				        "          </DataItem>\n"
				        "        </Geometry>\n"
				        "        <Attribute\n"
				        "            AttributeType=\"Scalar\"\n"
				        "            Name=\"%s\"\n"
				        "            Center=\"Cell\">\n"
				        "          <DataItem\n"
				        "              Dimensions=\"%ld 1\"\n"
				        "              Precision=\"%ld\"\n"
				        "              Format=\"Binary\">\n"
				        "            %s\n"
				        "          </DataItem>\n"
				        "        </Attribute>\n"
				        "      </Grid>\n",
				        f.time, f.ncells, f.ncells, f.geom.c_str(),
				        f.channel.c_str(), f.ncells, (long)sizeof(Real), f.attr.c_str());
			}
			fprintf(m,
			        "    </Grid>\n"
			        "  </Domain>\n"
			        "</Xdmf>\n");
			fclose(m);
			printf("wrote %s (temporal master, %zu frames)\n", path.c_str(), frames.size());
		}
	};
}
