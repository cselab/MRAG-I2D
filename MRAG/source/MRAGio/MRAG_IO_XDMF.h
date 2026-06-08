#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include "MRAGcore/MRAGCommon.h"

namespace MRAG
{
	template<typename E> inline Real _xdmf_get(E& e, int c, Real h) { return e.giveMe(c, h); }
	inline Real _xdmf_get(Real& e, int, Real) { return e; }

	template<typename TWavelets, typename TBlock, int nChannels = 1, int iChannelStart = 0>
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
			double time;
			long   ncells;
			std::string geom;
			std::vector<std::string> name;
			std::vector<std::string> file;
		};

		template<typename GridType>
		Frame Write(GridType& grid, const std::string& fileName,
		            double time = 0.0, int step = 0, const char* const* names = NULL)
		{
			std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
			const int  bx = TBlock::sizeX, by = TBlock::sizeY;
			const long ncells = (long)vInfo.size() * bx * by;

			Frame fr;
			fr.time = time; fr.ncells = ncells;
			fr.geom = _basename(fileName + ".xyz.raw");
			for (int c = 0; c < nChannels; c++) {
				char nm[64];
				if (names != NULL) snprintf(nm, sizeof nm, "%s", names[c]);
				else               snprintf(nm, sizeof nm, "channel%d", c);
				fr.name.push_back(nm);
				fr.file.push_back(_basename(fileName + "." + nm + ".raw"));
			}

			const std::string xdmfPath = fileName + ".xdmf2";
			FILE* xdmf = fopen(xdmfPath.c_str(), "w");
			if (xdmf == NULL) { printf("IO_XDMF: cannot open %s\n", xdmfPath.c_str()); return fr; }
			fprintf(xdmf,
			        "<Xdmf\n"
			        "    Version=\"2\">\n"
			        "  <Domain>\n"
			        "    <Grid>\n"
			        "      <Time\n"
			        "          Value=\"%.16e\"/>\n"
			        "      <Information\n"
			        "          Name=\"Step\"\n"
			        "          Value=\"%d\"/>\n"
			        "      <Topology\n"
			        "          TopologyType=\"Quadrilateral\"\n"
			        "          Dimensions=\"%ld\"/>\n"
			        "      <Geometry\n"
			        "          GeometryType=\"XY\">\n"
			        "        <DataItem\n"
			        "            Dimensions=\"%ld 2\"\n"
			        "            Format=\"Binary\">\n"
			        "          %s\n"
			        "        </DataItem>\n"
			        "      </Geometry>\n",
			        time, step, ncells, 4 * ncells, fr.geom.c_str());
			for (int c = 0; c < nChannels; c++)
				fprintf(xdmf,
				        "      <Attribute\n"
				        "          AttributeType=\"Scalar\"\n"
				        "          Name=\"%s\"\n"
				        "          Center=\"Cell\">\n"
				        "        <DataItem\n"
				        "            Dimensions=\"%ld 1\"\n"
				        "            Precision=\"%ld\"\n"
				        "            Format=\"Binary\">\n"
				        "          %s\n"
				        "        </DataItem>\n"
				        "      </Attribute>\n",
				        fr.name[c].c_str(), ncells, (long)sizeof(Real), fr.file[c].c_str());
			fprintf(xdmf,
			        "    </Grid>\n"
			        "  </Domain>\n"
			        "</Xdmf>\n");
			fclose(xdmf);

			const double off = TWavelets::CenteringOffset;
			FILE* fxyz = fopen((fileName + ".xyz.raw").c_str(), "wb");
			for (int i = 0; i < (int)vInfo.size(); i++) {
				BlockInfo& info = vInfo[i];
				const double h = info.h[0];
				for (int iy = 0; iy < by; iy++)
					for (int ix = 0; ix < bx; ix++) {
						double x[2]; info.pos(x, ix, iy);
						const float u0 = x[0] - off * h, v0 = x[1] - off * h, u1 = u0 + h, v1 = v0 + h;
						const float q[8] = { u0,v0, u0,v1, u1,v1, u1,v0 };
						fwrite(q, sizeof(float), 8, fxyz);
					}
			}
			fclose(fxyz);

			for (int c = 0; c < nChannels; c++) {
				FILE* fa = fopen((fileName + "." + fr.name[c] + ".raw").c_str(), "wb");
				for (int i = 0; i < (int)vInfo.size(); i++) {
					BlockInfo& info  = vInfo[i];
					TBlock&    block = grid.getBlockCollection()[info.blockID];
					const Real h = info.h[0];
					for (int iy = 0; iy < by; iy++)
						for (int ix = 0; ix < bx; ix++) {
							const Real v = _xdmf_get(block(ix, iy), c + iChannelStart, h);
							fwrite(&v, sizeof(Real), 1, fa);
						}
				}
				fclose(fa);
			}

			printf("wrote %s (+ geom + %d channel(s)), %ld cells\n", xdmfPath.c_str(), nChannels, ncells);
			return fr;
		}

		template<typename GridType, typename BInfo>
		Frame Write(GridType& grid, BInfo&, const std::string& fileName)
		{
			return Write(grid, fileName, 0.0, 0, NULL);
		}

		static void WriteTemporalMaster(const std::string& masterName, const std::vector<Frame>& frames)
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
				        "        </Geometry>\n",
				        f.time, f.ncells, 4 * f.ncells, f.geom.c_str());
				for (size_t c = 0; c < f.name.size(); c++)
					fprintf(m,
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
					        "        </Attribute>\n",
					        f.name[c].c_str(), f.ncells, (long)sizeof(Real), f.file[c].c_str());
				fprintf(m, "      </Grid>\n");
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
