/*
 *  I2D_ImposedAirfoil.cpp
 *  I2D_ROCKS
 *
 *  Created by Silvio Tödtli on 7/3/12.
 *  Copyright 2012 ETHZ. All rights reserved.
 *
 */

// TODO: implement calculation of m and J for discretized wing

#include "I2D_ImposedAirfoil.h"
#include <iostream>
#include <fstream>
#include <limits>

namespace AirfoilStuff {

	/**
	 * creates a non-rotated NACA profile at (0,0), the leading edge is the center of rotation
	 * geometry can be specified by the chord length and the NACA digits d1, d2, d3d4
	 *
	 */
	class NACA4gen
		{
		private:
			unsigned int nPanels;
			int maxDiffX, maxDiffEpsilon;
			vector<Real> x,yt,xu,yu,xl,yl,yC;
			// maxBoxX: stores physical coordinate of maximal base point in x; the others: analogue
			Real chord, maxBoxX, minBoxX, maxBoxY, minBoxY,leadingEdgeRadius,epsilon;
			Real leadingEdgeCenter[2];
			const int d1,d2,d3d4; // NACA digits: NACA d1-d2-d3d4

		public:

			Real getMaxX(void) const { return maxBoxX; }
			Real getMinX(void) const { return minBoxX; }
			Real getMaxY(void) const { return maxBoxY; }
			Real getMinY(void) const { return minBoxY; }
			Real getLeadingEdgeRadius(void){ return leadingEdgeRadius; }

			Real getNoseDistance(Real xx[2]) const
			{
				assert(xx[0]>=-2.0*epsilon && xx[0]<0.05*leadingEdgeRadius);

				Real dist = sqrt( (leadingEdgeCenter[0]-xx[0])*(leadingEdgeCenter[0]-xx[0]) + (leadingEdgeCenter[1]-xx[1])*(leadingEdgeCenter[1]-xx[1]) ) - leadingEdgeRadius;

				return dist;
			}


			/**
			 * - @ return: interpolates y coordinate of upper airfoil shape at position xx[0]
			 * - @ dist: calculates minimal distance between point xx and the upper sample points
			 *
			 * input:	xx[2]  	arbitrary point
			 *			dist	variable, where minimal distance is written to
			 */
			Real getYUpperProfile(Real xx[2], Real & dist) const
			{
				// nPanels = number of base points --> prerequisite: nPanels == xu.size == yu.size
				assert( xu.size() == yu.size() );
				assert( xu.size() == nPanels );
				assert(xx[0]>=0.0 && xx[0]<=chord);

				const Real dx = chord/((Real)nPanels-1.0); // interval length
				const int midIdx = floor(xx[0]/dx); // get base point before sought x-coordinate x[0]
				// base points to iterate over
				const int startIdx = ((midIdx-maxDiffX)<0) ? 0 : (midIdx-maxDiffX); // symmetrical profiles: maxDiff = 10
				const int endIdx = ((midIdx+maxDiffX)>(int)(nPanels))?nPanels:(midIdx+maxDiffX);

				assert(startIdx>=0);
				assert(endIdx>=0);
				assert(startIdx<=xu.size());
				assert(endIdx<=xu.size());
				assert(startIdx<endIdx);

				int lowIdx = 0;
				int highIdx = 0;
				for(int i=startIdx; i<endIdx; i++)
				{
					if(xu[i] > xx[0])
					{
						highIdx = i;
						lowIdx = i-1;
						break;
					}
				}

				assert(lowIdx>=0);
				assert(highIdx>=0);
				assert(lowIdx<=xu.size());
				assert(highIdx<=xu.size());
				assert(lowIdx<highIdx);

				// interpolate y value of the upper airfoil shape at position xx[0]
				const Real v = yu[lowIdx] + (yu[highIdx] - yu[lowIdx])/(xu[highIdx] - xu[lowIdx]) * (xx[0] - xu[lowIdx]);

				// base points to iterate over
				const int startIdxEps = ((lowIdx-maxDiffEpsilon)<0) ? 0 : (lowIdx-maxDiffEpsilon);
				const int endIdxEps = ((highIdx+maxDiffEpsilon)>(int)(nPanels))?nPanels:(highIdx+maxDiffEpsilon);

				assert(startIdxEps>=0);
				assert(endIdxEps>=0);
				assert(startIdxEps<=xu.size());
				assert(endIdxEps<=xu.size());
				assert(startIdxEps<endIdxEps);
				Real myDist = std::numeric_limits<Real>::max(); // ensure assignment in for loop (myDist for sure smaller than max value of Real)
				// calculate minimal distance between upper shape and point xx
				for(int i=startIdxEps; i<endIdxEps; i++)
				{
					const Real dxx = (xu[i]-xx[0]);
					const Real dyy = (yu[i]-xx[1]);
					myDist = min(myDist,(Real)sqrt(dxx*dxx+dyy*dyy));
				}
				dist = myDist;
				assert( dist!=std::numeric_limits<Real>::max() );

				return v;
			}


			/**
			 * - @ return: interpolates y coordinate of lower airfoil shape at position xx[0]
			 * - @ dist: calculates minimal distance between point xx and the lower sample points
			 *
			 * input:	xx[2]  	arbitrary point
			 *			dist	variable, where minimal distance is written to
			 */
			Real getYLowerProfile(Real xx[2], Real & dist) const
			{
				assert( xl.size() == yl.size() );
				assert( xl.size() == nPanels );
				assert(xx[0]>=0.0 && xx[0]<=chord);

				const Real dx = chord/((Real)nPanels-1.0);
				const int midIdx = floor(xx[0]/dx);
				const int startIdx = ((midIdx-maxDiffX)<0) ? 0 : (midIdx-maxDiffX);
				const int endIdx = ((midIdx+maxDiffX)>(int)(nPanels))?nPanels:(midIdx+maxDiffX);

				assert(startIdx>=0);
				assert(endIdx>=0);
				assert(startIdx<=xl.size());
				assert(endIdx<=xl.size());
				assert(startIdx<endIdx);

				int lowIdx = 0;
				int highIdx = 0;
				for(int i=startIdx; i<endIdx; i++)
				{
					if(xl[i] > xx[0])
					{
						highIdx = i;
						lowIdx = i-1;
						break;
					}
				}

				assert(lowIdx>=0);
				assert(highIdx>=0);
				assert(lowIdx<=xl.size());
				assert(highIdx<=xl.size());
				assert(lowIdx<highIdx);

				const Real v = yl[lowIdx] + (yl[highIdx] - yl[lowIdx])/(xl[highIdx] - xl[lowIdx]) * (xx[0] - xl[lowIdx]);

				const int startIdxEps = ((lowIdx-maxDiffEpsilon)<0) ? 0 : (lowIdx-maxDiffEpsilon);
				const int endIdxEps = ((highIdx+maxDiffEpsilon)>(int)(nPanels))?nPanels:(highIdx+maxDiffEpsilon);

				assert(startIdxEps>=0);
				assert(endIdxEps>=0);
				assert(startIdxEps<=xl.size());
				assert(endIdxEps<=xl.size());
				assert(startIdxEps<endIdxEps);
				Real myDist = std::numeric_limits<Real>::max();
				for(int i=startIdxEps; i<endIdxEps; i++)
				{
					const Real dxx = (xl[i]-xx[0]);
					const Real dyy = (yl[i]-xx[1]);
					myDist = min(myDist,(Real)sqrt(dxx*dxx+dyy*dyy));
				}
				dist = myDist;
				assert( dist!=std::numeric_limits<Real>::max() );

				return v;
			}


			/**
			 * @return: smallest distance between point xo and the airfoil silhouette
			 * 			distance negative if point inside body, positive if point outside body
			 *
			 * input:	xo[2]	arbitrary point
			 */
			Real sdf(Real xo[2]) const
			{
				// xo not within [0,chord]
				if (xo[0]<-2.0*epsilon){ return sqrt(pow(xo[0],2) + pow(xo[1], 2)); } // return distance to origin

				if (xo[0]>chord){ return sqrt(pow(xo[0]-chord,2) + pow(xo[1], 2)); } // return distance to trailing edge

				if (xo[0]>=-2.0*epsilon && xo[0]<0.05*leadingEdgeRadius){ return getNoseDistance(xo); } // return distance to leading edge

				// xo within [0,chord]
				Real distU = 0.0;
				Real distL = 0.0;
				const Real high = this->getYUpperProfile(xo,distU);
				const Real low = this->getYLowerProfile(xo,distL);

				assert(high >= low);

				const Real s = (Real)(xo[1]<low || xo[1]>high)*2.0 - 1.0; // s=-1, if xo inside body, s=1 if xo outside body
				const Real D = min(distU,distL);

				return s*D; // negative inside body, positive outside body
			}



			NACA4gen(Real chord_in,unsigned int nPanels_in,bool isFiniteTE,Real epsilon_in, const int _d1, const int _d2, const int _d3d4)
			:d1(_d1), d2(_d2), d3d4(_d3d4)
			{
				// 1. clear old basepoints
				x.clear();
				yt.clear();
				xu.clear();
				yu.clear();
				xl.clear();
				yl.clear();
				yC.clear();


				// 2. set up parameters
				printf("\nNACA-%d%d%d\n",d1,d2,d3d4);

				nPanels = nPanels_in;
				chord = chord_in;
				epsilon = epsilon_in;

				const Real t = (Real)d3d4/100.0;
				const Real m = (Real)d1/100.0;
				Real p = (Real)d2/10.0;

				// coefficients for symmetrical 4-digit NACA airfoil
				const Real a0 = 0.2969;
				const Real a1 = -0.1260;
				const Real a2 = -0.3516;
				const Real a3 = 0.2843;
				Real a4 = 0.0;

				if (isFiniteTE)
					a4=-0.1015; // For finite thick trailing edge
				else
					a4=-0.1036;  // For zero thick trailing edge


				// 3. calculate base points in [0,1]
				const Real dx = 1.0/((Real)nPanels-1.0);
				for(unsigned int i=0; i<nPanels; i++)
				{
					x.push_back( (Real)i*dx );
					yt.push_back( (t/0.2)*(a0*sqrt(x[i])+a1*x[i]+a2*x[i]*x[i]+a3*x[i]*x[i]*x[i]+a4*x[i]*x[i]*x[i]*x[i]) );
				}

				Real diffX = 0.0;

				// 3.a symmetrical airfoil
				if(p==0.0)
				{
					for(unsigned int i=0; i<nPanels; i++)
					{
						xu.push_back( x[i] );
						yu.push_back( yt[i] );

						xl.push_back( x[i] );
						yl.push_back( -yt[i] );

						yC.push_back(0.0);
					}
				}
				// 3.b non-symmetrical airfoil
				else
				{
					for(unsigned int i=0; i<nPanels; i++)
					{
						Real zc = 0.0;
						Real dyc_dx = 0.0;

						if(x[i]<=p)
						{
							zc=(m/(p*p))*((Real)2.0*p*x[i]-x[i]*x[i]);
							dyc_dx=(m/(p*p))*((Real)2.0*p-(Real)2.0*x[i]);
						}
						else
						{
							zc=(m/((1.0-p)*(1.0-p)))*((1.0-2.0*p)+2.0*p*x[i]-x[i]*x[i]);
							dyc_dx=(m/((1.0-p)*(1.0-p)))*(2.0*p-2.0*x[i]);
						}

						const Real theta=atan(dyc_dx);

						yC.push_back(zc);

						xu.push_back( x[i] - yt[i]*sin(theta) );
						diffX = max(diffX,(Real)fabs(x[i]-xu[i]));
						yu.push_back( zc + yt[i]*cos(theta) );

						xl.push_back( x[i] + yt[i]*sin(theta) );
						diffX = max(diffX,(Real)fabs(x[i]-xl[i]));
						yl.push_back( zc - yt[i]*cos(theta) );
					}
				}

				assert(nPanels==x.size());
				assert(nPanels==yt.size());
				assert(nPanels==xu.size());
				assert(nPanels==yu.size());
				assert(nPanels==xl.size());
				assert(nPanels==yl.size());
				assert(nPanels==yC.size());

				maxDiffX = ceil(diffX/dx)+10; // symmetrical profiles: maxDiffX = 10; ceil(x) returns the smallest int value not smaller than x
				maxDiffEpsilon = 2.0*ceil(epsilon/dx);


				// 4. scale basepoints from [0,1] to [0,chord]
				for(unsigned int i=0; i<nPanels; i++)
				{
					x[i] *= chord;
					yt[i] *= chord;
					xu[i] *= chord;
					yu[i] *= chord;
					xl[i] *= chord;
					yl[i] *= chord;
					yC[i] *= chord;
				}


				// 5. initialize some more stuff
				// assure assignment of one of the boxes (values are > 0 but smaller than max of Real for sure)
				maxBoxX = 0.0;
				minBoxX = std::numeric_limits<Real>::max(); // maximum value for a variable of type Real
				maxBoxY = 0.0;
				minBoxY = std::numeric_limits<Real>::max();

				leadingEdgeRadius = 1.1019*chord*t*t; // general formula

				const Real noseAngle = atan2((yC[1]-yC[0]),(x[1]-x[0]));
				leadingEdgeCenter[0] = leadingEdgeRadius*cos(noseAngle);
				leadingEdgeCenter[1] = leadingEdgeRadius*sin(noseAngle);

				for(unsigned int i=0; i<nPanels; i++)
				{
					maxBoxX = std::max((Real)maxBoxX,(Real)xu[i]);
					maxBoxX = std::max((Real)maxBoxX,(Real)xl[i]);

					minBoxX = std::min((Real)minBoxX,(Real)xu[i]);
					minBoxX = std::min((Real)minBoxX,(Real)xl[i]);

					maxBoxY = std::max((Real)maxBoxY,(Real)yu[i]);
					minBoxY = std::min((Real)minBoxY,(Real)yl[i]);
				}
			}

			~NACA4gen()
			{
			}
		};


	/**
	 * Provides transformation matrices between domain [unrotated airfoil at (0,0)] and codomain [rotated airfoil at (tx, ty)]
	 * basically each point of the codomain is mapped on the domain, where all functions involving the shape are evaluated
	 */
	class DiscretizedWing
	{
		Real O2W[3][3], W2O[3][3]; // W2O = inv(O2W)
		NACA4gen naca;
		const Real scaling, width, epsilon, s;
		Real ca, sa;


		/**
		 * ensures that O2W and W2O are inverse matrices
		 */
		void check_mapping() const // check if O2W and W2O are inverse matrices
		{
			Real A[3][3];

			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
				{
					A[i][j] = 0;

					for(int d=0; d<3; d++)
						A[i][j] += O2W[i][d]* W2O[d][j];
				}

			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					if (i==j)
						assert(fabs(A[i][j]-1)<1e-5);
					else
						assert(fabs(A[i][j]-0)<1e-5);
				}
			}
		}


		/**
		 * mapping codomain --> domain
		 * @ xo[2]: map point xw[2] in domain and write domain coordinates to it
		 *
		 * input:	xw[2]	coordinates of arbitrary point in codomain
		 * 			xo[2]	variable where coordinates of xw[2] in domain are written to
		 */
		void _w2o(const Real xw[2], Real xo[2]) const
		{
			xo[0] = W2O[0][0]*xw[0] + W2O[0][1]*xw[1] + W2O[0][2];
			xo[1] = W2O[1][0]*xw[0] + W2O[1][1]*xw[1] + W2O[1][2];
		}


		/**
		 * mapping domain --> codomain
		 * @ xw[2]: map point xo[2] in codomain and write codomain coordinates to it
		 *
		 * input:	xo[2]	coordinates of arbitrary point in domain
		 * 			xw[2]	variable where coordinates of xo[2] in codomain are written to
		 */
		void _o2w(const Real xo[2], Real xw[2]) const
		{
			xw[0] = O2W[0][0]*xo[0] + O2W[0][1]*xo[1] + O2W[0][2];
			xw[1] = O2W[1][0]*xo[0] + O2W[1][1]*xo[1] + O2W[1][2];
		}

	public:
		Real tx, ty, angle, vx, vy, angular_velocity, rho, m, J; // (tx,ty): position of airfoil center

		DiscretizedWing(Real xm, Real ym, Real width, Real angle_deg, Real epsilon, const int d1, const int d2, const int d3d4):
		naca(width,5000,false,epsilon,d1,d2,d3d4),
		scaling(1.0),
		epsilon(epsilon),
		width(width),
		tx(xm), ty(ym), angle(angle_deg), s(1.0), vx(0.0), vy(0.0), angular_velocity(0.0), rho(1.0), m(0.0), J(0.0)
		{
			ca = cos(angle*M_PI/180.);
			sa = sin(angle*M_PI/180.);

			// transformation matrix codomain --> domain
			W2O[0][0] = +1/s*ca;
			W2O[0][1] = +1/s*sa;
			W2O[0][2] = +1/s*(-tx*ca - ty*sa);
			W2O[1][0] = -1/s*sa;
			W2O[1][1] = +1/s*ca;
			W2O[1][2] = +1/s*(+tx*sa - ty*ca);
			W2O[2][0] = 0;
			W2O[2][1] = 0;
			W2O[2][2] = 1;

			// transformation matrix domain --> codomain
			O2W[0][0] = +s*ca;
			O2W[0][1] = -s*sa;
			O2W[0][2] = tx;
			O2W[1][0] = s*sa;
			O2W[1][1] = s*ca;
			O2W[1][2] = ty;
			O2W[2][0] = 0;
			O2W[2][1] = 0;
			O2W[2][2] = 1;

			check_mapping();
		}


		/**
		 * updates tx, ty, angle and the transformation matrices
		 */
		void update_all (const double dt, const double t)
		{
			tx += vx*dt;
			ty += vy*dt;
			angle += angular_velocity*dt;

			ca = cos(angle*M_PI/180.);
			sa = sin(angle*M_PI/180.);

			// transformation matrix codomain --> domain
			W2O[0][0] = +1/s*ca;
			W2O[0][1] = +1/s*sa;
			W2O[0][2] = +1/s*(-tx*ca - ty*sa);
			W2O[1][0] = -1/s*sa;
			W2O[1][1] = +1/s*ca;
			W2O[1][2] = +1/s*(+tx*sa - ty*ca);
			W2O[2][0] = 0;
			W2O[2][1] = 0;
			W2O[2][2] = 1;

			// transformation matrix domain --> codomain
			O2W[0][0] = +s*ca;
			O2W[0][1] = -s*sa;
			O2W[0][2] = tx;
			O2W[1][0] = s*sa;
			O2W[1][1] = s*ca;
			O2W[1][2] = ty;
			O2W[2][0] = 0;
			O2W[2][1] = 0;
			O2W[2][2] = 1;

			check_mapping();
		}


		/**
		 * @ return: smallest distance between point xw[2] and the airfoil silhouette
		 * 			 distance negative if point inside body, positive if point outside body
		 *
		 * input	 xw[2]	arbitrary point
		 */
		Real sdf(Real xw[2]) const
		{
			Real xo[2];
			_w2o(xw, xo);
			return scaling*naca.sdf(xo);
		}


		/**
		 * get edges of box around airfoil profile
		 *
		 * input:	eps		some kind of tolerance (?)
		 * 			xmin[2]	= (xmin, ymin) --> lower box edges
		 * 			xmax[2] = (xmax, ymax) --> upper box edges
		 */
		void bbox(const Real eps, Real xmin[2], Real xmax[2]) const
		{
			assert(eps>0);

			const Real v1[2] = {naca.getMinX(),naca.getMinY()};
			const Real v2[2] = {naca.getMinX(),naca.getMaxY()};
			const Real v3[2] = {naca.getMaxX(),naca.getMinY()};
			const Real v4[2] = {naca.getMaxX(),naca.getMaxY()};

			Real xw1[2], xw2[2], xw3[2], xw4[2];

			_o2w(v1, xw1);
			_o2w(v2, xw2);
			_o2w(v3, xw3);
			_o2w(v4, xw4);

			xmin[0] = min((Real)min((Real)min((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]); // xmin
			xmin[1] = min((Real)min((Real)min((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]); // ymin
			xmax[0] = max((Real)max((Real)max((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]); // xmax
			xmax[1] = max((Real)max((Real)max((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]); // ymax

			xmin[0] -= 2*eps;
			xmin[1] -= 2*eps;
			xmax[0] += 2*eps;
			xmax[1] += 2*eps;

			assert(xmin[0]<=xmax[0]);
			assert(xmin[1]<=xmax[1]);
		}

		void save(FILE * f) const
		{
			fprintf(f, "xm: %20.20e\n", tx);
			fprintf(f, "ym: %20.20e\n", ty);
			fprintf(f, "vx: %20.20e\n", vx);
			fprintf(f, "vy: %20.20e\n", vy);
			fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
			fprintf(f, "angle: %20.20e\n", angle);
		}

		void restart(FILE * f)
		{
			float val;

			fscanf(f, "xm: %e\n", &val);
			tx = val;
			printf("DiscretizedWing::restart(): xm is %e\n", tx);

			fscanf(f, "ym: %e\n", &val);
			ty = val;
			printf("DiscretizedWing::restart(): ym is %e\n", ty);

			fscanf(f, "vx: %e\n", &val);
			vx = val;
			printf("DiscretizedWing::restart(): vx is %e\n", vx);

			fscanf(f, "vy: %e\n", &val);
			vy = val;
			printf("DiscretizedWing::restart(): vy is %e\n", vy);

			fscanf(f, "angular_velocity: %e\n", &val);
			angular_velocity = val;
			printf("DiscretizedWing::restart(): angular_velocity is %e\n", angular_velocity);

			fscanf(f, "angle: %e\n", &val);
			angle = val;
			printf("DiscretizedWing::restart(): angle is %e\n", angle);
		}

		Real _mollified_heaviside(const Real x, const Real eps) const
		{
			const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
			return 0.5+0.5*cos(alpha);
		}


		/*
		 * calculate characteristic function at point p[2]
		 *
		 * @ return:	1 		inside body
		 * 				0		outside body
		 * 				[0,1]	at boundary (mollified heaviside)
		 */
		Real sample(Real p[2], const Real eps) const // calculate characteristic function at point (p[0],p[1])
		{
			Real dist = sdf(p);
			return _mollified_heaviside(dist,eps);
		}
	};


	/**
	 * Identical to FillBlocks of the other fluid mediated obstacles (Imposed Cylinder, Imposed Ellipse, ...)
	 */
	struct FillBlocks
	{
		Real eps;
		AirfoilStuff::DiscretizedWing * shape;

		FillBlocks(Real eps, AirfoilStuff::DiscretizedWing *_shape): eps(eps) { shape = _shape; }

		FillBlocks(const FillBlocks& c): eps(c.eps) { shape = c.shape; }

		static bool _is_touching(Real eps, const AirfoilStuff::DiscretizedWing * wing, const BlockInfo& info)
		{
			Real wing_min[2], wing_max[2];
			wing->bbox(eps, wing_min, wing_max); // store shapes' xmin, ymin in wing_min and shapes' xmax, ymax in wing_max

			Real min_pos[2], max_pos[2];
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);

			Real intersection[2][2] = {
				max(min_pos[0], wing_min[0]), min(max_pos[0],  wing_max[0]),
				max(min_pos[1], wing_min[1]), min(max_pos[1],  wing_max[1])
			};

			return
			intersection[0][1]-intersection[0][0]>0 &&
			intersection[1][1]-intersection[1][0]>0 ;
		}

		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			bool bEmpty = true;

			if(_is_touching(eps, shape, info))
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);

						b(ix,iy).tmp = max( shape->sample(p, eps), b(ix,iy).tmp );
					}

				bEmpty = false;
			}
		}
	};


	/**
	 * Identical to GetNonEmpty of the other fluid mediated obstacles (Imposed Cylinder, Imposed Ellipse, ...)
	 */
	struct GetNonEmpty
	{
		Real eps;
		AirfoilStuff::DiscretizedWing * shape;
		map<int, bool>& b2nonempty;

		GetNonEmpty(Real eps, AirfoilStuff::DiscretizedWing *_shape, map<int, bool>& b2nonempty): eps(eps), b2nonempty(b2nonempty)
		{
			shape = _shape;
		}

		GetNonEmpty(const GetNonEmpty& c): eps(c.eps), b2nonempty(c.b2nonempty)
		{
			shape = c.shape;
		}

		// specify in b which blocks are empty and which ones are not
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			bool bNonEmpty = false;

			if(FillBlocks::_is_touching(eps, shape, info)) // does shape lie (partially) within the block?
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);

						const Real Xs = shape->sample(p, eps); // calculate characteristic function (outside body 0, inside body 1, between: mollified heaviside
						bNonEmpty |= Xs>0; // bNonEmpty is true, if it was true before or if Xs>0 or if both are true
					}

				assert(b2nonempty.find(info.blockID) != b2nonempty.end());

				b2nonempty[info.blockID] = bNonEmpty;
			}
		}
	};


	/**
	 * Identical to FillVelblocks of the other fluid mediated obstacles (Imposed Cylinder, Imposed Ellipse, ...)
	 */
	struct FillVelblocks
	{
		double eps;
		AirfoilStuff::DiscretizedWing * shape;
		vector<pair< BlockInfo, VelocityBlock *> >& workitems;

		FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double eps, AirfoilStuff::DiscretizedWing *_shape):
		workitems(workitems), eps(eps)
		{
			shape = _shape;
		}

		FillVelblocks(const FillVelblocks& c): workitems(c.workitems), eps(c.eps)
		{
			shape = c.shape;
		}

		void operator()(blocked_range<int> range) const
		{
			for(int iblock=range.begin(); iblock<range.end(); iblock++)
			{
				BlockInfo info = workitems[iblock].first;
				VelocityBlock * u_desired = workitems[iblock].second;

				const Real xm = shape->tx;
				const Real ym =	shape->ty;

				const Real av = shape->angular_velocity;
				const Real vx = shape->vx;
				const Real vy = shape->vy;

				if(FillBlocks::_is_touching(eps, shape, info))
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);

							const Real Xs = shape->sample(p, eps);

							if (Xs > 0)
							{
								u_desired->u[0][iy][ix] = - av*(p[1]-ym) + vx;
								u_desired->u[1][iy][ix] = + av*(p[0]-xm) + vy;
							}
							else
							{
								u_desired->u[0][iy][ix] = 0.0;
								u_desired->u[1][iy][ix] = 0.0;
							}
						}
			}
		}
	};
}


// TODO: why is D and not _D used within constructor, why is Uinf assigned again?
I2D_ImposedAirfoil::I2D_ImposedAirfoil(Grid<W,B>& grid, ArgumentParser& parser, const Real _xm, const Real _ym, const Real _D, const Real _angle, const Real _vx, const Real _vy, const int d1, const int d2, const int d3d4, const Real _eps, const Real Uinf[2], I2D_PenalizationOperator& penalization)
:I2D_FloatingObstacleOperator(parser, grid, D, _eps, Uinf, penalization), vx_imposed(_vx), vy_imposed(_vy)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
	wing = new AirfoilStuff::DiscretizedWing(_xm,_ym,_D,_angle,eps,d1,d2,d3d4); // get rotated NACA profile
}

I2D_ImposedAirfoil::~I2D_ImposedAirfoil()
{
	assert(wing!=NULL);
	if(wing!=NULL)
	{
		delete wing;
		wing = NULL;
	}
}



//******************************************************************//
// below the pure virtual functions provided by ObstacleOperator    //
// or FloatingObstacleOperator are implemented                      //
//******************************************************************//

void I2D_ImposedAirfoil::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	AirfoilStuff::FillBlocks fill(eps,wing);
	block_processing.process(vInfo, coll, fill);
}

void I2D_ImposedAirfoil::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	wing->update_all(dt,t);

	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedAirfoil.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}

	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", this->dimT, t, wing->tx, wing->ty, wing->vx, wing->vy, wing->angle, wing->angular_velocity, wing->J, wing->m, wing->rho, this->Cd);
	fflush(ppFile);
	// Close file
	fclose(ppFile);
}

void I2D_ImposedAirfoil::save(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_ImposedAirfoil.txt", "w");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "w");
		assert(ppFile!=NULL);
	}

	// Actual save
	wing->save(ppFile);

	// Close file
	fclose(ppFile);
}

void I2D_ImposedAirfoil::restart(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_ImposedAirfoil.txt", "r");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "r");
		assert(ppFile!=NULL);
	}

	// Actual restart
	wing->restart(ppFile);

	// Close file
	fclose(ppFile);
}

void I2D_ImposedAirfoil::refresh(const double t, string filename)
{
	// Open file stream
	ifstream filestream;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		filestream.open("restart_I2D_ImposedAirfoil.txt");
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}
	else
	{
		filestream.open(filename.c_str());
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}

	// Open file
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();

	FILE * f = NULL;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_ImposedAirfoil.txt", "r");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "r");
		assert(f!=NULL);
	}

	//data stored in dataVect until t=t_restart
	vector < vector<float> > dataVect;
	for(int i=0; i<c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f,"%f",&variable);
		const Real dimensionlessTime = variable;
		fscanf(f,"%f",&variable);
		const Real time = variable;
		if(time <= t)
		{
			row.push_back(dimensionlessTime);
			row.push_back(time);
			fscanf(f,"%f",&variable);
			const Real xm = variable;
			row.push_back(xm);
			fscanf(f,"%f",&variable);
			const Real ym = variable;
			row.push_back(ym);
			fscanf(f,"%f",&variable);
			const Real vx = variable;
			row.push_back(vx);
			fscanf(f,"%f",&variable);
			const Real vy = variable;
			row.push_back(vy);
			fscanf(f,"%f",&variable);
			const Real angle = variable;
			row.push_back(angle);
			fscanf(f,"%f",&variable);
			const Real angular_velocity = variable;
			row.push_back(angular_velocity);
			fscanf(f,"%f",&variable);
			const Real J = variable;
			row.push_back(J);
			fscanf(f,"%f",&variable);
			const Real m = variable;
			row.push_back(m);
			fscanf(f,"%f",&variable);
			const Real rho = variable;
			row.push_back(rho);
			fscanf(f,"%f",&variable);
			const Real cD = variable;
			row.push_back(cD);
		}
		else
		{
			break;
		}
		dataVect.push_back(row);
	}
	fclose(f);
	f = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_ImposedAirfoil.txt", "w");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "w");
		assert(f!=NULL);
	}

	// Print stuff before t<restart
	for(unsigned int i=0; i<dataVect.size(); i++)
	{
		for(unsigned int j=0; j<dataVect[i].size(); j++)
			fprintf(f, "%e ", dataVect[i][j]);

		fprintf(f, "\n");
	}

	fclose(f);
}

void I2D_ImposedAirfoil::_setMotionPattern(const Real t)
{
	wing->vx = vx_imposed;
	wing->vy = vy_imposed;
	wing->angular_velocity = 0.0;
}

void I2D_ImposedAirfoil::computeDesiredVelocity(const double t)
{
	_setMotionPattern(t);

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(1);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	// Get non-empty blocks
	AirfoilStuff::GetNonEmpty getNonEmpty(eps, wing, nonempty);
	block_processing.process(vInfo, coll, getNonEmpty);

	// Set desired velocities
	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
		it->second = NULL;
	}
	desired_velocity.clear();

	vector<pair< BlockInfo, VelocityBlock *> > velblocks;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		if(nonempty[it->blockID] == true)
		{
			VelocityBlock * velblock = VelocityBlock::allocate(1);
			desired_velocity[it->blockID] = velblock;
			velblocks.push_back(pair< BlockInfo, VelocityBlock *>(*it, velblock));
		}
	}

	AirfoilStuff::FillVelblocks fillvelblocks(velblocks, eps, wing);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}
