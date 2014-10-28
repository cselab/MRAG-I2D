/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_WingObstacleOperator.h"
#include <limits>

namespace WingStuff {

	class NACA4gen
		{
		private:
			unsigned int nPanels;
			int maxDiffX, maxDiffEpsilon;
			vector<Real> x,yt,xu,yu,xl,yl,yC;
			// maxBoxX: stores physical coordinate of maximal base point in x; the others: analogue
			Real chord, maxBoxX, minBoxX, maxBoxY, minBoxY,leadingEdgeRadius,epsilon;
			Real leadingEdgeCenter[2];
			
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
				const Real D = min(distU,distL); // positive inside body, negative outside body
				
				return s*D;
			}
			
			
			
			NACA4gen(Real chord_in,unsigned int nPanels_in,bool isFiniteTE,Real epsilon_in)
			{
				x.clear();
				yt.clear();
				xu.clear();
				yu.clear();
				xl.clear();
				yl.clear();
				yC.clear();
				
				/*
				char cstrD1;
				char cstrD2;
				char cstrD3D4[3];
				string strD1;
				string strD2;
				string strD3D4;
				
				cstrD1 = wingType[0];
				cstrD2 = wingType[1];				
				cstrD3D4[0] = wingType[2];
				cstrD3D4[1] = wingType[3];
				
				strD1 = cstrD1;
				strD2 = cstrD2;
				strD3D4 = cstrD3D4;
								
				int d1 = atoi(strD1.c_str());
				int d2 = atoi(strD2.c_str());
				int d3d4 = atoi(strD3D4.c_str());
				 */
				
				const int d1 = 4;
				const int d2 = 4;
				const int d3d4 = 15;
								
				assert(d1>=0 && d1<=9);
				assert(d2>=0 && d2<=9);				
				assert(d3d4>=0 && d3d4<=99);
				
				printf("\nNACA-%d%d%d\n",d1,d2,d3d4);
				
				nPanels = nPanels_in;				
				chord = chord_in;
				epsilon = epsilon_in;
					
				const Real t = (Real)d3d4/100.0;
				const Real m = (Real)d1/100.0;
				Real p = (Real)d2/10.0;
				
				const Real a0 = 0.2969;
				const Real a1 = -0.1260;
				const Real a2 = -0.3516;
				const Real a3 = 0.2843;
				Real a4 = 0.0;
				
				if (isFiniteTE)
					a4=-0.1015; // For finite thick TE
				else
					a4=-0.1036;  // For zero thick TE
				
				const Real dx = 1.0/((Real)nPanels-1.0);
				for(unsigned int i=0; i<nPanels; i++)
				{
					x.push_back( (Real)i*dx );
					yt.push_back( (t/0.2)*(a0*sqrt(x[i])+a1*x[i]+a2*x[i]*x[i]+a3*x[i]*x[i]*x[i]+a4*x[i]*x[i]*x[i]*x[i]) );
				}
				
				Real diffX = 0.0;
				if(p==0.0) // symmetrical airfoil
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
				else // non-symmetrical airfoil
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
				
				// assure assignment of one of the boxes (values are > 0 but smaller than max of Real for sure)
				maxBoxX = 0.0;
				minBoxX = std::numeric_limits<Real>::max(); // maximum value for a variable of type Real
				maxBoxY = 0.0;
				minBoxY = std::numeric_limits<Real>::max();
				
				leadingEdgeRadius = 1.1019*chord*t*t; // general formula
				
				//std::cout << x[0] << std::endl;
				//std::cout << x[1] << std::endl;
				//std::cout << yC[0] << std::endl;
				//std::cout << yC[1] << std::endl;
				const Real noseAngle = atan2((yC[1]-yC[0]),(x[1]-x[0]));
				leadingEdgeCenter[0] = leadingEdgeRadius*cos(noseAngle);
				leadingEdgeCenter[1] = leadingEdgeRadius*sin(noseAngle);
				//std::cout << noseAngle << std::endl;
				//std::cout << leadingEdgeCenter[0] << std::endl;
				//std::cout << leadingEdgeCenter[1] << std::endl;
				
				for(unsigned int i=0; i<nPanels; i++)
				{
					maxBoxX = std::max((Real)maxBoxX,(Real)xu[i]);
					maxBoxX = std::max((Real)maxBoxX,(Real)xl[i]);
					
					minBoxX = std::min((Real)minBoxX,(Real)xu[i]);
					minBoxX = std::min((Real)minBoxX,(Real)xl[i]);
					
					maxBoxY = std::max((Real)maxBoxY,(Real)yu[i]);					
					minBoxY = std::min((Real)minBoxY,(Real)yl[i]);
				}
				
				//printf("maxBoxX %f\n",maxBoxX);
				//printf("minBoxX %f\n",minBoxX);
				//printf("maxBoxY %f\n",maxBoxY);
				//printf("minBoxY %f\n",minBoxY);				
			}
			
			~NACA4gen()
			{
			}
		};
	
	class DiscretizedWing
	{
		Real O2W[3][3], W2O[3][3]; // W2O = inv(O2W)
		NACA4gen naca;
		const Real scaling, width, epsilon;
		
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
		
		
		void _w2o(const Real xw[2], Real xo[2]) const
		{
			xo[0] = W2O[0][0]*xw[0] + W2O[0][1]*xw[1] + W2O[0][2];
			xo[1] = W2O[1][0]*xw[0] + W2O[1][1]*xw[1] + W2O[1][2];
		}
		
		void _o2w(const Real xo[2], Real xw[2]) const
		{
			xw[0] = O2W[0][0]*xo[0] + O2W[0][1]*xo[1] + O2W[0][2];
			xw[1] = O2W[1][0]*xo[0] + O2W[1][1]*xo[1] + O2W[1][2];
		}
		
	public:
		DiscretizedWing(Real xm, Real ym, Real width, Real angle_deg, Real epsilon):
		naca(width,5000,false,epsilon),
		scaling(1.0),
		epsilon(epsilon),
		width(width)
		{
			const Real tx = xm;
			const Real ty = ym;
			const Real ca = cos(angle_deg*M_PI/180.);
			const Real sa = sin(angle_deg*M_PI/180.);
			const Real s = 1.0;
			
			W2O[0][0] = +1/s*ca;
			W2O[0][1] = +1/s*sa;
			W2O[0][2] = +1/s*(-tx*ca - ty*sa);
			W2O[1][0] = -1/s*sa;
			W2O[1][1] = +1/s*ca;
			W2O[1][2] = +1/s*(+tx*sa - ty*ca);
			W2O[2][0] = 0;
			W2O[2][1] = 0;
			W2O[2][2] = 1;
			
			// rotation matrix
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
		
		void getAerodynamicCenter(Real cor[2]) const
		{
			cor[0] = 0.25*width;
			cor[1] = 0.0;
			
			Real corw[2] = {0.0,0.0};
			_o2w(cor, corw); // rotate cor according to angle
			
			cor[0] = corw[0];
			cor[1] = corw[1];
		}
		
		Real sdf(Real xw[2]) const
		{
			Real xo[2];			
			_w2o(xw, xo);
			return scaling*naca.sdf(xo);
		}
		
		void getTrailingBleedsLocations(int bleeIDX, Real & lenghtBleed, Real & widthBleed, Real & angleBleed, Real cmBleed[2], Real outletStart[2])
		{	
			Real inlet01 = 0.596*width;
			Real inlet02 = 0.816*width;
			
			Real outlet01 = (13.0/200.0)*width;
			Real outlet02 = (19.0/200.0)*width;
			
			Real inlet[2] = {0.0,0.0};
			Real outlet[2] = {0.0,0.0};
			widthBleed = 0.0165*width;
			
			assert(bleeIDX>=1 && bleeIDX<=2);
			
			if(bleeIDX==1)
			{
				inlet[0] = inlet01;
				outlet[0] = outlet01;
			}
			if(bleeIDX==2)
			{
				inlet[0] = inlet02;
				outlet[0] = outlet02;
			}
			
			Real dist = 0.0;
			Real low = naca.getYLowerProfile(inlet,dist);
			Real high = naca.getYUpperProfile(outlet,dist);
			
			inlet[1] = low;
			outlet[1] = high;
			
			Real inletw[2] = {0.0,0.0};
			Real outletw[2] = {0.0,0.0};			
			_o2w(inlet, inletw); // rotate according to angle
			_o2w(outlet, outletw);
			
			cmBleed[0] = inletw[0];
			cmBleed[1] = inletw[1];
			
			Real extension = 10.0*epsilon;
			
			Real lenghtTube = sqrt( (inletw[0]-outletw[0])*(inletw[0]-outletw[0]) + (inletw[1]-outletw[1])*(inletw[1]-outletw[1]) );
			
			lenghtBleed = lenghtTube;// + 2.0*extension;
			
			angleBleed = atan2( (outletw[1]-inletw[1]), (outletw[0]-inletw[0]) );
			
			cmBleed[0] = inletw[0];//-extension*cos(angleBleed);
			cmBleed[1] = inletw[1];//-extension*sin(angleBleed);
			
			outletStart[0] = outletw[0];// + (lenghtTube/2.0)*cos(angleBleed);
			outletStart[1] = outletw[1];// + (lenghtTube/2.0)*sin(angleBleed);
		}
		
		void getRearBleedsLocations(int bleeIDX, Real & lenghtBleed, Real & widthBleed, Real & angleBleed, Real cmBleed[2], Real outletStart[2])
		{	
			Real inlet01 = 0.596*width;
			Real inlet02 = 0.632*width;
			Real inlet03 = 0.665*width;
			Real inlet04 = 0.697*width;
			Real inlet05 = 0.729*width;
			Real inlet06 = 0.761*width;
			Real inlet07 = 0.789*width;
			Real inlet08 = 0.816*width;
			
			Real outlet01 = 0.764*width;
			Real outlet02 = 0.802*width;
			Real outlet03 = 0.837*width;
			
			Real inlet[2] = {0.0,0.0};
			Real outlet[2] = {0.0,0.0};
			widthBleed = 0.0165*width;
			
			assert(bleeIDX>=1 && bleeIDX<=8);
			
			if(bleeIDX==1)
			{
					inlet[0] = inlet01;
					outlet[0] = outlet01;
			}
			if(bleeIDX==2)
			{
				inlet[0] = inlet02;
				outlet[0] = outlet01;
			}
			if(bleeIDX==3)
			{
				inlet[0] = inlet03;
				outlet[0] = outlet02;
			}
			if(bleeIDX==4)
			{
				inlet[0] = inlet04;
				outlet[0] = outlet02;
			}
			if(bleeIDX==5)
			{
				inlet[0] = inlet05;
				outlet[0] = outlet02;
			}
			if(bleeIDX==6)
			{
				inlet[0] = inlet06;
				outlet[0] = outlet03;
			}
			if(bleeIDX==7)
			{
				inlet[0] = inlet07;
				outlet[0] = outlet03;
			}
			if(bleeIDX==8)
			{
				inlet[0] = inlet08;
				outlet[0] = outlet03;
			}

			Real dist = 0.0;			
			Real low = naca.getYLowerProfile(inlet,dist);
			Real high = naca.getYUpperProfile(outlet,dist);
			
			inlet[1] = low;
			outlet[1] = high;
						
			Real inletw[2] = {0.0,0.0};
			Real outletw[2] = {0.0,0.0};			
			_o2w(inlet, inletw);
			_o2w(outlet, outletw);
			
			cmBleed[0] = inletw[0];
			cmBleed[1] = inletw[1];
			
			Real extension = 10.0*epsilon;
			
			Real lenghtTube = sqrt( (inletw[0]-outletw[0])*(inletw[0]-outletw[0]) + (inletw[1]-outletw[1])*(inletw[1]-outletw[1]) );
			
			lenghtBleed = lenghtTube;// + 2.0*extension;
			
			angleBleed = atan2( (outletw[1]-inletw[1]), (outletw[0]-inletw[0]) );
			
			cmBleed[0] = inletw[0];//-extension*cos(angleBleed);
			cmBleed[1] = inletw[1];//-extension*sin(angleBleed);
			
			outletStart[0] = outletw[0];// + (lenghtTube/2.0)*cos(angleBleed);
			outletStart[1] = outletw[1];// + (lenghtTube/2.0)*sin(angleBleed);
		}
		
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
	};
	
	DiscretizedWing * wing;
	
	struct FillCylBleed
	{
		const Real smooth_radius, smoothing_length, support_radius;
		Real sphere_position[2];
		Real sphere_box[2][2];
		
		void _find_sphere_box()
		{
			sphere_box[0][0] = sphere_position[0] - support_radius;
			sphere_box[0][1] = sphere_position[0] + support_radius;
			sphere_box[1][0] = sphere_position[1] - support_radius;
			sphere_box[1][1] = sphere_position[1] + support_radius;
		}
		
		FillCylBleed(Real support_radius, Real smooth_radius, Real smoothing_length,  Real sphere_position[2]): 
		support_radius(support_radius), smooth_radius(smooth_radius), smoothing_length(smoothing_length)
		{
			this->sphere_position[0] = sphere_position[0];
			this->sphere_position[1] = sphere_position[1];
			
			_find_sphere_box();
		}
		
		FillCylBleed(const FillCylBleed& c): 
		support_radius(c.support_radius), smooth_radius(c.smooth_radius), smoothing_length(c.smoothing_length)
		{
			sphere_position[0] = c.sphere_position[0];
			sphere_position[1] = c.sphere_position[1];
			
			_find_sphere_box();
		}
		
		Real mollified_heaviside(const Real x, const Real eps) const
		{
			const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));			
			return 0.5+0.5*cos(alpha);
		}
		
		bool _is_touching(const BlockInfo& info) const // true if (at least a part of) the shape lies within the block
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0); // minimal physical coordinates of block
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1); // maximal physical coordinates of block
			
			Real intersection[2][2] = {
				max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
				max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
			};
			
			return 
			intersection[0][1]-intersection[0][0]>0 && 
			intersection[1][1]-intersection[1][0]>0 ;
		}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{		
			if(_is_touching(info))
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);
						
						const Real r = sqrt(pow(p[0]-sphere_position[0],2) + 
											pow(p[1]-sphere_position[1],2) );
						
						b(ix, iy).tmp *= 1.0-this->mollified_heaviside(r-smooth_radius, smoothing_length);
					}
			}
		}
	};
	
	struct FillBlocksBleed
	{
		Real m_length;
		Real m_width;
		Real m_epsilon;
		Real m_cm[2], v1[2], v2[2], v3[2], v4[2], nor1[2], nor2[2];
		Real m_angle;
		Real rectangle_box[2][2];
		
		FillBlocksBleed(Real length, Real width, Real epsilon,  Real cm[2], Real angle)
		{
			m_length = length;
			m_width = width;
			m_epsilon = epsilon;
			m_angle = angle;
			
			m_cm[0] = cm[0];
			m_cm[1] = cm[1];
			
			_compute_vertices(v1,v2,v3,v4,nor1,nor2,m_cm,m_length,m_width,m_angle);
			_find_rectangle_box();
		}
		
		FillBlocksBleed(const FillBlocksBleed& c):
		m_length(c.m_length), m_width(c.m_width), m_epsilon(c.m_epsilon), m_angle(c.m_angle)
		{
			m_cm[0] = c.m_cm[0];
			m_cm[1] = c.m_cm[1];
			
			v1[0] = c.v1[0];
			v1[1] = c.v1[1];
			v2[0] = c.v2[0];
			v2[1] = c.v2[1];
			v3[0] = c.v3[0];
			v3[1] = c.v3[1];
			v4[0] = c.v4[0];
			v4[1] = c.v4[1];
			nor1[0] = c.nor1[0];
			nor1[1] = c.nor1[1];
			nor2[0] = c.nor2[0];
			nor2[1] = c.nor2[1];
			
			rectangle_box[0][0] = c.rectangle_box[0][0];
			rectangle_box[0][1] = c.rectangle_box[0][1];
			rectangle_box[1][0] = c.rectangle_box[1][0];
			rectangle_box[1][1] = c.rectangle_box[1][1];
		}
		
		bool _is_touching(const BlockInfo& info) const
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real intersection[2][2] = {
				max(min_pos[0], rectangle_box[0][0]), min(max_pos[0], rectangle_box[0][1]),
				max(min_pos[1], rectangle_box[1][0]), min(max_pos[1], rectangle_box[1][1]),
			};
			
			return 
			intersection[0][1]-intersection[0][0]>0 && 
			intersection[1][1]-intersection[1][0]>0 ;
		}
		
		Real mollified_heaviside(const Real x, const Real eps) const
		{
			const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));	// x negative inside box, positive outside box
			return 0.5+0.5*cos(alpha); // outside: 0, inside: 1, boundary: mollified heavyside
		}
		
		void _find_rectangle_box()
		{
			rectangle_box[0][0] = min(min(min(v1[0],v2[0]),v3[0]),v4[0]) - m_epsilon; // xmin
			rectangle_box[0][1] = max(max(max(v1[0],v2[0]),v3[0]),v4[0]) + m_epsilon; // xmax
			rectangle_box[1][0] = min(min(min(v1[1],v2[1]),v3[1]),v4[1]) - m_epsilon; // ymin
			rectangle_box[1][1] = max(max(max(v1[1],v2[1]),v3[1]),v4[1]) + m_epsilon; // ymax
		}
		
		void _compute_norm(const Real *v2,const Real *v1,Real *nor1,Real *nor2)
		{
			nor1[0] = 1.0;
			nor1[1] = 0.0;
			Real v[2] = {0.0,0.0};
			v[0] = v2[0]-v1[0];
			v[1] = v2[1]-v1[1];
			Real mod2 = sqrt(v[0]*v[0]+v[1]*v[1]);
			nor2[0] = v[0]/mod2;
			nor2[1] = v[1]/mod2; // nor2: normalized vector in direction v1-v2

			// compute direction normal to v1-v2
			if(v[0]==0.0)
			{
				nor1[0] = 1.0;
				nor1[1] = 0.0;
				return;
			}
			if(v[1]==0.0)
			{
				nor1[0] = 0.0;
				nor1[1] = 1.0;
				return;
			}
			nor1[1] = -(v[0]*nor1[0])/v[1];
			Real mod1 = sqrt(nor1[0]*nor1[0]+nor1[1]*nor1[1]);
			nor1[0] /= mod1;
			nor1[1] /= mod1; // nor1: normalized vector in direction orthogonal to v1-v2
		}
		
		void _compute_vertices(Real *v1,Real *v2,Real *v3,Real *v4,Real *nor1,Real *nor2,const Real *cm,const Real length,const Real width,const Real angle)
		{
			nor1[0] = 0.0;
			nor1[1] = 0.0;
			nor2[0] = 0.0;
			nor2[1] = 0.0;
			Real cm1[2] = {cm[0],cm[1]};
			cm1[0] += length*cos(angle); // cm: "leading edge", cm1: "trailing edge"
			cm1[1] += length*sin(angle);
			_compute_norm(cm1,cm,nor1,nor2); // calculate specific block unit vectors nor1 and nor2
			v1[0] = cm[0] + nor1[0]*width/2.0;
			v1[1] = cm[1] + nor1[1]*width/2.0;
			v2[0] = cm[0] - nor1[0]*width/2.0;
			v2[1] = cm[1] - nor1[1]*width/2.0;			
			v3[0] = cm1[0] + nor1[0]*width/2.0;
			v3[1] = cm1[1] + nor1[1]*width/2.0;
			v4[0] = cm1[0] - nor1[0]*width/2.0;
			v4[1] = cm1[1] - nor1[1]*width/2.0;
		}
		
		Real _compute_char_func_box(const Real *v1,const Real *v2,const Real *v3,const Real *v4,const Real *nor1,const Real *nor2,const Real *pos,const Real length,const Real width) const
		{
			Real dist = 0.0;
			Real a[2] = {v1[0]-pos[0],v1[1]-pos[1]};
			Real b[2] = {v2[0]-pos[0],v2[1]-pos[1]};
			Real c[2] = {v2[0]-pos[0],v2[1]-pos[1]};
			Real d[2] = {v4[0]-pos[0],v4[1]-pos[1]};
			Real distA = fabs(a[0]*nor1[0]+a[1]*nor1[1]);
			Real distB = fabs(b[0]*nor1[0]+b[1]*nor1[1]);
			Real distC = fabs(c[0]*nor2[0]+c[1]*nor2[1]);
			Real distD = fabs(d[0]*nor2[0]+d[1]*nor2[1]);
			if(distA<=width && distB<=width) // y-pos inside y-range rectangular box
			{
				dist = min(distA,distB);
			}
			else // y-pos outside y-range rectangular box
			{
				dist = -min(distA,distB);
			}
			if(distC>length || distD>length) // x-pos outside x-range rectangular box
			{
				dist = -1e30;
			}
			return dist; // positive inside box, negative outside box
		}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const // calculates characteristic function
		{		
			if(_is_touching(info))
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);
						
						Real dist = _compute_char_func_box(v1,v2,v3,v4,nor1,nor2,p,m_length,m_width);
						
						b(ix,iy).tmp *= 1.0-this->mollified_heaviside(-dist,m_epsilon); // outside:1, inside: 0
					}
			}
		}
	};
	
	struct FillChiBlocks
	{
		Real eps;
		BlockInfo * vInfo;
		
		FillChiBlocks(){}
		FillChiBlocks(const FillChiBlocks& c): eps(c.eps), vInfo(c.vInfo){}
		
		Real mollified_heaviside(const Real x, const Real eps) const
		{
			const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
			
			return 0.5+0.5*cos(alpha);
		}
		
		bool _is_touching(const BlockInfo& info) const
		{
			Real wing_min[2], wing_max[2];
			WingStuff::wing->bbox(eps, wing_min, wing_max); // store shapes' xmin, ymin in wing_min and shapes' xmax, ymax in wing_max
			
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
			
			return true;
		}
		
		void operator()(blocked_range<int> range) const // computes characteristic function
		{
			for(int iblock=range.begin(); iblock<range.end(); iblock++)
			{
				BlockInfo info = vInfo[iblock];
				ChiBlock& chi_block = *(ChiBlock*)info.ptrBlock;
				
				if(_is_touching(info))
				{
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							chi_block.chi[iy][ix] = this->mollified_heaviside(WingStuff::wing->sdf(p), eps);							
						}
				}
				else
				{
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							chi_block.chi[iy][ix] = 0;
							
							assert( this->mollified_heaviside(WingStuff::wing->sdf(p), eps) < 1e-10 );
						}
				}
			}
		}	
	};
	
	struct UpdateFluidBlocks
	{
		map<int, ChiBlock *>& cached_blocks;
		
		UpdateFluidBlocks(map<int, ChiBlock *>& cached_blocks):cached_blocks(cached_blocks){}
		UpdateFluidBlocks(const UpdateFluidBlocks& c): cached_blocks(c.cached_blocks){}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			map<int, ChiBlock *>::iterator it = cached_blocks.find(info.blockID);
			assert(it != cached_blocks.end());
			assert(it->second != NULL);
			
			ChiBlock& chi_block = *(ChiBlock*)it->second;
			
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					b(ix, iy).tmp = chi_block.chi[iy][ix];
		}
	};
}




I2D_WingObstacleOperator::I2D_WingObstacleOperator(Grid<W,B>& grid, ArgumentParser& parser, Real cm[2], const Real smoothing_length)
:grid(grid), smoothing_length(smoothing_length)
{
	parser.set_strict_mode();
	
	const Real angle = parser("-NACA_angle").asDouble();
	const Real D = parser("-D").asDouble();
	nBleeds = parser("-NACA_nRearBleeds").asInt();
	nTrailingBleeds = parser("-NACA_nTrailingBleeds").asInt();
	plainAirfoil = parser("-NACA_plain").asBool();
	assert(nBleeds>=0 && nBleeds<=8);
	assert(nTrailingBleeds>=0 && nTrailingBleeds<=2);
	WingStuff::wing = new WingStuff::DiscretizedWing(cm[0],cm[1],D,angle,smoothing_length); // get rotated NACA profile
}

void I2D_WingObstacleOperator::getObstacleInfo(vector<Real> & infoObstacle)
{
	infoObstacle.clear();
	
	Real cor[2] = {0.0,0.0};
	WingStuff::wing->getAerodynamicCenter(cor);
	
	infoObstacle.push_back(cor[0]);
	infoObstacle.push_back(cor[1]);	
}

void I2D_WingObstacleOperator::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	// Airfoil
	set<int> valid_blockids;
	vector<BlockInfo> misses;
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		if(cached_blocks.find(it->blockID) == cached_blocks.end())
		{
			BlockInfo info = *it;
			
			ChiBlock * chi = new ChiBlock;
			
			cached_blocks[it->blockID] = chi;
			info.ptrBlock = chi;
			
			misses.push_back(info);
		}
		else
			it->ptrBlock = cached_blocks[it->blockID];
			
		valid_blockids.insert(it->blockID);
	}
	
	vector<int> toerase;
	for(map<int, ChiBlock *>::iterator it=cached_blocks.begin(); it != cached_blocks.end(); it++)
		if(valid_blockids.find(it->first) == valid_blockids.end())
			toerase.push_back(it->first);
	
	for(vector<int>::iterator e = toerase.begin(); e!=toerase.end(); ++e)
	{
		map<int, ChiBlock *>::iterator it = cached_blocks.find(*e);
		
		assert(it!=cached_blocks.end());
		delete it->second;
		it->second = NULL;
		
		cached_blocks.erase(*e);
	}
	
	WingStuff::FillChiBlocks fill;
	fill.eps = smoothing_length;
	fill.vInfo = &misses.front();
	
	tbb::parallel_for(blocked_range<int>(0, misses.size()), fill, auto_partitioner());
	
	WingStuff::UpdateFluidBlocks update(cached_blocks);
	block_processing.process(vInfo, coll, update);

	// Upper indentation
	if(plainAirfoil==false)
	{
		vector<int> bleedToOut;
		bleedToOut.push_back(1);
		bleedToOut.push_back(3);
		bleedToOut.push_back(6);
		for(unsigned int i=0; i<3; i++)
		{
			Real lenghtBleed = 0.0;
			Real widthBleed = 0.0;
			Real angleBleed = 0.0;
			Real cmBleed[2] = {0.0,0.0};	
			Real outletStart[2] = {0.0,0.0};	
			WingStuff::wing->getRearBleedsLocations(bleedToOut[i],lenghtBleed,widthBleed,angleBleed,cmBleed,outletStart);
			WingStuff::FillCylBleed fillCylBleed(widthBleed/2.0+smoothing_length, widthBleed/2.0, smoothing_length, outletStart);
			block_processing.process(vInfo, coll, fillCylBleed);
		}
	}
	
	// Real bleeds
	vector<unsigned int> activatedBleeds;
	
	for(int i=0;i<nBleeds;i++)
		activatedBleeds.push_back(i+1);
	
	for(unsigned int i=0; i<activatedBleeds.size(); i++)
	{
		Real lenghtBleed = 0.0;
		Real widthBleed = 0.0;
		Real angleBleed = 0.0;
		Real cmBleed[2] = {0.0,0.0};	
		Real outletStart[2] = {0.0,0.0};	
		WingStuff::wing->getRearBleedsLocations(activatedBleeds[i],lenghtBleed,widthBleed,angleBleed,cmBleed,outletStart);
		WingStuff::FillBlocksBleed fillBleed(lenghtBleed,widthBleed,smoothing_length,cmBleed,angleBleed);
		block_processing.process(vInfo, coll, fillBleed);
		WingStuff::FillCylBleed fillCylBleedLower(widthBleed/2.0+smoothing_length, widthBleed/2.0, smoothing_length, cmBleed);
		block_processing.process(vInfo, coll, fillCylBleedLower);
		
		if(plainAirfoil==true)
		{
			WingStuff::FillCylBleed fillCylBleedUpper(widthBleed/2.0+smoothing_length, widthBleed/2.0, smoothing_length, outletStart);
			block_processing.process(vInfo, coll, fillCylBleedUpper);
		}
	}
	
	// Trailing bleeds
	vector<unsigned int> activatedTrailingBleeds;
	
	for(int i=0;i<nTrailingBleeds;i++)
		activatedTrailingBleeds.push_back(i+1);
	
	for(unsigned int i=0; i<activatedTrailingBleeds.size(); i++)
	{
		Real lenghtBleed = 0.0;
		Real widthBleed = 0.0;
		Real angleBleed = 0.0;
		Real cmBleed[2] = {0.0,0.0};	
		Real outletStart[2] = {0.0,0.0};	
		WingStuff::wing->getTrailingBleedsLocations(activatedTrailingBleeds[i],lenghtBleed,widthBleed,angleBleed,cmBleed,outletStart);
		WingStuff::FillBlocksBleed fillBleed(lenghtBleed,widthBleed,smoothing_length,cmBleed,angleBleed);
		block_processing.process(vInfo, coll, fillBleed);
		WingStuff::FillCylBleed fillCylBleedLower(widthBleed/2.0+smoothing_length, widthBleed/2.0, smoothing_length, cmBleed);
		block_processing.process(vInfo, coll, fillCylBleedLower);
		WingStuff::FillCylBleed fillCylBleedUpper(widthBleed/2.0+smoothing_length, widthBleed/2.0, smoothing_length, outletStart);
		block_processing.process(vInfo, coll, fillCylBleedUpper);
	}	
}
