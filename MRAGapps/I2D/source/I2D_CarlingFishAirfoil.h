/*
 * I2D_CarlingFishAirfoil.h
 *
 *  Created on: Feb 16, 2012
 *      Author: mgazzola
 */

#ifndef I2D_CARLINGFISHAIRFOIL_H_
#define I2D_CARLINGFISHAIRFOIL_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_CarlingFish.h"

class I2D_CarlingFishAirfoil : public I2D_CarlingFish
{
public:

	class CarlingFishAirfoil: public Fish
		{
		protected:
			double _getWidth(const double & ss) const;

		public:
			string naca;

			CarlingFishAirfoil(string _naca, double xm, double ym, double _D, double _T, double phase, double angle_rad, double angleInSpace_rad, double eps, const int LMAX);
			virtual ~CarlingFishAirfoil(){};
		};

	I2D_CarlingFishAirfoil(ArgumentParser & parser, Grid<W,B>& grid, string _naca, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real angle, const Real angleInSpace, const Real eps, const Real Uinf[2],
			I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** policy = NULL);
	virtual ~I2D_CarlingFishAirfoil();
};

#endif /* I2D_CARLINGFISHAIRFOIL_H_ */
