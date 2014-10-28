/**
 * @file I2D_ThreeLinkFish.h
 * @date Mar 21, 2012
 * @author Andrew Tchieu
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"

/**
 *
 * Frames of references are are as follows:
 *
 * Frame 0 -> inertial frame, i.e. simulation frame
 * Frame 1 -> attached to ellipse0 center, with ellipse0 horizonal
 * Frame 2 -> attached to center of mass, moving w/ COM, no rotation
 * Frame 3 -> attached to center of mass, moving w/ COM, rotates to conserve angular momentum
 * Frame 4 -> attached to center of mass, moving w/ COM, rotates to conserve angular momentum, rotated by thetav instantaneously
 */
class I2D_ThreeLinkFish: public I2D_FloatingObstacleOperator
{
private:

protected:
	double xm_corr, ym_corr, vx_corr, vy_corr, omega_corr;

public:

	class ThreeEllipse
	{
	private:
		Real _mollified_heaviside(const double dist, const double eps) const;

	protected:
		void _cubicInterpolation(double t0, double t1, double t, double k0, double k1, double & k, double & dkdt);

	public:
		/// Variables that will always be in the inertial frame, Frame0
		double xm, ym, um, vm, thetam, omegam;

		/// Variables that will change frame of reference from Frame1 -> Frame4
		double x0, y0, u0, v0, theta0, omega0;
		double x1, y1, u1, v1, theta1, omega1;
		double x2, y2, u2, v2, theta2, omega2;

		/// Variables to keep track of orientation in Frame4
		double thetav, omegav; /// swimmer orientation WRT its COM (vacuum motion)

		/// User defined parameters
		double D;
		double aspectRatio;
		double separation; /// separation from ellipse to joint, c, (see Kanso)
		double rho; /// density of the object
		double Uinf[2];

		/// User defined motion parameters
		double amplitude1, period1, phase1;
		double amplitude2, period2, phase2;

		/// Dependent Parameters
		double a, b; /// major and minor axis calculated from aspectRatio
		double m0, m; /// mass of a single ellipse, " " whole fish
		double J0, J; /// moment of inertial of a single ellipse, " " whole fish
		double ellipseLength; /// length of each ellipse

		ThreeEllipse();
		ThreeEllipse(Real xm_, Real ym_, Real D, Real angle_rad, Real aspectRatio);

		void computeAndSetShapeChangeWrtFrame1(const double t);
		void computeComCurrent(double& xcom, double& ycom) const;
		void computeComCurrentNumerical(double& xcom, double& ycom) const;
		void computeVelComCurrent(double& ucom, double& vcom) const;
		void computeVelComCurrentNumerical(double& ucom, double& vcom) const;
		double computeAngMomCurrent() const;
		double computeAngMomCurrentNumerical() const;
		void computeAndSetJ();
		void computeAndSetJNumerical();
		void computeAngVelWrtCurrent(double& u, double& v, double omega, double x, double y) const;

		void shiftFromFrame1To2(const double x, const double y, const double u, const double v);
		void shiftFromFrame2To3(const double L);
		void shiftFromFrame3To4();
		void shiftFromFrame4To0();

		void updateShapeAndMotionInVacuum(const double t);
		void updateShapeAndMotionInVacuumNumerical(const double t);
		void updateIntertialQuantities(const double dt);

		void rotateVecWrtCurrent(double& x, double& y, double beta);

		void printExactMomentumDiagnostics(const double t) const;
		void printNumericalMomentumDiagnostics(const double t) const;
		void printShapeLinks(const double t, const string filename) const;
		void printShape(const double t, const string filename) const;

		Real sample(const Real x_, const Real y_, const Real eps) const;
		void sampleDeformationVel(const Real x_, const Real y_, Real& udef, Real& vdef) const;

		///------ OLD ASS SHIT -----------------------------------------------------

		void computeShapeChange(double t);

		void restart(FILE * f);
		void save(FILE * f) const;

		void rotate(Real v[2], Real angle) const;
		void bbox(const Real eps, Real xmin[2], Real xmax[2]) const;
	};

	ThreeEllipse * shape;

	I2D_ThreeLinkFish(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real aspectRatio, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_ThreeLinkFish();

	void characteristic_function();
	Real getD() const {return D;}

	void create(const double t);
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());

	void printGeneralStuff(const double t, const string filename, vector<double> stuffToPrint);

	Real getModulusMaxVel() {return max(shape->D/shape->period1, shape->D/shape->period2);}
};
