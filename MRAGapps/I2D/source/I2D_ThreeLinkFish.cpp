/**
 * @file IF2DThreeLinkFish.cpp
 * @date Mar 21, 2012
 * @author Andrew Tchieu
 *
 * @todo Fix mass for all other ObjectFactory objects that include ThreeEllipses.
 */

#include "I2D_ThreeLinkFish.h"
#include <iostream>
#include <fstream>
/**
 * Default constructor for ThreeEllipse structure
 */
I2D_ThreeLinkFish::ThreeEllipse::ThreeEllipse(): thetam(0), D(0.1), xm(0.5), ym(0.5), aspectRatio(2), omegam(0), um(0), vm(0), rho(1.0)
{
	/// @todo Remove hard coded stuff
	thetam = 0.5*M_PI;
	omegam = 0.0;

	x0 = y0 = u0 = v0 = theta0 = omega0 = 0.0;
	x1 = y1 = u1 = v1 = theta1 = omega1 = 0.0;
	x2 = y2 = u2 = v2 = theta2 = omega2 = 0.0;
	thetav = omegav = 0.0;

	/// @todo Remove hard coded parameters
	separation = 0.025*D;
	rho = 1.0;

	/// Motion
	period1 = 1.0;
	period2 = 1.0;
	amplitude1 = 0.5*M_PI;
	amplitude2 = 0.5*M_PI;
	phase1 = 0.0;
	phase2 = 0.0;

	/// Dependent parameters
	ellipseLength = (D-4.0*separation)/3.0;
	a = 0.5*ellipseLength;
	b = aspectRatio*a;

	m0 = rho*M_PI*a*b;
	m = 3.0*m0;
	J0 = 0.25*M_PI*a*b*(a*a+b*b);

	updateShapeAndMotionInVacuum(0.0);

	const double dist1Sqr = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
	const double dist2Sqr = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0);
	J = 3*J0 + m0*(dist1Sqr + dist2Sqr); /// Use parallel axis theorem, will be updated, don't worry about details
}

/**
 * Overloaded constructor for ThreeEllipse structure
 */
I2D_ThreeLinkFish::ThreeEllipse::ThreeEllipse(Real xm_, Real ym_, Real D, Real angle_rad, Real aspectRatio): thetam(angle_rad), D(D), xm(xm_), ym(ym_), aspectRatio(aspectRatio), omegam(0), um(0), vm(0), rho(1.0)
{
	/// @todo Remove hard coded stuff
	thetam = 0.5*M_PI;
	omegam = 0.0;

	x0 = y0 = u0 = v0 = theta0 = omega0 = 0.0;
	x1 = y1 = u1 = v1 = theta1 = omega1 = 0.0;
	x2 = y2 = u2 = v2 = theta2 = omega2 = 0.0;
	thetav = omegav = 0.0;

	/// @todo Remove hard coded parameters
	separation = 0.025*D;
	rho = 1.0;

	/// Motion
	period1 = 1.0;
	period2 = 1.0;
	amplitude1 = 0.5*M_PI;
	amplitude2 = 0.5*M_PI;
	phase1 = 0.0;
	phase2 = 0.0;

	/// Dependent parameters
	ellipseLength = (D-4.0*separation)/3.0;
	a = 0.5*ellipseLength;
	b = aspectRatio*a;

	m0 = rho*M_PI*a*b;
	m = 3.0*m0;
	J0 = 0.25*M_PI*a*b*(a*a+b*b);

	updateShapeAndMotionInVacuum(0.0);

	const double dist1Sqr = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
	const double dist2Sqr = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0);
	J = 3*J0 + m0*(dist1Sqr + dist2Sqr); /// Use parallel axis theorem, will be updated, don't worry about details
}

/**
 * Defines the smoothed heavyside function given the distance from the edge
 * and its smoothing length. Returns 1 inside (negative distance) and 0
 * outside (positive distance).
 *
 * @param dist
 * @param eps
 * @return
 */
Real I2D_ThreeLinkFish::ThreeEllipse::_mollified_heaviside(const double dist, const double eps) const
{
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));
	return 0.5+0.5*cos(alpha);
}

/**
 * Interpolates at time t between points [t0, t1], given the values
 * [k0(t0), k1(t1)]. Returns the value and its derivative.
 *
 * @param t0
 * @param t1
 * @param t
 * @param k0
 * @param k1
 * @param k
 * @param dkdt
 */
void I2D_ThreeLinkFish::ThreeEllipse::_cubicInterpolation(double t0, double t1, double t, double k0, double k1, double & k, double & dkdt)
{
	k = 0.0;
	dkdt = 0.0;

	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;

	const double dT = t1-t0;
	const double dt = t-t0;

	// CUBIC
	d = k0;
	c = 0.0;
	b = 3.0*(k1-k0)/(dT*dT);
	a = -2.0*(k1-k0)/(dT*dT*dT);
	k = a*dt*dt*dt + b*dt*dt + c*dt + d;
	dkdt = 3.0*a*dt*dt + 2.0*b*dt + c;
}


/**
 * Scan in restart conditions.
 *
 * CHANGE!!!
 *
 * @param f
 */
void I2D_ThreeLinkFish::ThreeEllipse::restart(FILE * f)
{
	float val;

	fscanf(f, "xm: %e\n", &val);
	xm = val;
	printf("ThreeEllipse::restart(): xm is %e\n", xm);

	fscanf(f, "ym: %e\n", &val);
	ym = val;
	printf("ThreeEllipse::restart(): ym is %e\n", ym);

	fscanf(f, "um: %e\n", &val);
	um = val;
	printf("ThreeEllipse::restart(): um is %e\n", um);

	fscanf(f, "vm: %e\n", &val);
	vm = val;
	printf("ThreeEllipse::restart(): vm is %e\n", vm);

	fscanf(f, "omegam: %e\n", &val);
	omegam = val;
	printf("ThreeEllipse::restart(): omegam is %e\n", omegam);

	fscanf(f, "thetam: %e\n", &val);
	thetam = val;
	printf("ThreeEllipse::restart(): thetam is %e\n", thetam);
}

/**
 * Save state of object to file for restarts.
 *
 * CHANGE!!!
 *
 * @param f
 */
void I2D_ThreeLinkFish::ThreeEllipse::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "um: %20.20e\n", um);
	fprintf(f, "vm: %20.20e\n", vm);
	fprintf(f, "omegam: %20.20e\n", omegam);
	fprintf(f, "thetam: %20.20e\n", thetam);
}

/**
 * Calculates distance from boundary and then returns back the molified
 * heavyside function at location x_, y_.
 *
 * CHANGE!!!
 *
 * @param x_
 * @param y_
 * @param eps
 * @return
 *
 * @todo Probably a better way to do this than sampling three times
 */
Real I2D_ThreeLinkFish::ThreeEllipse::sample(const Real x_, const Real y_, const Real eps) const
{
	/// ThreeEllipse0
	const double alpha0 = atan2(y_-y0,x_-x0) - theta0;
	const double radius0 = a*b/sqrt(b*b*cos(alpha0)*cos(alpha0) + a*a*sin(alpha0)*sin(alpha0));
	const double dist0 = sqrt((x_-x0)*(x_-x0) + (y_-y0)*(y_-y0)) - radius0;
	const Real ThreeEllipse0 = _mollified_heaviside(dist0, eps);
	Real myMax = ThreeEllipse0;

	const double alpha1 = atan2(y_-y1,x_-x1) - theta1;
	const double radius1 = a*b/sqrt(b*b*cos(alpha1)*cos(alpha1) + a*a*sin(alpha1)*sin(alpha1));
	const double dist1 = sqrt((x_-x1)*(x_-x1) + (y_-y1)*(y_-y1)) - radius1;
	const Real ThreeEllipse1 = _mollified_heaviside(dist1, eps);
	myMax = max(myMax, ThreeEllipse1);

	const double alpha2 = atan2(y_-y2,x_-x2) - theta2;
	const double radius2 = a*b/sqrt(b*b*cos(alpha2)*cos(alpha2) + a*a*sin(alpha2)*sin(alpha2));
	const double dist2 = sqrt((x_-x2)*(x_-x2) + (y_-y2)*(y_-y2)) - radius2;
	const Real ThreeEllipse2 = _mollified_heaviside(dist2, eps);
	myMax = max(myMax, ThreeEllipse2);

	return myMax;
}

void I2D_ThreeLinkFish::ThreeEllipse::sampleDeformationVel(const Real x_, const Real y_, Real& udef, Real& vdef) const
{
	const Real dist0Sqr = (x_-x0)*(x_-x0) + (y_-y0)*(y_-y0);
	const Real dist1Sqr = (x_-x1)*(x_-x1) + (y_-y1)*(y_-y1);
	const Real dist2Sqr = (x_-x2)*(x_-x2) + (y_-y2)*(y_-y2);
	if ((dist0Sqr < dist1Sqr) && (dist0Sqr < dist2Sqr)) /// ellipse 0
	{
		const Real urot = -omega0*(y_-y0);
		const Real vrot = +omega0*(x_-x0);
		udef = u0 + urot;
		vdef = v0 + vrot;
	}
	else if ((dist1Sqr < dist0Sqr) && (dist1Sqr < dist2Sqr)) /// ellipse 1
	{
		const Real urot = -omega1*(y_-y1);
		const Real vrot = +omega1*(x_-x1);
		udef = u1 + urot;
		vdef = v1 + vrot;
	}
	else if ((dist2Sqr < dist0Sqr) && (dist2Sqr < dist1Sqr))/// ellipse 2
	{
		const Real urot = -omega2*(y_-y2);
		const Real vrot = +omega2*(x_-x2);
		udef = u2 + urot;
		vdef = v2 + vrot;
	}
	else
	{
		printf("WHAT THE FUCK ARE YOU DOING HERE?!? YOU SHOULD QUIT YOUR JOB.\n");
		exit(0);
	}
}

/**
 * Rotates a given two dimensional vector by alpha about the ThreeEllipses COM.
 *
 * @param[in,out] v
 * @param angle
 */
void I2D_ThreeLinkFish::ThreeEllipse::rotate(Real v[2], Real angle) const
{
	const Real a00 = cos(angle);
	const Real a01 = -sin(angle);
	const Real a10 = sin(angle);
	const Real a11 = cos(angle);

	const Real xv = v[0];
	const Real yv = v[1];

	v[0] = a00*(xv-xm) + a01*(yv-ym);
	v[1] = a10*(xv-xm) + a11*(yv-ym);

	v[0] += xm;
	v[1] += ym;
}

void I2D_ThreeLinkFish::ThreeEllipse::computeAndSetShapeChangeWrtFrame1(const double t)
{
	assert(separation > 0.0);
	assert(a > 0.0);
	const double apc = a+separation;
	const double freq1 = 2.0*M_PI/period1;
	const double freq2 = 2.0*M_PI/period2;

	/// Update angles
	theta0 = 0.0;
	omega0 = 0.0;
	theta1 = amplitude1*cos(freq1*t+phase1);
	omega1 = -amplitude1*freq1*sin(freq1*t+phase1);
	theta2 = amplitude2*cos(freq2*t+phase2);
	omega2 = -amplitude2*freq2*sin(freq2*t+phase2);

	/// Update positions and velocities
	x0 = 0.0;
	y0 = 0.0;
	u0 = 0.0;
	v0 = 0.0;
	x1 = apc*(1+cos(theta1));
	y1 = apc*sin(theta1);
	u1 = -omega1*y1;
	v1 = omega1*(x1-apc);
	x2 = -apc*(1+cos(theta2));
	y2 = -apc*sin(theta2);
	u2 = -omega2*y2;
	v2 = omega2*(x2+apc);
}

void I2D_ThreeLinkFish::ThreeEllipse::computeComCurrent(double& xcom, double& ycom) const
{
	assert(m0 != 0.0);
	assert(m != 0.0);
	xcom = m0*(x0+x1+x2)/m;
	ycom = m0*(y0+y1+y2)/m;
}

void I2D_ThreeLinkFish::ThreeEllipse::computeComCurrentNumerical(double& xcom, double& ycom) const
{
	assert(m0 != 0.0);
	assert(m != 0.0);
	xcom = m0*(x0+x1+x2)/m;
	ycom = m0*(y0+y1+y2)/m;
}

void I2D_ThreeLinkFish::ThreeEllipse::computeVelComCurrent(double& ucom, double& vcom) const
{
	assert(m0 != 0.0);
	assert(m != 0.0);
	ucom = m0*(u0+u1+u2)/m;
	vcom = m0*(v0+v1+v2)/m;
}

void I2D_ThreeLinkFish::ThreeEllipse::computeVelComCurrentNumerical(double& ucom, double& vcom) const
{
	assert(m0 != 0.0);
	assert(m != 0.0);
	ucom = m0*(u0+u1+u2)/m;
	vcom = m0*(v0+v1+v2)/m;
}

double I2D_ThreeLinkFish::ThreeEllipse::computeAngMomCurrent() const
{
	assert(m0 != 0.0);
	assert(J0 != 0.0);
	const double L0 = m0*(x0*v0-y0*u0);
	const double L1 = m0*(x1*v1-y1*u1);
	const double L2 = m0*(x2*v2-y2*u2);
	const double L0J = J0*omega0;
	const double L1J = J0*omega1;
	const double L2J = J0*omega2;
	return (L0+L1+L2+L0J+L1J+L2J);
}

double I2D_ThreeLinkFish::ThreeEllipse::computeAngMomCurrentNumerical() const
{
	assert(m0 != 0.0);
	assert(J0 != 0.0);
	const double L0 = m0*(x0*v0-y0*u0);
	const double L1 = m0*(x1*v1-y1*u1);
	const double L2 = m0*(x2*v2-y2*u2);
	const double L0J = J0*omega0;
	const double L1J = J0*omega1;
	const double L2J = J0*omega2;
	return (L0+L1+L2+L0J+L1J+L2J);
}


void I2D_ThreeLinkFish::ThreeEllipse::computeAndSetJ()
{
	assert(J0 != 0.0);
	const double dist0Sqr = x0*x0+y0*y0;
	const double dist1Sqr = x1*x1+y1*y1;
	const double dist2Sqr = x2*x2+y2*y2;
	J = 3*J0+m0*(dist0Sqr+dist1Sqr+dist2Sqr);
}

void I2D_ThreeLinkFish::ThreeEllipse::computeAndSetJNumerical()
{
	assert(J0 != 0.0);
	const double dist0Sqr = x0*x0+y0*y0;
	const double dist1Sqr = x1*x1+y1*y1;
	const double dist2Sqr = x2*x2+y2*y2;
	J = 3*J0+m0*(dist0Sqr+dist1Sqr+dist2Sqr);
}

void I2D_ThreeLinkFish::ThreeEllipse::shiftFromFrame1To2(const double x, const double y, const double u, const double v)
{
	x0 -= x;
	y0 -= y;
	u0 -= u;
	v0 -= v;
	x1 -= x;
	y1 -= y;
	u1 -= u;
	v1 -= v;
	x2 -= x;
	y2 -= y;
	u2 -= u;
	v2 -= v;
}

void I2D_ThreeLinkFish::ThreeEllipse::shiftFromFrame2To3(const double L)
{
	assert(J > 0);
	omegav = -L/J;
	omega0 += omegav;
	omega1 += omegav;
	omega2 += omegav;

	double uCorrect, vCorrect;
	computeAngVelWrtCurrent(uCorrect, vCorrect, omegav, x0, y0);
	u0 += uCorrect;
	v0 += vCorrect;
	computeAngVelWrtCurrent(uCorrect, vCorrect, omegav, x1, y1);
	u1 += uCorrect;
	v1 += vCorrect;
	computeAngVelWrtCurrent(uCorrect, vCorrect, omegav, x2, y2);
	u2 += uCorrect;
	v2 += vCorrect;
}

void I2D_ThreeLinkFish::ThreeEllipse::shiftFromFrame3To4()
{
	theta0 += thetav;
	theta1 += thetav;
	theta2 += thetav;
	rotateVecWrtCurrent(x0, y0, thetav);
	rotateVecWrtCurrent(u0, v0, thetav);
	rotateVecWrtCurrent(x1, y1, thetav);
	rotateVecWrtCurrent(u1, v2, thetav);
	rotateVecWrtCurrent(x2, y2, thetav);
	rotateVecWrtCurrent(u2, v2, thetav);
}

void I2D_ThreeLinkFish::ThreeEllipse::shiftFromFrame4To0()
{
	theta0 += thetam;
	theta1 += thetam;
	theta2 += thetam;
	rotateVecWrtCurrent(x0, y0, thetam);
	rotateVecWrtCurrent(u0, v0, thetam);
	rotateVecWrtCurrent(x1, y1, thetam);
	rotateVecWrtCurrent(u1, v2, thetam);
	rotateVecWrtCurrent(x2, y2, thetam);
	rotateVecWrtCurrent(u2, v2, thetam);
	x0 += xm;
	y0 += ym;
	x1 += xm;
	y1 += ym;
	x2 += xm;
	y2 += ym;
	assert(x0 > 0 && x0 < 1);
	assert(y0 > 0 && y0 < 1);
	assert(x1 > 0 && x1 < 1);
	assert(y1 > 0 && y1 < 1);
	assert(x2 > 0 && x2 < 1);
	assert(y2 > 0 && y2 < 1);
}



void I2D_ThreeLinkFish::ThreeEllipse::updateShapeAndMotionInVacuum(const double t) {

	double xcom, ycom, ucom, vcom;
	double extraAngMom; /// amount of extra angular momentum

	computeAndSetShapeChangeWrtFrame1(t);
	printShapeLinks(t, "ThreeLinkFrame1.txt");
	computeComCurrent(xcom, ycom);
	computeVelComCurrent(ucom, vcom);
	shiftFromFrame1To2(xcom, ycom, ucom, vcom);
	printShapeLinks(t, "ThreeLinkFrame2.txt");
	extraAngMom = computeAngMomCurrent();
	computeAndSetJ();
	shiftFromFrame2To3(extraAngMom);
	printShapeLinks(t, "ThreeLinkFrame2.txt");
	shiftFromFrame3To4();
	printShapeLinks(t, "ThreeLinkFrame3.txt");
	shiftFromFrame4To0();
	printShapeLinks(t, "ThreeLinkFrame4.txt");
	printExactMomentumDiagnostics(t);
}

void I2D_ThreeLinkFish::ThreeEllipse::updateShapeAndMotionInVacuumNumerical(const double t)
{
	double xcom, ycom, ucom, vcom;
	double extraAngMom; /// amount of extra angular momentum

	computeAndSetShapeChangeWrtFrame1(t);
	printShapeLinks(t, "ThreeLinkFrame1.txt");
	computeComCurrent(xcom, ycom); /// move it to where you expect using exact results, for an estimate
	shiftFromFrame1To2(xcom, ycom, 0.0, 0.0);
	computeComCurrentNumerical(xcom, ycom); /// move it a little bit if necessary
	shiftFromFrame1To2(xcom, ycom, 0.0, 0.0);
	computeVelComCurrentNumerical(ucom, vcom);
	shiftFromFrame1To2(0.0, 0.0, ucom, vcom);
	printShapeLinks(t, "ThreeLinkFrame2.txt");
	extraAngMom = computeAngMomCurrentNumerical();
	computeAndSetJNumerical();
	shiftFromFrame2To3(extraAngMom);
	printShapeLinks(t, "ThreeLinkFrame2.txt");
	shiftFromFrame3To4();
	printShapeLinks(t, "ThreeLinkFrame3.txt");
	shiftFromFrame4To0();
	printShapeLinks(t, "ThreeLinkFrame4.txt");
	printExactMomentumDiagnostics(t);
}


void I2D_ThreeLinkFish::ThreeEllipse::updateIntertialQuantities(const double dt)
{
	xm += um*dt;
	ym += vm*dt;
	thetam += omegam*dt;
	thetav += omegav*dt; // Remember to update the internal angle
}

void I2D_ThreeLinkFish::ThreeEllipse::computeAngVelWrtCurrent(double& u, double& v, double omega, double x, double y) const
{
	u = -omega*y;
	v = omega*x;
}

void I2D_ThreeLinkFish::ThreeEllipse::rotateVecWrtCurrent(double& x, double& y, double beta)
{
	const double a00 = cos(beta);
	const double a01 = -sin(beta);
	const double a10 = sin(beta);
	const double a11 = cos(beta);
	const double xc = x;
	const double yc = y;

	x = a00*(xc) + a01*(yc);
	y = a10*(xc) + a11*(yc);
}

void I2D_ThreeLinkFish::ThreeEllipse::printExactMomentumDiagnostics(const double t) const
{
	double L;
	double xcom, ycom;
	double ucom, vcom;
	computeComCurrent(xcom, ycom);
	computeVelComCurrent(ucom, vcom);
	L = computeAngMomCurrent();

	FILE * outfile = NULL;
	outfile = fopen("ThreeLinkExactDiagnostics.txt", t == 0.0 ? "w" : "a");
	assert(outfile != NULL);
	fprintf(outfile, "%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n", t, xcom, ycom, ucom, vcom, L);
	fclose(outfile);
}

void I2D_ThreeLinkFish::ThreeEllipse::printNumericalMomentumDiagnostics(const double t) const
{
	double L;
	double xcom, ycom;
	double ucom, vcom;
	computeComCurrent(xcom, ycom);
	computeVelComCurrent(ucom, vcom);
	L = computeAngMomCurrent();

	FILE * outfile = NULL;
	outfile = fopen("ThreeLinkExactDiagnostics.txt", t == 0.0 ? "w" : "a");
	assert(outfile != NULL);
	fprintf(outfile, "%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n", t, xcom, ycom, ucom, vcom, L);
	fclose(outfile);
}

void I2D_ThreeLinkFish::ThreeEllipse::printShapeLinks(const double t, const string filename) const
{
	FILE * outfile = NULL;
	outfile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
	assert(outfile != NULL);

	fprintf(outfile, "%12.5e ", t);
	fprintf(outfile, "%12.5e %12.5e %12.5e " , x0, y0, theta0);
	fprintf(outfile, "%12.5e %12.5e %12.5e " , x1, y1, theta1);
	fprintf(outfile, "%12.5e %12.5e %12.5e\n", x2, y2, theta2);
	fclose(outfile);
}


void I2D_ThreeLinkFish::ThreeEllipse::printShape(const double t, const string filename) const
{
	FILE* outfile = NULL;

	if(filename==std::string())
		outfile = fopen("update_I2D_ThreeLinkFish.txt", t == 0.0 ? "w" : "a");
	else
		outfile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");

	assert(outfile != NULL);

	fprintf(outfile, "%12.5e ", t);
	fprintf(outfile, "%12.5e %12.5e " , xm, ym);
	fprintf(outfile, "%12.5e %12.5e " , um, vm);
	fprintf(outfile, "%12.5e %12.5e " , thetam, omegam);
	fprintf(outfile, "%12.5e %12.5e " , thetav, omegav);
	fprintf(outfile, "%12.5e %12.5e\n" , m, J);
	fclose(outfile);
}


/**
 * Find the bounding box around the ThreeEllipse.
 *
 * CHANGE!!!
 *
 * @param eps
 * @param xmin
 * @param xmax
 * @todo Make this a tighter box
 */
void I2D_ThreeLinkFish::ThreeEllipse::bbox(const Real eps, Real xmin[2], Real xmax[2]) const
{
	xmin[0] = x0 - D/2.0; /// length of entire fish is D
	xmax[0] = x0 + D/2.0;
	xmin[1] = y0 - (separation+ellipseLength); /// maximum of one ThreeEllipse length
	xmax[1] = y0 + (separation+ellipseLength);

	Real v1[2] = { xmin[0], xmin[1] }; /// bottom left corner
	Real v2[2] = { xmax[0], xmax[1] }; /// upper right corner
	Real v3[2] = { xmin[0], xmax[1] }; /// upper left corner
	Real v4[2] = { xmax[0], xmin[1] }; /// bottom right corner

	rotate(v1, thetam); /// rotate the box by body orientation
	rotate(v2, thetam);
	rotate(v3, thetam);
	rotate(v4, thetam);

	xmin[0] = min((Real)min((Real)min(v1[0], v2[0]), v3[0]), v4[0]); /// min x
	xmax[0] = max((Real)max((Real)max(v1[0], v2[0]), v3[0]), v4[0]); /// max x
	xmin[1] = min((Real)min((Real)min(v1[1], v2[1]), v3[1]), v4[1]); /// min y
	xmax[1] = max((Real)max((Real)max(v1[1], v2[1]), v3[1]), v4[1]); /// max y

	assert(eps >= 0);

	xmin[0] -= 2*eps; /// add smoothing length
	xmin[1] -= 2*eps;
	xmax[0] += 2*eps;
	xmax[1] += 2*eps;

	assert(xmin[0] <= xmax[0]);
	assert(xmin[1] <= xmax[1]);
}

/**
 * Create a special namespace to avoid confusion when using MRAG and overloaded
 * operations.
 */
namespace ThreeLinkFish
{
  /**
  * Used when defining the characteristic function. Helpers are used to see
  * whether the block is touching any part where the characteristic function
  * is nonzero.
  *
  * @struct FillBlocks
  */
	struct FillBlocks
	{
		Real eps;
		I2D_ThreeLinkFish::ThreeEllipse * shape;

		/**
		 * Default constructor
		 * @param eps
		 * @param _shape
		 */
		FillBlocks(Real eps, I2D_ThreeLinkFish::ThreeEllipse *_shape): eps(eps) {shape = _shape;}

		/**
		 * Copy constructor
		 * @param c
		 * @return
		 */
		FillBlocks(const FillBlocks& c): eps(c.eps) {shape = c.shape;}

		/**
		 * Checks to see if the block is touching any part of the body or not. This
		 * is used to figure whether computing things in this block is necessary.
		 *
		 * @param eps
		 * @param wheel
		 * @param info
		 * @return
		 */
		static bool _is_touching(Real eps, const I2D_ThreeLinkFish::ThreeEllipse * wheel, const BlockInfo& info)
		{
			Real min_pos[2], max_pos[2];

			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);

			Real bbox[2][2];
			wheel->bbox(eps, bbox[0], bbox[1]);

			Real intersection[2][2] = {
				max(min_pos[0], bbox[0][0]), min(max_pos[0], bbox[1][0]),
				max(min_pos[1], bbox[0][1]), min(max_pos[1], bbox[1][1]),
			};

			return
			intersection[0][1]-intersection[0][0]>0 &&
			intersection[1][1]-intersection[1][0]>0 ;
		}

		/**
		 * Overload used in block_processing.process(vInfo, coll, fill). This
		 * specific version gets the position and then calls sample, which in turn
		 * calculates the characteristic function. Then it "returns" the processed
		 * block.
		 *
		 * @param info
		 * @param b
		 */
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			if(_is_touching(eps, shape, info))
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				{
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);
						b(ix,iy).tmp = max( shape->sample(p[0], p[1], eps), b(ix,iy).tmp );
					}
				}
			}
		}
	};

	struct ComputeCenterOfMass
	{
		Real Uinf[2];
		double eps;
		I2D_ThreeLinkFish::ThreeEllipse * shape;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;

		ComputeCenterOfMass(double eps, I2D_ThreeLinkFish::ThreeEllipse *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
		{
			this->Uinf[0] = Uinf[0];
			this->Uinf[1] = Uinf[1];

			shape = _shape;
		}

		ComputeCenterOfMass(const ComputeCenterOfMass& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
		{
			Uinf[0] = c.Uinf[0];
			Uinf[1] = c.Uinf[1];

			shape = c.shape;
		}

		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			bool bNonEmpty = false;

			if(FillBlocks::_is_touching(eps, shape, info))
			{
				double mass = 0;
				double xbar = 0;
				double ybar = 0;

				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						double p[2];
						info.pos(p, ix, iy);

						const double Xs = shape->sample(p[0], p[1], eps);
						bNonEmpty |= Xs > 0;

						mass += Xs;
						xbar += Xs*p[0];
						ybar += Xs*p[1];
					}

				assert(b2sum.find(info.blockID) != b2sum.end());
				assert(b2nonempty.find(info.blockID) != b2nonempty.end());

				b2sum[info.blockID][0] = mass*info.h[0]*info.h[0];
				b2sum[info.blockID][1] = xbar*info.h[0]*info.h[0];
				b2sum[info.blockID][2] = ybar*info.h[0]*info.h[0];

				b2nonempty[info.blockID] = bNonEmpty;
			}
		}
	};

	/**
	 *
	 */
	struct ComputeAllCorrections
	{
		Real Uinf[2];
		double eps;
		double xcm, ycm;
		I2D_ThreeLinkFish::ThreeEllipse* shape;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;

		ComputeAllCorrections(double xcm, double ycm, double eps, I2D_ThreeLinkFish::ThreeEllipse* _shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): xcm(xcm), ycm(ycm), eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
		{
			this->Uinf[0] = Uinf[0];
			this->Uinf[1] = Uinf[1];

			shape = _shape;
		}

		ComputeAllCorrections(const ComputeAllCorrections& c): xcm(c.xcm), ycm(c.ycm), eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
		{
			Uinf[0] = c.Uinf[0];
			Uinf[1] = c.Uinf[1];

			shape = c.shape;
		}

		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			bool bNonEmpty = false;

			if(FillBlocks::_is_touching(eps, shape, info))
			{
				double mass = 0.0;
				double J = 0.0;
				double vxbar = 0.0;
				double vybar = 0.0;
				double omegabar = 0.0;

				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						double p[2];
						info.pos(p, ix, iy);

						Real Xs = 0.0;
						Real vdefx = 0.0;
						Real vdefy = 0.0;
						Xs = shape->sample(p[0], p[1], eps);
						if (Xs > 0)
							shape->sampleDeformationVel(p[0], p[1], vdefx, vdefy);
						else
						{
							vdefx = 0.0;
							vdefy = 0.0;
						}
						bNonEmpty |= Xs>0;

						mass += Xs;
						J += Xs*((p[0]-xcm)*(p[0]-xcm) + (p[1]-ycm)*(p[1]-ycm));
						vxbar += Xs*vdefx;
						vybar += Xs*vdefy;
						omegabar += Xs*(vdefy*(p[0]-xcm)-vdefx*(p[1]-ycm));
					}

				assert(b2sum.find(info.blockID) != b2sum.end());
				assert(b2nonempty.find(info.blockID) != b2nonempty.end());

				b2sum[info.blockID][0] = mass*info.h[0]*info.h[0];
				b2sum[info.blockID][1] = J*info.h[0]*info.h[0];
				b2sum[info.blockID][2] = vxbar*info.h[0]*info.h[0];
				b2sum[info.blockID][3] = vybar*info.h[0]*info.h[0];
				b2sum[info.blockID][4] = omegabar*info.h[0]*info.h[0];

				b2nonempty[info.blockID] = bNonEmpty;
			}
		}
	};

	/**
	 * Structure for computing all integrals, i.e. for mass, moment of inertia,
	 * Eq. (11), and Eq. (12) in the PoF paper.
	 *
	 * @struct ComputeAll
	 */
	struct ComputeAll
	{
		Real Uinf[2];
		double eps;
		I2D_ThreeLinkFish::ThreeEllipse * shape;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;

		/**
		 * Constructor
		 *
		 * @param eps
		 * @param _shape
		 * @param b2sum
		 * @param Uinf
		 * @param b2nonempty
		 */
		ComputeAll(double eps, I2D_ThreeLinkFish::ThreeEllipse *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
		{
			this->Uinf[0] = Uinf[0];
			this->Uinf[1] = Uinf[1];
			shape = _shape;
		}

		/**
		 * Copy constructor
		 *
		 * @param c
		 */
		ComputeAll(const ComputeAll& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
		{
			Uinf[0] = c.Uinf[0];
			Uinf[1] = c.Uinf[1];

			shape = c.shape;
		}

		/**
		 * Overload for use in block_processing.process(vInfo, coll, computeAll).
		 * This function calculates the integrals to determine the mass, moment of
		 * inertia, velocity, and angular velocity of the body as given in Eq. (11)
		 * and (12) in the PoF paper.
		 *
		 * @param info
		 * @param b
		 */
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{
			bool bNonEmpty = false;

			if(FillBlocks::_is_touching(eps, shape, info))
			{
				double mass = 0;
				double J = 0;
				double vxbar = 0;
				double vybar = 0;
				double omegabar = 0;

				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						double p[2];
						info.pos(p, ix, iy);

						const double Xs = shape->sample(p[0], p[1], eps);
						bNonEmpty |= Xs>0;

						mass += Xs;
						J += Xs*((p[0]-shape->xm)*(p[0]-shape->xm) + (p[1]-shape->ym)*(p[1]-shape->ym));
						vxbar += Xs*(b(ix, iy).u[0]+Uinf[0]);
						vybar += Xs*(b(ix, iy).u[1]+Uinf[1]);
						omegabar += Xs*((b(ix, iy).u[1]+Uinf[1])*(p[0]-shape->xm)-(b(ix, iy).u[0]+Uinf[0])*(p[1]-shape->ym));
					}

				assert(b2sum.find(info.blockID) != b2sum.end());
				assert(b2nonempty.find(info.blockID) != b2nonempty.end());

				b2sum[info.blockID][0] = mass*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][1] = J*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][2] = vxbar*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][3] = vybar*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][4] = omegabar*info.h[0]*info.h[0]*shape->rho;
				b2nonempty[info.blockID] = bNonEmpty;
			}
		}
	};

	/**
	 * Fill the velocity blocks with the velocity and angular velocities
	 * calculated from ComputeAll.
	 *
	 * @struct FillVelblocks
	 */
	struct FillVelblocks
	{
		double eps, xcm, ycm, vxcorr, vycorr, omegacorr;
		I2D_ThreeLinkFish::ThreeEllipse * shape;
		vector<pair< BlockInfo, VelocityBlock *> >& workitems;

		/**
		 * Constructor
		 *
		 * @param workitems
		 * @param eps
		 * @param _shape
		 */
		FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double xcm, double ycm, double vxcorr, double vycorr, double omegacorr, double eps, I2D_ThreeLinkFish::ThreeEllipse *_shape)
		: workitems(workitems), xcm(xcm), ycm(ycm), vxcorr(vxcorr), vycorr(vycorr), omegacorr(omegacorr), eps(eps)
		{
			shape = _shape;
		}

		/**
		 * Copy constructor
		 *
		 * @param c
		 */
		FillVelblocks(const FillVelblocks& c)
		: workitems(c.workitems), xcm(c.xcm), ycm(c.ycm), vxcorr(c.vxcorr), vycorr(c.vycorr), omegacorr(c.omegacorr), eps(c.eps)
		{
			shape = c.shape;
		}

		/**
		 * Overload for use in tbb::parallel_for(...); Sets the velocity according
		 * to Eq. (13), specifically (mind the angular rotation).
		 *
		 * @param range
		 */
		inline void operator()(blocked_range<int> range) const
		{
			for(int iblock=range.begin(); iblock<range.end(); iblock++)
			{
				BlockInfo info = workitems[iblock].first;
				VelocityBlock * u_desired = workitems[iblock].second;

				const Real xm = xcm;
				const Real ym = ycm;
				const Real omegam = shape->omegam - omegacorr;
				const Real um = shape->um - vxcorr;
				const Real vm = shape->vm - vycorr;

				if(FillBlocks::_is_touching(eps, shape, info))
				{
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					{
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);

							const Real Xs = shape->sample(p[0], p[1], eps);

							if (Xs > 0)
							{
								Real udef = 0.0;
								Real vdef = 0.0;
								shape->sampleDeformationVel(p[0], p[1], udef, vdef);

								u_desired->u[0][iy][ix] = - omegam*(p[1]-ym) + um + udef;
								u_desired->u[1][iy][ix] = + omegam*(p[0]-xm) + vm + vdef;
							}
							else
							{
								u_desired->u[0][iy][ix] = 0.0;
								u_desired->u[1][iy][ix] = 0.0;
							}
						}
					}
				}
			}
		}
	};
}


/**
 * Constructor
 *
 * @param grid
 * @param _xm
 * @param _ym
 * @param D
 * @param aspectRatio
 * @param eps
 * @param Uinf
 * @param penalization
 */
I2D_ThreeLinkFish::I2D_ThreeLinkFish(ArgumentParser & parser, Grid<W,B>& grid,
		const Real _x0, const Real _y0, const Real D, const Real aspectRatio,
		const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization)
: I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization),
  xm_corr(0.0), ym_corr(0.0), vx_corr(0.0), vy_corr(0.0), omega_corr(0.0)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
	shape = new ThreeEllipse(_x0, _y0, D, 0.0, aspectRatio);
}

/**
 * Destructor
 */
I2D_ThreeLinkFish::~I2D_ThreeLinkFish()
{
	assert(shape != NULL);
	if (shape != NULL)
	{
		delete shape;
		shape = NULL;
	}
}

/**
 * Calculates CHI, the characteristic function, in the paper.
 */
void I2D_ThreeLinkFish::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	ThreeLinkFish::FillBlocks fill(eps,shape);
	block_processing.process(vInfo, coll, fill);
}

/**
 * Restart from previous file.
 *
 * @param t
 * @param filename
 */
void I2D_ThreeLinkFish::restart(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string()) // If the string is not set I write my own file
	{
		ppFile = fopen("restart_I2D_ThreeLinkFish.txt", "r");
		assert(ppFile!=NULL);
	}
	else // If string is set I open the corresponding file
	{
		ppFile = fopen(filename.c_str(), "r");
		assert(ppFile!=NULL);
	}

	shape->restart(ppFile);
	fclose(ppFile);
}

/**
 * Save state to file for restart later.
 *
 * @param t
 * @param filename
 */
void I2D_ThreeLinkFish::save(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string()) // If the string is not set I write my own file
	{
		ppFile = fopen("restart_I2D_ThreeLinkFish.txt", "w");
		assert(ppFile!=NULL);
	}
	else // If string is set I open the corresponding file
	{
		ppFile = fopen(filename.c_str(), "w");
		assert(ppFile!=NULL);
	}

	shape->save(ppFile);
	fclose(ppFile);
}

/**
 * When doing a restart, I need to back up to the correct time step for which
 * I have a restart file.
 *
 * @param t
 * @param filename
 */
void I2D_ThreeLinkFish::refresh(const double t, string filename)
{
	ifstream filestream;
	if(filename==std::string()) // Not set, open the file previously created
	{
		filestream.open("restart_I2D_ThreeLinkFish.txt");
		if (!filestream.good())
		{
			cout << "File not found. Exiting now." << endl;
			exit(-1);
		}
	}
	else
	{
		filestream.open(filename.c_str());
		if (!filestream.good())
		{
			cout << "File not found. Exiting now." << endl;
			exit(-1);
		}
	}

	// Open file and count number of lines
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();

	FILE * f = NULL;
	if(filename==std::string()) // If the string is not set I write my own file
	{
		f = fopen("restart_I2D_ThreeLinkFish.txt", "r");
		assert(f!=NULL);
	}
	else // If string is set I open the corresponding file
	{
		f = fopen(filename.c_str(), "r");
		assert(f!=NULL);
	}

	// Data is stored in dataVect until t=t_restart
	int N=0;
	vector < vector<float> > dataVect;
	for(int i=0; i<c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f,"%f",&variable);
		if (variable <= t)
		{
			const Real time = variable;
			row.push_back(time);
			fscanf(f,"%f",&variable);
			const Real xm = variable;
			row.push_back(xm);
			fscanf(f,"%f",&variable);
			const Real ym = variable;
			row.push_back(ym);
			fscanf(f,"%f",&variable);
			const Real um = variable;
			row.push_back(um);
			fscanf(f,"%f",&variable);
			const Real vm = variable;
			row.push_back(vm);
			fscanf(f,"%f",&variable);
			const Real thetam = variable;
			row.push_back(thetam);
			fscanf(f,"%f",&variable);
			const Real omegam = variable;
			row.push_back(omegam);
			fscanf(f,"%f",&variable);
			const Real J = variable;
			row.push_back(J);
			fscanf(f,"%f",&variable);
			const Real m = variable;
			row.push_back(m);
			fscanf(f,"%f",&variable);
			const Real rho = variable;
			row.push_back(rho);
			N=i;
		}
		else
		{
			break;
		}
		dataVect.push_back(row);
	}
	fclose(f);
	f = NULL;

	// dataVect is copied in a new "shape_001" file
	if(filename==std::string()) // If the string is not set I write my own file
	{
		f = fopen("restart_I2D_ThreeLinkFish.txt", "w");
		assert(f!=NULL);
	}
	else // If string is set I open the corresponding file
	{
		f = fopen(filename.c_str(), "w");
		assert(f!=NULL);
	}
	for (int i=0; i<N; i++)
	{
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e \n",
				dataVect[i][0], dataVect[i][1], dataVect[i][2], dataVect[i][3],
				dataVect[i][4], dataVect[i][5], dataVect[i][6], dataVect[i][7],
				dataVect[i][8], dataVect[i][9]);
	}
	fclose(f);
}

/**
 * Compute the velocity from Eq. (11) and (12) of the PoF paper. This is the
 * heavy lifting part of this class.
 *
 * @param t
 */
void I2D_ThreeLinkFish::computeDesiredVelocity(const double t)
{
	const int NQUANTITIES = 5;

	double mass = 0.0;
	double J = 0.0;
	double vxbar = 0.0;
	double vybar = 0.0;
	double omegabar = 0.0;

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(NQUANTITIES);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	/// Compute all
	ThreeLinkFish::ComputeAll computeAll(eps, shape, integrals, Uinf, nonempty);

	/// For each block, compute integrals for mass, moment of inertia, and
	/// Eqs. (11) and (12)
	block_processing.process(vInfo, coll, computeAll);

	/// Sum up contrubutions of each block
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		// ATCHIEU: first and second refer to <first, second> of map
		mass += (it->second)[0];
		J += (it->second)[1];
		vxbar += (it->second)[2]; // Eq. (11) of PoF
		vybar += (it->second)[3]; // Eq. (11) of PoF
		omegabar += (it->second)[4]; // Eq. (12) of PoF
	}

	/// @todo Put it super high for now
//	J = 10000.0;
//	mass = 10000.0;

	// Normalization (divide by mass or moment of inertia - in the most
	// complicated cases J varies and must be recomputed on the fly)
	vxbar /= mass;
	vybar /= mass;
	omegabar /= J;

	assert(vxbar < 1.0);
	assert(vybar < 1.0);
	assert(omegabar < 10);

	// Set the right um, vm and angular velocity for each single object
	shape->um = vxbar;
	shape->vm = vybar;
	shape->omegam = omegabar;
	shape->m = mass;
	shape->J = J;

	// Set desired velocities
	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
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

	ThreeLinkFish::FillVelblocks fillvelblocks(velblocks, xm_corr, ym_corr, vx_corr, vy_corr, omega_corr, eps, shape);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}

/**
 * Update and write to shape file at each timestep.
 *
 * @param dt
 * @param t
 * @param filename
 * @param _data
 */
void I2D_ThreeLinkFish::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	shape->updateIntertialQuantities(dt);
	shape->printShape(t, filename);
}

/**
 * Currently calculating uDEF on the fly given the exact deformation in a
 * vacuum.
 *
 * @param t
 */
void I2D_ThreeLinkFish::create(const double t)
{
	shape->updateShapeAndMotionInVacuum(t);
//	shape->updateShapeAndMotionInVacuumNumerical(t);

  ///---------------------------------------------------------------------------
	///----- CALCULATING THINGS NUMERICALLY --------------------------------------
	///---------------------------------------------------------------------------

	const int NQUANTITIES = 5;

	double mass = 0.0;
	double J = 0.0;
	double xbar = 0.0;
	double ybar = 0.0;
	double ubar = 0.0;
	double vbar = 0.0;
	double omegabar = 0.0;

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(NQUANTITIES);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	/// Corrections for the center of mass
	mass = J = xbar = ybar = ubar = vbar = omegabar = 0.0;
	ThreeLinkFish::ComputeCenterOfMass computeCenterOfMass(eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeCenterOfMass);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		xbar += (it->second)[1];
		ybar += (it->second)[2];
	}
	xbar /= mass;
	ybar /= mass;
	xm_corr = xbar;
	ym_corr = ybar;

	vector<double> printMe;
	printMe.push_back(xm_corr);
	printMe.push_back(ym_corr);
	printGeneralStuff(t, "CorrectionsMass.txt", printMe);


	/// Corrections to velocity and angular velocity
	/// (N.B. you have to calculate mass and J over again since your COM shifted)
	mass = J = xbar = ybar = ubar = vbar = omegabar = 0.0;
	ThreeLinkFish::ComputeAllCorrections computeAllCorrections(xm_corr, ym_corr, eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeAllCorrections);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		J += (it->second)[1];
		ubar += (it->second)[2];
		vbar += (it->second)[3];
		omegabar += (it->second)[4];
	}
	ubar /= mass;
	vbar /= mass;
	omegabar /= J;
	vx_corr = ubar;
	vy_corr = vbar;
	omega_corr = omegabar;

	printMe.clear();
	printMe.push_back(vx_corr);
	printMe.push_back(vy_corr);
	printMe.push_back(omega_corr);
	printGeneralStuff(t, "CorrectionsVelocity.txt", printMe);
}


void I2D_ThreeLinkFish::printGeneralStuff(const double t, const string filename, vector<double> stuffToPrint)
{
	int n = stuffToPrint.size();
	FILE * outfile = NULL;
	outfile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
	assert(outfile != NULL);
	for (int i = 0; i < n; i++)
		fprintf(outfile, "%16.6e", stuffToPrint[i]);
	fprintf(outfile, "\n");
	fclose(outfile);
}

