/*
 * PF_SwarmCylinders.h
 *
 *  Created on: May 29, 2012
 *      Author: atchieu
 */

#ifndef PF_SWARMCYLINDERS_H_
#define PF_SWARMCYLINDERS_H_

#include "RL_Environment.h"
#include <vector>
#include <cmath>

class SwarmCylinders
{
public:
	SwarmCylinders(int n_, double radius_, double rmin_, double A_, double K1_, double K2_, double alpha_);

	class Spline
	{
	public:
		Spline() {};
		Spline(double desiredArea, double k1, double k2, double alpha);
		double evaluate(const double theta) const;
		void computeArea();

	protected:
		double m_a0;
		double m_b0;
		double m_c0;
		double m_d0;
		double m_a1;
		double m_b1;
		double m_c1;
		double m_d1;
		double m_alpha;
		double m_rmin;
	};

	static void convertToPolar(vector<double> &x, vector<double> &y);
	static void convertToCartesian(vector<double> &x, vector<double> &y);

	void printSpline() const;
	void printPositions(const char* filename) const;
	void readPositions(const char* filename);

	void placeCylinders();
	void findEquilibrium();
	void checkCollision() const;
	void checkValidityOfSpline() const;
	void rotateClockwiseQuarterTurn();

	vector<double> getX() const
	{
		return x;
	}
	vector<double> getY() const
	{
		return y;
	}

protected:
	Spline spline;

	static const int M = 512; /// number of points on boundary
	int n; /// number of cylinders
	const double radius; /// cylinder radius defines min distance to wall of spline
	const double dt, threshold; /// steady state simulation parameters
	const double A, K1, K2, alpha, rmin; /// spline paramters

	const double factorCylinder; /// hard coded, scales force from each body
	const double factorBoundary; /// hard coded, scales force from boundary

	vector<double> x;
	vector<double> y;
};

#endif /* PF_SWARMCYLINDERS_H_ */
