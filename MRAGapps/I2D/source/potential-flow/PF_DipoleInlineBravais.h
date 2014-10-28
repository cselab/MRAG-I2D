#pragma once

#include "RL_Environment.h"
#include <vector>
#include <cmath>

namespace PF
{

class DipoleInlineBravais
{
public:
	DipoleInlineBravais(double a1, double a2, double phi, double theta, int nDesired);

	int getn() const { return n; }
	vector<double> getX() const { return x; }
	vector<double> getY() const { return y; }
	vector<bool> getIsInterior() const { return isInterior; }

	void printPositions(const char* filename) const;
	void readPositions(const char* filename);

	void computeRadius();
	int computeN(double rGuess);
	double computeCoordinates(double rGuess, bool isFillingVectors);

protected:
	int n; /// actual number of swimmers
	int nDesired; /// number of desired swimmers
	const double a1;
	const double a2;
	const double phi;
	const double theta; /// angle to incline the lattice
	double r;

	vector<double> x;
	vector<double> y;
	vector<bool> isInterior;
};

}
