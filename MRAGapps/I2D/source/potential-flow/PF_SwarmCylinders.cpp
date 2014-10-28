/*
 * PF_SwarmCylinders.cpp
 *
 *  Created on: May 29, 2012
 *      Author: atchieu
 */

#include "PF_SwarmCylinders.h"

/**
 * Setup the swarm.
 *
 * @param n_
 * @param radius_
 * @param rmin_
 * @param A_
 * @param K1_
 * @param K2_
 * @param alpha_
 */
SwarmCylinders::SwarmCylinders(int n_, double radius_, double rmin_, double A_, double K1_, double K2_, double alpha_) :
		n(n_), radius(radius_), rmin(rmin_), A(A_), K1(K1_), K2(K2_), alpha(alpha_), spline(A_, K1_, K2_, alpha_), dt(1e-3), threshold(5e-6), factorCylinder(0.001), factorBoundary(0.05)
{
	/// Note the default 'sim' parameters given in the initializer list.
	/// These numbers are configured for the case where we have n = 100.
}

/**
 * Setup the spline.
 *
 * @param desiredArea
 * @param k1
 * @param k2
 * @param alpha
 */
SwarmCylinders::Spline::Spline(double desiredArea, double k1, double k2, double alpha)
{
	assert(alpha >= 0.0 && alpha <= M_PI);

	double r = sqrt(0.35e2) * sqrt(0.1e1 / ((double) ((9 * k1 - 13 * k2 * k2 - 9 * k1 * k2 + 13) * alpha) + 0.13e2 * 0.3141592654e1 * (double) (k1 * k1) + 0.13e2 * 0.3141592654e1 * (double) (k2 * k2) + 0.9e1 * 0.3141592654e1 * (double) k1 * (double) k2) * desiredArea);
	double r1 = r * k1;
	double r2 = r * k2;

	assert(r > 0.0);
	assert(r1 > 0.0);
	assert(r2 > 0.0);

	m_a0 = r;
	m_b0 = 0;
	m_c0 = 0.3e1 * pow(alpha, -0.2e1) * r * k1 - 0.3e1 * pow(alpha, -0.2e1) * r;
	m_d0 = -0.2e1 * pow(alpha, -0.3e1) * r * k1 + 0.2e1 * pow(alpha, -0.3e1) * r;
	m_a1 = 0.3141592654e1 * 0.3141592654e1 * (0.3141592654e1 - 0.3e1 * alpha) / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k1 + alpha * alpha * (0.3e1 * 0.3141592654e1 - alpha) / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k2;
	m_b1 = 0.6e1 * 0.3141592654e1 * alpha / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k1 - 0.6e1 * 0.3141592654e1 * alpha / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k2;
	m_c1 = -0.3e1 * (0.3141592654e1 + alpha) / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k1 + 0.3e1 * (0.3141592654e1 + alpha) / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k2;
	m_d1 = 0.2e1 / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k1 - 0.2e1 / (0.3141592654e1 - alpha) / (0.3141592654e1 * 0.3141592654e1 - 0.2e1 * 0.3141592654e1 * alpha + alpha * alpha) * r * k2;

	m_alpha = alpha;
	m_rmin = min(min(r, r1), r2);
}

/**
 * Evaluate the spline at a given angle theta. Returns back the value of r.
 * Remember to convert to Cartesian if you need to.
 *
 * @param theta
 * @return
 */
double SwarmCylinders::Spline::evaluate(const double theta) const
{
	double a = theta;
	if (theta > M_PI)
		a = 2.0 * M_PI - theta;

	if (a <= m_alpha)
		return m_a0 + m_b0 * a + m_c0 * a * a + m_d0 * a * a * a;
	else
		return m_a1 + m_b1 * a + m_c1 * a * a + m_d1 * a * a * a;
}

/**
 * Compute the area as a diagnostic and print it stdout.
 */
void SwarmCylinders::Spline::computeArea()
{
	double da = 1e-3;
	vector<double> x;
	vector<double> y;

	for (int i = 0; i < 1.0 / da; i++)
	{
		vector<double> a;

		a.push_back(i * da * M_PI);
		vector<double> r;

		r.push_back(evaluate(a[0]));
		convertToCartesian(r, a);
		x.push_back(r[0]);
		y.push_back(a[0]);
		a.clear();
		r.clear();
	}

	double area = 0;

	const int xsize = x.size();
	for (int i = 0; i < xsize - 1; i++)
		area += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);

	printf("** Area: %1.10f", -area);
}

/**
 * Converts from Cartesian to polar coordinates. In place transformation.
 *
 * @param x
 * @param y
 */
void SwarmCylinders::convertToPolar(vector<double>& x, vector<double>& y)
{
	const int xsize = x.size();
	for (int i = 0; i < xsize; i++)
	{
		double r = sqrt(x[i] * x[i] + y[i] * y[i]);
		double theta = 0.0;

		if (x[i] && y[i] == 0.0)
			theta = 0.0;
		else if (x[i] > 0.0 && y[i] >= 0.0)
			theta = atan(y[i] / x[i]);
		else if (x[i] > 0.0 && y[i] < 0.0)
			theta = atan(y[i] / x[i]) + 2.0 * M_PI;
		else if (x[i] < 0.0)
			theta = atan(y[i] / x[i]) + M_PI;
		else if (x[i] == 0.0 && y[i] > 0.0)
			theta = 0.5 * M_PI;
		else if (x[i] == 0.0 && y[i] < 0.0)
			theta = 1.5 * M_PI;

		x[i] = r;
		y[i] = theta;
	}
}

/**
 * Convert from polar coordinates to Cartesian. In place transformation.
 *
 * @param x
 * @param y
 */
void SwarmCylinders::convertToCartesian(vector<double>& x, vector<double>& y)
{
	const int xsize = x.size();
	for (int i = 0; i < xsize; i++)
	{
		double r = x[i];
		double theta = y[i];

		x[i] = r * cos(theta);
		y[i] = r * sin(theta);
	}
}

/**
 * Print the spline in (x,y) coordinates to a file named spline.txt.
 */
void SwarmCylinders::printSpline() const
{
	double dalpha = 2.0 * M_PI / (double) M;

	vector<double> xs;
	vector<double> ys;

	for (int i = 0; i < M; i++)
	{
		xs.push_back(spline.evaluate(((double) i + 0.5) * dalpha));
		ys.push_back(((double) i + 0.5) * dalpha);
	}
	convertToCartesian(xs, ys);

	FILE* fid;
	fid = fopen("spline.txt", "w");
	for (int i = 0; i < (int) xs.size(); i++)
	{
		fprintf(fid, "%10.6f %10.6f\n", xs[i], ys[i]);
	}
	fclose(fid);
}

/**
 * Print the position of the cylidners to a file.
 *
 * @param filename
 */
void SwarmCylinders::printPositions(const char* filename) const
{
	FILE* fid;
	fid = fopen(filename, "w");
	for (int i = 0; i < n; i++)
	{
		fprintf(fid, "%e %e\n", x[i], y[i]);
	}
	fclose(fid);
}

/**
 * Read positions from a previously created file.
 *
 * @param filename
 */
void SwarmCylinders::readPositions(const char* filename)
{
	FILE* fid;
	fid = fopen(filename, "r");
	for (int i = 0; i < n; i++)
	{
		Real var1, var2;
		fscanf(fid, "%f %f", &var1, &var2);
		x.push_back(var1);
		y.push_back(var2);
	}
	fclose(fid);
}

/**
 * Initially place the cylinders somewhere inside the region bounded by the
 * spline.
 */
void SwarmCylinders::placeCylinders()
{
	printf("** Number of agents: %d\n", n);
	printf("** Placing agents %1.4f inside of the boundary! \n", rmin);

	double dalpha = 2.0 * M_PI / (double) n;
	for (int i = 0; i < n; i++)
	{
		x.push_back(spline.evaluate(((double) i + 0.5) * dalpha) - rmin);
		y.push_back(((double) i + 0.5) * dalpha);
	}
	convertToCartesian(x, y);
}

/**
 * Find the equilibrium after placing the cylinders in the domain somewhere.
 * Iterates until tolerence is met.
 */
void SwarmCylinders::findEquilibrium()
{
	double dalpha = 2.0 * M_PI / (double) M;

	vector<double> m;

	for (int i = 0; i < M; i++)
	{
		m.push_back(spline.evaluate(i * dalpha) * dalpha); /// assuming constant density charge on boundary
	}

	vector<double> xs;
	vector<double> ys;

	for (int i = 0; i < M; i++)
	{
		xs.push_back(spline.evaluate(((double) i + 0.5) * dalpha));
		ys.push_back(((double) i + 0.5) * dalpha);
	}

	convertToCartesian(xs, ys);

	int nIter = 0;
	bool steadyState = false;
	while (!steadyState)
	{
		vector<double> xNew;
		vector<double> yNew;
		const int xsize = x.size();
		for (int i = 0; i < xsize; i++)
		{
			xNew.push_back(0.0);
			yNew.push_back(0.0);
		}

		for (int i = 0; i < xsize; i++)
		{
			double u = 0.0;
			double v = 0.0;

			for (int j = 0; j < xsize; j++) /// from each body
			{
				if (i != j)
				{
					double norm = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
					double factor = factorCylinder / (norm * norm * norm);
					u += (x[i] - x[j]) * factor;
					v += (y[i] - y[j]) * factor;
				}
			}

			for (int j = 0; j < M; j++) /// from boundary
			{
				double norm = sqrt((x[i] - xs[j]) * (x[i] - xs[j]) + (y[i] - ys[j]) * (y[i] - ys[j]));
				double factor = factorBoundary * m[j] / (norm * norm * norm);
				u += (x[i] - xs[j]) * factor;
				v += (y[i] - ys[j]) * factor;
			}

			xNew[i] = x[i] + dt * u;
			yNew[i] = y[i] + dt * v;
		}

		convertToPolar(xNew, yNew);

		for (int i = 0; i < xsize; i++)
		{
			const double rMax = spline.evaluate(yNew[i]) - radius;
			if (rMax < xNew[i])
				xNew[i] = 0.1 * rMax; /// move it toward (0,0)

//			xNew[i] = min(spline.evaluate(yNew[i])-radius, xNew[i]); /// make it sit on the boundary if it is not repelled enough (not exactly correct)

		}

		convertToCartesian(xNew, yNew);

		double maxNorm = 0.0;
		for (int i = 0; i < xsize; i++)
		{
			double norm = sqrt((xNew[i] - x[i]) * (xNew[i] - x[i]) + (yNew[i] - y[i]) * (yNew[i] - y[i])); /// relative difference
			maxNorm = max(norm, maxNorm);
		}

		for (int i = 0; i < xsize; i++)
		{
			x[i] = xNew[i];
			y[i] = yNew[i];
		}

		// TODO Remove this testing shit, writing the positions each iteration

		if (nIter % 10 == 0)
		{
			FILE* fidx;
			FILE* fidy;

			if (nIter == 0)
			{
				fidx = fopen("x.txt", "w");
				fidy = fopen("y.txt", "w");
			}
			else
			{
				fidx = fopen("x.txt", "a");
				fidy = fopen("y.txt", "a");
			}

			const int xsize = x.size();
			for (int i = 0; i < xsize; ++i)
			{
				fprintf(fidx, "%e ", x[i]);
				fprintf(fidy, "%e ", y[i]);
			}

			fprintf(fidx, "\n");
			fprintf(fidy, "\n");

			fclose(fidx);
			fclose(fidy);
		}
		nIter++;

		if (maxNorm <= threshold)
			break;

		if (nIter % 5000 == 0)
			printf("** Norm at iteration %d: %e\n", nIter, maxNorm);

		if (nIter > 50000)
		{
			printf("** Something is wrong, swarm shape is not converging quickly enough... aborting... \n");
			abort();
		}
	}
	printf("** Number of iterations to steady state formation: %d\n", nIter);
}

/**
 * Check for a collision based on radius, if there is one then quit the
 * simulation.
 */
void SwarmCylinders::checkCollision() const
{
	const int xsize = x.size();
	for (int i = 0; i < xsize; i++)
		for (int j = 0; j < xsize; j++)
		{
			if (i != j)
			{
				double dist = sqrt(pow(x[i] - x[j], 2.0) + pow(y[i] - y[j], 2.0));
				if (dist < 2.0 * radius)
				{
					printf("** For K1, K2, alpha: %e, %e, %e\n", K1, K2, alpha);
					printf("** There is a collisions between cylinders somewhere!\n");
					exit(-1);
				}
			}
		}
}

/**
 * Check if the spline crosses the y-axis, if so, abort.
 */
void SwarmCylinders::checkValidityOfSpline() const
{
	const int nSpline = 512;
	const double ymin = 0.10; // if 0, prevents you from creating splines that cross the y-axis
	const double dalpha = 2.0 * M_PI / (double) nSpline;

	vector<double> xTest;
	vector<double> yTest;

	for (int i = 0; i < nSpline; i++)
	{
		xTest.push_back(spline.evaluate((double) i + 0.5) * dalpha);
		yTest.push_back(((double) i + 0.5) * dalpha);
	}
	convertToCartesian(xTest, yTest);

	for (int i = 0; i < nSpline; i++)
	{
		if (yTest[i] < ymin)
		{
			printf("** The spline crosses or is too close to the the mirror plane! (ymin = %0.4f, k1 = %0.4f k2 = %0.4f, alpha = %0.4f)\n", ymin, K1, K2, alpha);
			abort();
		}
	}
}

/**
 * Rotate agents to point up as commonly used in our simulations.
 */
void SwarmCylinders::rotateClockwiseQuarterTurn()
{
	vector<double> xtemp;
	vector<double> ytemp;

	/// Rotate to point swarm up
	for (int i = 0; i < n; i++)
	{
		xtemp.push_back(y[i]);
		ytemp.push_back(-x[i]);
	}

	x = xtemp;
	y = ytemp;
}

