#include "PF_DipoleInlineBravais.h"

using namespace PF;

DipoleInlineBravais::DipoleInlineBravais(double a1, double a2, double phi, double theta, int nDesired)
		: n(0), a1(a1), a2(a2), phi(phi), theta(theta), nDesired(nDesired), r(1)
{
}

/**
 * Print the positions to a file.
 *
 * @param filename
 */
void DipoleInlineBravais::printPositions(const char* filename) const
{
	FILE* fid;
	fid = fopen(filename, "w");
	fprintf(fid, "%d\n", n);
	fprintf(fid, "%f\n", r);
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
void DipoleInlineBravais::readPositions(const char* filename)
{
	FILE* fid;
	fid = fopen(filename, "r");
	fscanf(fid, "%d", &n);
	fscanf(fid, "%f", &r);
	for (int i = 0; i < n; i++)
	{
		Real var1, var2;
		fscanf(fid, "%f %f", &var1, &var2);
		x.push_back(var1);
		y.push_back(var2);
	}
	fclose(fid);
}

void PF::DipoleInlineBravais::computeRadius()
{
	double rGuess = 0.5;
	double rOld = 0;
	double dr = rGuess;
	int nActual;

	printf("** Computing the smallest radius needed for nDesired = %d\n", nDesired);

	nActual = computeCoordinates(rGuess, false);
	if (nActual < nDesired)
	{
		printf("** At this size, swimmers cannot be fit in the (0,1) x (0,1) domain.\n");
		abort();
	}

	// Iterate to get it right
	const int maxIter = 100;
	int iter = 0;
	bool isDone = false;
	while (~isDone)
	{
		if (nActual == nDesired || iter > maxIter) break;

		dr = 0.5 * abs(rGuess - rOld);
		rOld = rGuess;
		if (nActual < nDesired)
		{
			rGuess += dr;
		}
		else if (nActual > nDesired)
		{
			rGuess -= dr;
		}

		nActual = computeCoordinates(rGuess, false);
		iter++;
//		printf("** nActual at iteration # %d (rGuess = %f): %d\n", iter, rGuess, nActual);
	}

	n = nActual;
	r = rGuess;

	nActual = computeCoordinates(r, true);
	printf("** Number of dipoles created in %d iterations (radius %f): %d\n", iter, r, n);
}

double PF::DipoleInlineBravais::computeCoordinates(double rAssign, bool isFillingVectors)
{
	r = rAssign;
	const double rInterior = 0.5*r;

	if (isFillingVectors)
	{
		x.clear();
		y.clear();
		isInterior.clear();
	}

	int count = 0;
	int countTrue = 0;

	const int lastRow = ceil(r / (a2 * sin(phi)));
	const int firstRow = -lastRow;

	for (int i = firstRow; i < lastRow; i++)
	{
		const double yi = i * a2 * sin(phi);
		const double shift = i * a2 * cos(phi);
		const double startx = shift - a1 * floor(shift / a1);

		bool isPastCircle = false;
		int j = 0;
		while (!isPastCircle) // go to left
		{
			const double xi = startx - j * a1;
			const double r2 = xi * xi + yi * yi;
			if (r2 >= r * r)
			{
				isPastCircle = true;
			}
			else if (r2 < r * r)
			{
				if (isFillingVectors)
				{
					x.push_back(xi * cos(theta) - yi * sin(theta));
					y.push_back(xi * sin(theta) + yi * cos(theta));

					// Determine interior swimmers
					if (r2 <= rInterior * rInterior)
					{
						isInterior.push_back(true);
						countTrue++;
					}
					else
					{
						isInterior.push_back(false);
					}
				}
				count++;
			}
			j++;
		}

		isPastCircle = false;
		j = 1;
		while (!isPastCircle) // go to right
		{
			const double xi = startx + j * a1;
			const double r2 = xi * xi + yi * yi;
			if (r2 >= r * r) isPastCircle = true;
			else if (r2 < r * r)
			{
				if (isFillingVectors)
				{
					x.push_back(xi * cos(theta) - yi * sin(theta));
					y.push_back(xi * sin(theta) + yi * cos(theta));

					if (r2 <= rInterior * rInterior)
					{
						isInterior.push_back(true);
						countTrue++;
					}
					else
					{
						isInterior.push_back(false);
					}
				}
				count++;
			}
			j++;
		}
	}

	if (isFillingVectors) printf("** Number in averaging window = %d\n", countTrue);


	return count;
}
