/*
 * RL_MultiTable.cpp
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "RL_MultiTable.h"

namespace RL
{

RL_MultiTable::RL_MultiTable(vector<int> dim):dim(dim)
{
}

RL_MultiTable::~RL_MultiTable()
{
}

bool RL_MultiTable::_check_bounds(vector<int> idx)const
{
	const unsigned int DIM = dim.size();
	assert( DIM > 0 );
	assert( idx.size() == DIM );

	for(unsigned int i=0; i<DIM; i++)
		assert( idx[i]>=0 && idx[i]<dim[i] );

	return true;
}

double & RL_MultiTable::operator()(const vector<int> & idx)
{
	assert(_check_bounds(idx));
	return data[idx];
}

double RL_MultiTable::read(const vector<int> & idx)
{
	assert(_check_bounds(idx));
	return data[idx];
}

double RL_MultiTable::usage() const
{
	const int n = data.size();
	int maxn = 1;
	for(int i=0; i<(int)dim.size(); i++)
		maxn *= dim[i];

	return (double)n/(double)maxn;
}

int RL_MultiTable::_lines(const char * const filename)
{
	ifstream filestream(filename);

	int c = 0;
	string line;

	//count the lines
	while( getline(filestream, line) )
		c++;

	return c;
}

void RL_MultiTable::save(string name)
{
	printf("save %s\n", name.c_str());

	const int DIM = dim.size();
	string nameBackup = name + "_backup";
	ofstream out(nameBackup.c_str());
	out.precision(20);

	out << DIM << "  ";
	for(int i=0; i<DIM; i++)
		out << dim[i] << "  ";

	out << endl;

	unsigned int counter = 0;
	for(map< vector<int>, double >::iterator it = data.begin(); it!=data.end(); it++)
	{
		for(int i=0; i<DIM; i++)
			out << (*it).first[i] << "  ";

		out << scientific << (*it).second << endl;
		out.flush();
		counter++;
	}
	out.flush();
	out.close();

	const int ndata = data.size();
	const int nlinesBackup = _lines(nameBackup.c_str())-1;

	// Prepare copying command
	std::string command = "cp ";
	std::string nameOriginal = name;
	command = command + nameBackup + " " + nameOriginal;

	// Submit the command to the system
	FILE *ptr = popen(command.c_str(), "r");
	//if( ptr != NULL){ std::cout << command << std::endl; std::cout << "successful policy backup!" << std::endl; }
	//else{ std::cout << command << std::endl; std::cout << "policy backup failed" << std::endl; }
	pclose( ptr );

	const int nlines = _lines(name.c_str())-1;

	// Check consistency
	printf("counter=%d, ndata=%d, nlines=%d, nlinesBackup=%d\n", counter, ndata, nlines, nlinesBackup);
	if(ndata!=nlines || nlines!=counter || counter!=ndata || nlinesBackup!=nlines)
	{
		printf("Something went fucking wrong!\n");
		abort();
	}
}

void RL_MultiTable::restart(string name)
{
	const int DIM = dim.size();
	if(DIM==0){ printf("Policy dimension was not set!\n"); abort(); }

	string nameBackup = name + "_backup";

	const int nlinesBackup = _lines(nameBackup.c_str())-1;

	if(nlinesBackup>=1)
	{
		ifstream in(nameBackup.c_str());
		printf("%s\n", nameBackup.c_str());

		if(in.good())
		{
			int dimSignature = 0;
			in >> dimSignature;
			if(DIM!=dimSignature){ printf("Saved policy and current policy do not match in dimensionality!\n"); abort(); }
			else{ printf("dimSignature=%d\n",dimSignature); }

			for(int i=0; i<DIM; i++)
			{
				int dummy = 0;
				in >> dummy;
				if(dummy!=dim[i]){ printf("Saved policy and current policy do not match in dimensionality!\n"); abort(); }
				else{ printf("dim[i]=%d\n",dim[i]); }
			}

			unsigned counter = 0;
			for(int k = 0; k<nlinesBackup; k++)
			{
				vector<int> key(DIM,0.0);

				for(int i=0; i<DIM; i++)
					in >> key[i];

				double dummy = 0.0;
				in >> dummy;
				data[key] = dummy;
				counter++;

				for(int i=0; i<DIM; i++)
					printf("%d ", key[i]);

				printf("---> %10.10e\n", dummy);
			}

			// Check consistency
			const int ndata = data.size();
			printf("counter=%d, ndata=%d, nlinesBackup=%d\n", counter, ndata, nlinesBackup);
			if(ndata!=nlinesBackup || nlinesBackup!=counter || counter!=ndata)
			{
				printf("Something went fucking wrong!\n");
				abort();
			}
		}
		else
		{
			printf("WTF couldnt open file (ok keep going mofo)!\n");
		}

		in.close();
	}
}

}

