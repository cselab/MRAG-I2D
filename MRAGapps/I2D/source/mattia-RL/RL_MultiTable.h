/*
 * RL_MultiTable.h
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */

#ifndef RL_MULTITABLE_H_
#define RL_MULTITABLE_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

namespace RL
{

class RL_MultiTable
{
	vector<int> dim;
	map< vector<int>, double > data;

	bool _check_bounds(vector<int> idx) const;
	int _lines(const char * const filename);

public:
	// Costructor-Destructor
	RL_MultiTable(){}
	RL_MultiTable(vector<int> dim);
	~RL_MultiTable();

	// Methods
	inline void setdim(vector<int> _dim){ data.clear(); dim.clear(); dim = _dim; }
	double & operator()(const vector<int> & idx);
	double read(const vector<int> & idx);
	double usage() const;
	void save(string name = "savedQ");
	void restart(string name = "savedQ");
};

}

#endif /* RL_MULTITABLE_H_ */
