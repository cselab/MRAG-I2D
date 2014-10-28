/*
 * RL_PolicyBasClass.h
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */

#ifndef RL_TABULARPOLICY_H_
#define RL_TABULARPOLICY_H_

#include <assert.h>
#include <vector>
#include <string>
#include "RL_MultiTable.h"

namespace RL
{

class RL_TabularPolicy
{
protected:
	vector<int> dim;
	RL_MultiTable Q;
	RL_MultiTable Qusage;
public:
	RL_TabularPolicy(){};
	virtual ~RL_TabularPolicy(){};

	virtual int selectAction(const vector<int> & state) = 0;
	virtual void setReward(double reward) = 0;
	virtual void setStateActionStart( const vector<int> & idx ) = 0;
	virtual void setStateEnd( const vector<int> & idx ) = 0;
	virtual void update(string name = "learning") = 0;

	void save(string name = "savedQ")
	{
		assert(name!=string());
		Q.save(name);

		string usage(name+"_usage");
		Qusage.save(usage);

	};

	void restart(string name = "savedQ")
	{
		assert(name!=string());
		Q.restart(name);

		//string usage(name+"_usage");
		//Qusage.restart(usage);
	};

	inline void setdim(vector<int> _dim)
	{
		dim.clear();
		dim = _dim;
		Q.setdim(dim);
		Qusage.setdim(dim);
	}

	inline bool samedim(vector<int> _dim) const
	{
		if( _dim.size()!=dim.size() )
			return false;

		bool rightdim = true;
		for(int i=0; i<(int)dim.size(); i++)
			rightdim *= (_dim[i]==dim[i]);

		return rightdim;
	}
};

}

#endif /* RL_TABULARPOLICY_H_ */
