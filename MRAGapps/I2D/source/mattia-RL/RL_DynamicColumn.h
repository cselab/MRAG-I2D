/*
 * RL_DynamicColumn.h
 *
 *  Created on: Apr 5, 2012
 *      Author: dalmassg
 */

#ifndef RL_DYNAMICCOLUMN_H_
#define RL_DYNAMICCOLUMN_H_

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Agent.h"

namespace RL
{

class RL_DynamicColumn : public RL_Agent
{
protected:
	unsigned int count;

public:
	const Real D;
	Real x,y;
	RL_DynamicColumn(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID = 0, RL_TabularPolicy ** _policy = NULL);
	virtual ~RL_DynamicColumn();

	virtual void update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_DYNAMICCOLUMN_H_ */
