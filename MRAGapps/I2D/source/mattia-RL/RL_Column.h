/*
 * RL_Column.h
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */

#ifndef RL_COLUMN_H_
#define RL_COLUMN_H_

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Agent.h"

namespace RL
{

class RL_Column : public RL_Agent
{
public:
	const Real D,x,y;
	RL_Column(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID = 0, RL_TabularPolicy ** _policy = NULL);
	virtual ~RL_Column();

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_COLUMN_H_ */
