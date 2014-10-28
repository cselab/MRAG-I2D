/*
 * RL_CircularWall.h
 *
 *  Created on: Mar 23, 2012
 *      Author: mgazzola
 */

#ifndef RL_CIRCULARWALL_H_
#define RL_CIRCULARWALL_H_

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Agent.h"

namespace RL
{

class RL_CircularWall : public RL_Agent
{
public:
	const Real D,innerD,x,y;
	RL_CircularWall(MRAG::ArgumentParser & parser);
	virtual ~RL_CircularWall();

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_CIRCULARWALL_H_ */
