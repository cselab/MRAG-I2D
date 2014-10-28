/*
 * RL_ObjectFactory.h
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */

#ifndef RL_OBJECTFACTORY_H_
#define RL_OBJECTFACTORY_H_

#include <map>
#include <vector>
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Agent.h"
#include "RL_TabularPolicy.h"

namespace RL
{

class RL_ObjectFactory
{
	const Real charLength, XCM, YCM;

	int _lines(string filename);

public:
	RL_ObjectFactory(const Real charLength, const Real XCM, const Real YCM);
	~RL_ObjectFactory();

	void create(MRAG::ArgumentParser & parser, map< string, vector<RL_Agent *> >& shapesMap, RL_TabularPolicy ** policy = NULL);
};

} /* namespace RL */
#endif /* RL_OBJECTFACTORY_H_ */
