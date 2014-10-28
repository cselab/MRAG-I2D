/*
 *  PF_ObjectFactory.h
 *  DipoleCode
 *
 *  Created by Alexia on 4/4/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <map>
#include <vector>
#include "RL_Environment.h"
#include "RL_TabularPolicy.h"
#include "PF_Agent.h"

namespace PF
{

class PF_ObjectFactory
{
public:
	PF_ObjectFactory(const Real charLength, const Real XCM, const Real YCM);
	~PF_ObjectFactory(){};

	void create(MRAG::ArgumentParser & parser, map< string, vector<PF_Agent*> >& shapesMap, const Real lr, const Real gamma, RL::RL_TabularPolicy ** policy = NULL);
	Real getUnitVelocity();

private:
	const Real charLength, XCM, YCM;
	Real unitVelocity;

	int _lines(string filename);
	void maximumVelocity(Real v);
};

} /*namespace PF*/
