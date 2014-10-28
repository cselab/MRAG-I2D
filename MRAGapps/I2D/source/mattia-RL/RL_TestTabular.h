#pragma once

#include <string>
#include "RL_Environment.h"
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "IF2D_Test.h"
#include "FieldViewer.h"
#include "RL_TabularPolicy.h"
#include "RL_Agent.h"

namespace RL
{

class RL_TestTabular: public IF2D_Test
{
protected:
	Real CHARLENGTH, XPOS, YPOS;
	bool RESTART;
	int SAVEFREQ;

	MRAG::ArgumentParser parser;
	RL_TabularPolicy * policy;
	RL_Agent * agent;
	MRAG::Profiler profiler;

	void _save();
	void _dispose();
	void _refresh();
	void _prepareAgents();

#ifdef _RL_VIZ
	//FieldViewer field_viewer;
	void _paint();
#endif

public:
	
	RL_TestTabular(const int argc, const char ** argv);
	~RL_TestTabular();
	
	void run();
	void paint(){};
};

}
