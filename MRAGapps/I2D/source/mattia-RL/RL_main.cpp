#include "RL_Environment.h"
#include "RL_TestTabular.h"

#ifdef _RL_VIZ
#ifdef __APPLE__
#include "GLUT/glut.h"
#endif
#endif

using namespace MRAG;
using namespace RL;
using namespace std;

IF2D_Test * test = NULL;

#ifdef _RL_VIZ
struct VisualSupport
{
	static void display()
	{
	}

	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		test->run();
		glutSwapBuffers();
	}

	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;

		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}

		glutDisplayFunc(display);
		glutIdleFunc(idle);

		glutMainLoop();
	}

	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		//glutInitWindowSize(1024,1024);
		glutInitWindowSize(700,700);
		glutCreateWindow("School");
		glutDisplayFunc(display);
		//glClearColor(1,1,1,1);
		glClearColor(0,0,0,1);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, 1.0, 0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
};
#endif

int main (int argc, const char ** argv)
{
	ArgumentParser parser(argc, argv);

	Environment::setup(max(1, parser("-nthreads").asInt()));

	if( parser("-study").asString() == "RL" )
		test = new RL_TestTabular(argc, argv);
	else
	{
		printf("Study case not defined!\n");
		abort();
	}

#ifdef _RL_VIZ
	VisualSupport::run(argc, argv);
#else
	test->run();
#endif

	return 0;
}
