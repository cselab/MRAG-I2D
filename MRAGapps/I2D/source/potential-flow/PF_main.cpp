#include "RL_Environment.h"
#include "PF_Solver.h"

#ifdef _RL_VIZ
#ifdef __APPLE__
#include "GLUT/glut.h"
#endif
#endif

using namespace MRAG;
using namespace RL;
using namespace PF;
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

	static void setup(int argc, const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(800, 800); //1024 //512
		glutCreateWindow("School");
		glutDisplayFunc(display);
		glClearColor(1, 1, 1, 1); // (0, 0, 0, 1) for black

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0, 1.0, 0, 1.0, 1.0, -1.0);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
};
#endif

int main(int argc, const char ** argv)
{
	MRAG::ArgumentParser parser(argc, argv);

	Environment::setup(max(1, parser("-nthreads").asInt()));

	if (parser("-study").asString() == "LEARNING_DIPOLE") test = new PF_Solver(argc, argv);
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

	return (0);
}
