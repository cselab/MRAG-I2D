//#include "I2D_DaFish.h"
//#include "I2D_FlowPastObstacleRK.h"
#include "I2D_FlowPastObstacle_Gudonov.h"

using namespace MRAG;
using namespace std;

I2D_Test * test = NULL;

#ifdef _MRAG_GLUT_VIZ 
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
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("Fish");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(0.0, 1.0, 0.0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);
		
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}
};
#endif

int main (int argc, const char ** argv) 
{
	ArgumentParser parser(argc, argv);
	
	Environment::setup(max(1, parser("-nthreads").asInt()));
	
	//test = new I2D_FlowPastObstacleRK(argc, argv);
	test = new I2D_FlowPastObstacle_Gudonov(argc, argv);	
	//test = new I2D_DaFish(argc, (char **)argv);
	
	tbb::tick_count t1,t0;
	
	{
		t0=tbb::tick_count::now();
		
#ifdef _MRAG_GLUT_VIZ 
		VisualSupport::run(argc, argv);
#else
		test->run();
#endif
		
		t1=tbb::tick_count::now();
	}
	
	printf("we spent: %2.2f \n",(t1-t0).seconds());
	
    return 0;
}
