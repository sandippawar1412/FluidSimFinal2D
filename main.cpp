#include "commonData.h"
#include "keyboard.h"
#include "main.h"
#include <time.h>
#include <pthread.h>
#include <omp.h>
#include <sys/time.h>


//#include "viennacl/vector.hpp"

#define WG  //WG - With Graphics
//#undef WG

int grid_size;
int nthreads;
#ifndef WG
int main(int argc, char** argv)
{
	grid_size = GRID_SIZE;
	nthreads = NTHREADS;	
	if(argc>=2){
	nthreads = atoi(argv[1]);
	grid_size = atoi(argv[2]);
	}
 	omp_set_num_threads(nthreads);
	struct timeval tt1, tt2;
	initMain();
	int i=pthread_getconcurrency();
	int it=0;
	int ni= 0;
	while(ni++<NO_OF_ITERATIONS){
		gettimeofday(&tt1, NULL);
		animate();
		gettimeofday(&tt2, NULL);
		int milliSeconds = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
		cout<<"Iteration "<<it<<" : "<<milliSeconds<<"ms"<<endl<<endl;
		it++;
	}
	return 0;
}
#endif


void animate()
{
	fluidSim->simulate(timestep);
} 

void initMain(void) 
{
   sGrid->initGridStag();
   fluidSim->init(sGrid);
   fluidSim->initFluidBody(fBT);// 2: indicates dam break center
   print->init(sGrid);
   render->init(sGrid);
   
}

void initParticles()
{
	sGrid->fluidParticles.clear();
}

/*------------------Code Below this is used when rendered is ON---------*/
void display ( void ) ;
void preDisplay ( );
void postDisplay ( );
void openGlutWindow ( char* windowName) ;
void reshape ( int w, int h ) ;
void idleFun();

#ifdef WG
int main(int argc, char** argv)
{
	grid_size = GRID_SIZE;
	nthreads = NTHREADS;	
	if(argc>=2){
	nthreads = atoi(argv[1]);
	grid_size = atoi(argv[2]);
	}
   	omp_set_num_threads(nthreads);
   
   /*#pragma omp parallel
   {
		cout<<"Hello"<<endl;
   }*/
   //int i=0;
   /*matrix<double > temp;
   temp.resize(10,10);
   #pragma omp parallel for
   	for(int i=0;i<10;i++)
		for(int j=0;j<10;j++)
			cout<<temp(i,j);
   */
  
   glutInit(&argc, argv);
   initMain ();
   int i=pthread_getconcurrency();

   char windowName[]="   Liquid_Simulator-LevelSet+Surface" ;
   openGlutWindow(windowName);

   glutMainLoop();
   
   return 0; 
}
#endif
	
void display(void){
	preDisplay();

	static bool flag[10]={false,false,false,false,false,true};//boundary grid particles surface vector mat

	char output1 = ' ';
	bool anyUpdation = false;

	switch(output1){
		case'~':
			flag[0] = !flag[0];
			//render->renderBoundary
			break;
		case '!':
			flag[1] = !flag[1];
			//render->renderGrid();
			break;
		case '@':
			flag[2] = !flag[2];
			//render->renderParticles();
			break;
		case '#':
			flag[3] = !flag[3];
			//render->renderSurfaceBoundary();
			break;
		case '$':
			flag[4] = !flag[4];
			//render->renderVector2D(sGrid->u,sGrid->v);
			break;
		case '%':
			flag[5] = !flag[5];
			//render->renderMat(sGrid->distanceLevelSet,1);
			break;
		/*case '^':
			flag[6] = !flag[6];
			//render->renderMat(sGrid->isFluidBoundary,1);
			break;*/
	}
	if(anyUpdation){
		cout<<"Flags :"<<" ~"<<flag[0]<<" !"<<flag[1]<<" @"<<flag[2]<<" #"<<flag[3]<<
			         " $"<<flag[4]<<" %"<<flag[5]/*<<" ^"<<flag[6]*/<<"   +"<<endl;
		anyUpdation = false;
	}
	extern int swich;
	if(swich==0){
		
		//render->renderBoundary();
		render->renderGrid();
		render->renderParticles();
		render->renderSurfaceBoundary();
	}
	else if(swich==1){
		render->renderMat(sGrid->p,2);
	}
	else if(swich==2){
		render->renderMat(sGrid->distanceLevelSet,2);
	}
	else if(swich==3){
		render->renderBoundary();
		render->renderGrid();
		render->renderVector2D(sGrid->u,sGrid->v); 
		}
	
	
	/*if(flag[0])
		render->renderBoundary();
	if(flag[1])
		render->renderGrid();
	if(flag[2])
		render->renderParticles();
	if(flag[3])
		render->renderSurfaceBoundary();
	if(flag[4])
		render->renderVector2D(sGrid->u,sGrid->v); */
	//if(flag[5])
		//render->renderMat(sGrid->distanceLevelSet,2);
		//render->renderMat(sGrid->p,2);
	/*if(flag[6])
		render->renderMat(sGrid->isFluidBoundary,1);
*/
	output1 = ' ';
	postDisplay();
}
void idleFun ( void )
{
	extern bool isPause ;
	if(!isPause){
	struct timeval tt1, tt2;
	static int it = 0;
	static struct timeval tt3, tt4;
	gettimeofday(&tt4, NULL);

	gettimeofday(&tt1, NULL);
	animate();
	gettimeofday(&tt2, NULL);
	int milliSeconds2 = (tt4.tv_sec - tt3.tv_sec) * 1000 + (tt4.tv_usec - tt3.tv_usec)/1000;
	cout<<"Rander Time Frame   "<<it-1<<" : "<<milliSeconds2<<"ms"<<endl;

	int milliSeconds = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
	cout<<"Iteration "<<it<<" : "<<milliSeconds<<"ms"<<endl<<endl;
	gettimeofday(&tt3, NULL);
	
	it++;
	glutSetWindow ( winId );
	glutPostRedisplay ( );
	}
}

void reshape (int w, int h)
{
   glutSetWindow(winId);
   glutReshapeWindow(w ,h);
   winSizeX = w;
   winSizeY = h;
}

void preDisplay()
{
	glViewport ( 0, 0, winSizeX, winSizeY);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho( 0, zoomFactor, 0, zoomFactor	, 0,1 ); //better to use ortho..
	glClearColor( 0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

void postDisplay()
{
   glutSwapBuffers();
}

void openGlutWindow(char* windowName)
{
   glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
   glutInitWindowSize (winSizeX, winSizeY);
   glutInitWindowPosition (570, 100);

   winId = glutCreateWindow (&windowName[2]);

   glClearColor( 0.0f, 0.0f, 0.0f, 1.0f);
   glClear(GL_COLOR_BUFFER_BIT);
   glutSwapBuffers();
   glClear(GL_COLOR_BUFFER_BIT);
   glutSwapBuffers();
   glutDisplayFunc(display);
//   glutSpecialFunc(&SpecialKeyPressed);
   glutKeyboardFunc(&KeyPressed);
   glutReshapeFunc(reshape);
   glutIdleFunc(idleFun);
}

