#ifndef _KEYBOARD_H
#define _KEYBOARD_H
#include <stdlib.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
//#include "main.h"
//#include "commonData.h"
///*
//#define VEL 1
//#define DEN 2
//#define PAR 3
//#define PRE 4
//*/
//enum renderFlag{VEL=4, DEN, PAR, PRE, DUB, BOU, LES};//choose what to render..
//renderFlag output=LES; //default is velocity
//char output1='~';//extra.. to renderFlag.. it takes keys SHIFT + < ~ to ^ >  and ons offs alternatively different render
////enum ProjectionRoutines{ GSR,PCG,PCGBridson} ProjectFLAG ;
//bool idleOn = true;// to run the simulator on idle
//bool anyUpdation = true;
//
//double nearFactor = 0;
//double farFactor = zoomFactor;
//
//void animate();
//void init();
///*--checkIt--*/
//
////bool GSrFlag = false; //flag to decide the project routine to be used..
//bool AdvectFlag = true;
//bool AdvectFlagP = true;
//bool extrapolateFlag = true;
//
//bool AdvectFlagCell = false;
//bool AdvectFlagParticle = false;
//
//#define GSRflag 0
//#define PCGflag 1
//#define PCGBridsonflag 2
//
//
//
//int ProjectFLAG  = PCGBridsonflag;
//
//void initParticles();
//extern int winId;
//
//
enum fluidBody{ DAM_BREAK=1,DAM_CENTER=2,STATIC_BED=3,DOUBLE_DAM=4} fBT=DAM_CENTER ; //fluidBodyType
//
bool isPause = false;
int swich = 0;

void idleFun();
void initMain();
void initParticles();
void KeyPressed (unsigned char key, int x, int y)
{
if (key==27)
	exit(0);
switch( (key) )
{

	case '1' :
				fBT = DAM_BREAK;
				initMain();
				glutPostRedisplay();
				break;
	case 'i' :
				//initMain();
				//glutPostRedisplay();
				initParticles();
				break;
	case '2':
				fBT = DAM_CENTER;
				initMain();
				glutPostRedisplay();
				break;
case '3':
				fBT = STATIC_BED;
				initMain();
				glutPostRedisplay();
				break;
	
case '4':
				fBT = DOUBLE_DAM;
				initMain();
				glutPostRedisplay();
				break;
		
	case 'p' :
			isPause = !isPause;
			std::cout<<isPause<<std::endl;
			break;
	case 's':
			swich = (swich+1)%4;
			break;
	case 'n' :
				idleFun();
				break;
	/*case 'f':
			render3DFlag = (render3DFlag+1)	% 2;
			glutPostRedisplay();
			break;*/
//		default:
//			output1 = key;
//			anyUpdation = true;
//			break;
//		case '+':
//			idleOn = !idleOn;
//			anyUpdation = true;
//			break;
//		case 'n' :
//				animate();
//				//fluidSim->simulate();
//				//glutSetWindow ( winId );
//				//traceCellFlag = OFF;
////				glutPostRedisplay();
//				break;
//
//		case 'r' :
//				output = PRE;
////				glutPostRedisplay();
//				break;
//
//		case 't' :
//				output = VEL;
////				glutPostRedisplay();
//				break;
//		/*case 'y':
//				output = DEN	;
//				glutPostRedisplay();
//				break;
//		*/
//		case 'u':
//				output = DUB	;
////				glutPostRedisplay();
//				break;
//		case 'o':
//				output = BOU	;
////				glutPostRedisplay();
//				break;
//		case 'y':
//				output = LES	;
////				glutPostRedisplay();
//				break;
//		case 'p':
//				output = PAR	;
////				glutPostRedisplay();
//				break;
//		case 'm':
////				fBT= DAM_CENTER;
//				init();
//				flagC=2;
////				glutPostRedisplay();
//				break;
//		case 'l':
////				fBT= DAM_BREAK;
//				init();
//				flagC=2;
////				glutPostRedisplay();
//				break;
//		case 's':
////				fBT= fBT==STATIC_BED?DOUBLE_DAM:STATIC_BED;
//				init();
//				flagC=2;
//				break;
//		case 'g':
//				ProjectFLAG = (ProjectFLAG+1)%4;
//				{char str[][15] = {"GSR","PCG","PCGBridson","NoPressure"};
//				cout<<"\nProjectFlag : "<<str[ProjectFLAG]<<endl;}
//				break;
//
//		case 'q':
//				AdvectFlagCell = !AdvectFlagCell;
//				AdvectFlagParticle = false;
//				break;
//		case 'w':
//				AdvectFlagCell = false;
//				AdvectFlagParticle = !AdvectFlagParticle;
//				break;
//		case 'a':
//				AdvectFlag = !AdvectFlag;
//				{char str[][15] = {"RK2","Normal"};
//				cout<<"\nAdvectFlag : "<<str[(int)AdvectFlag]<<endl;}
//				break;
//		case 'e':
//				extrapolateFlag = !extrapolateFlag;
//				//glutPostRedisplay();
//				{char str[][15] = {"Off","On"};
//				cout<<"\nExtrapolation : "<<str[(int)extrapolateFlag]<<endl;}
//
//
//				break;
//		case 'I':
//				initParticles();
//				break;
//
//
//
//
}
}
//
//
/*
float xRotated, yRotated, zRotated;
void SpecialKeyPressed(int key, int x, int y)
{
 switch (key)
      {
       case GLUT_KEY_UP:
						yRotated += 11;
						break;
       case GLUT_KEY_DOWN:
						yRotated -= 11;
						break;
       case GLUT_KEY_LEFT:
  						   xRotated += 11;
						break;
       case GLUT_KEY_RIGHT:
						xRotated -= 11;
						break;
       case GLUT_KEY_PAGE_DOWN:
							zRotated += 11;
                    	    break;
       case GLUT_KEY_PAGE_UP:
						zRotated -= 11;
                	    break;
       default:
			break;
      }
     glutPostRedisplay();
}*/




#endif

