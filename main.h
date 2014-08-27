/*
 * main.h
 *
 *  Created on: 23-Jul-2014
 *      Author: srp
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "Printer.h"
#include "GridStag.h"
#include "FluidSim.h"

#include "Renderer.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
unsigned int winSizeX = 600;
unsigned int winSizeY = 600;
int winId;
Renderer* render = new Renderer;

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;



double timestep = 0.01;

GridStag* sGrid = new GridStag;
Printer * print = new Printer;
FluidSim* fluidSim = new FluidSim;


void initMain ( void ) ;
void animate () ;




#endif /* MAIN_H_ */
