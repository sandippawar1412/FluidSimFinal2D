CC=g++

GL_INC_DIR=/usr/include
GL_LIB_DIR=/usr/lib

GL_LIBS=-L$(GL_LIB_DIR) -lglut -lGLU -lGL -lm -lpthread

CFLAGS = -Wall 
OPTFLAG = -O3
# -dNDEBUG

.PHONY:clean

OBJ = main.o GridStag.o  FluidSim.o  Particles.o Printer.o  Renderer.o
LIBS = $(GL_LIBS)

all: $(OBJ) liquid2D

$(OBJ): %.o : %.cpp
	$(CC)  -I./ -c $(CFLAGS) $(OPTFLAG) $< -o $@ -fopenmp

#-Wno-write-strings : to ignore warning related to deprecated conversion from const char* to char*
liquid2D: $(OBJ) 
	$(CC)    $(CFLAGS) $(OBJ) $(LIBS) $(OPTFLAG) -o $@ -fopenmp

clean:
	rm -f  ./*~ ./core $(OBJ) liquid2D

#DEPENDENCIES- DO NOT DELETE
GridStag.o : GridStag.h ParameterFLAGS.hpp
Renderer.o : Renderer.h
FluidSim.o : FluidSim.h ./pcgsolver/*.h ParameterFLAGS.hpp 
main.o : main.h commonData.h ParameterFLAGS.hpp keyboard.h
Printer.o : Printer.h
Particles.o : Particles.h


