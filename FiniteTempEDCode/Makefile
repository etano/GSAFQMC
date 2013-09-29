SHELL = /bin/sh
OS = Linux

Linux_LIB = -lblas -llapack 
Linux_MYLIB = 
Linux_CLIB =  -lm
Linux_FCLIB = 	
OSF1_LIB = -ldxml
OSF1_MYLIB = 
OSF1_CLIB =  -lm
OSF1_FCLIB = 	
IRIX_LIB = -lcomplib.sgimath
IRIX_MYLIB = 
IRIX_CLIB =  -lfastm -lm
IRIX_FCLIB = 	
# -lI77 -lF77 -lU77 -lisam -lm
# -lf2c

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
Linux_OPTFLAGS = -O3 -ffast-math
Linux_CFLAGS = -Wall $(${OS}_OPTFLAGS) 
Linux_FFLAGS = $(${OS}_OPTFLAGS)
OSF1_OPTFLAGS = -O5
OSF1_CFLAGS = $(${OS}_OPTFLAGS) 
OSF1_FFLAGS = $(${OS}_OPTFLAGS)

IRIX_OPTFLAGS = -O3 -mips4 -64 -r12000 -OPT:IEEE_arithmetic=3
IRIX_CFLAGS = $(${OS}_OPTFLAGS) 
IRIX_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = mix

#compiler ...
Linux_CC = cc
Linux_FF = gfortran
Linux_LD = cc

OSF1_CC = cc
OSF1_FF = f77
OSF1_LD = cc

IRIX_CC = cc
IRIX_FF = f90
IRIX_LD = cc

OBJECTS = \
	main.o matrix.o zero.o energy.o init.o \
        neighbors.o state.o print.o hubbard.o  

# compilation ...

# (use gmake for IRIX).
.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) $(${OS}_CFLAGS) -c  $*.c  
 
$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) $(${OS}_OPTFLAGS) -o $(MAINNAM).x $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
