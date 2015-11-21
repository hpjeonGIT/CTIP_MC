#
# Makefile for SIMOX serial version
# 
#
.SUFFIXES: .o .f90
#
F90 =  ifort # gfortran
LD  = xild
#FLAGS =  -O3 -m64 -ipo -no-prec-div #-heap-arrays  
#FLAGS = -O3 -m64 -axSSE4.1,SSSE3 -xSSSE3 -mp1 -prec-div -pc80 -pad -ip
#FLAGS = -xHost -O3 -ipo -no-prec-div  -pc80 -pad -ip
FLAGS = -g -xHost -O3
F77FLAG = ${FLAGS}
#FLAGS =  -march=native -ffast-math -funroll-loops -O3 
#F77FLAG = -O3 -m64 -ipo -no-prec-div #-heap-arrays  # 
#F77FLAG = -march=native -ffast-math -funroll-loops -O3 -fomit-frame-pointer -I..
OBJ = datafmt.o main.o force.o vverlet.o ctip.o util.o eam.o \
      fft235.o kernel.o mfft235.o zfft3d.o
#LIB = -L/opt/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas
TARGET = simox_opt
${TARGET}:${OBJ}
	${F90} -o ${TARGET} ${OBJ} ${LIB} 

.f90.o:
	${F90} ${FLAGS} -c $< ${INC}
.f.o:
	${F90} ${F77FLAG} -c $<
clean:
	rm -rf *.mod *.o *.f90~ core ${TARGET}
