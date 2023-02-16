#Makefile para compilar Programa de Proteinas
.SUFFIXES : .o .f90
OBJ1 =  Ziggurat.o main.o secuencia.o 
#OBJ1 =  Ziggurat.o main.o initialize.o energy.o updates.o secuencia.o 

exec   : ${OBJ1}  Makefile
	gfortran -o $@ ${OBJ1}

.f90.o :
	gfortran -c $<

clean :
	rm *.o *.mod  K*.xyz exec
