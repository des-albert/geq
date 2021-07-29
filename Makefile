FC = gfortran
LD = gfortran

FFLAGS = -cpp

SOURCES = solver_mod.f90 solver.f90 bndmat.f90 eqsil.f90 flux.f90 splnco.f90 intpol.f90 contour.f90 \
	startt.f90 compar.f90 xcur.f90 gelg.f90 saddle.f90 curf.f90 plotit.f90 topol.f90 smootl.f90

OBJECTS= $(SOURCES:.f90=.o)

.SUFFIXES: .f90 .o .mod 
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

solver: $(OBJECTS)
	$(LD) -o solver $(OBJECTS) -ldislin

clean: 
	rm -f *.o solver *.mod *.met
