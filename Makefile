
FORTRAN=ifort

FFLAGS= -mcmodel=large -shared-intel -heap-arrays  -132

A_OUT=apex_prog

OBJECTS=\
	generate_apex_coordinates_lowres_v2.o \
	divve.o \
	apxntrpb4lf.o \
	apex.o \
	ggrid.o \
	magfld.o \
        calc_apex_params_2d_2_new2.o
 
$(A_OUT): $(OBJECTS)
	$(FORTRAN) -mcmodel=large -shared-intel -heap-arrays -o  $(A_OUT) $(OBJECTS) $(LFLAGS) 

.f.o:
	$(FORTRAN) $(FFLAGS) -c $<

calc_apex_params_2d_2_new2.o: calc_apex_params_2d_2_new2.f90
	$(FORTRAN) $(FFLAGS) -c calc_apex_params_2d_2_new2.f90

clean:
	- $(RM) $(A_OUT) *.o *.mod *.MOD *.lst *.a *.x
