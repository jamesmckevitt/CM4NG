FC=ifort

ACML_DIR=.
ACML_LIB=$(ACML_DIR)/libacml_mp.a -lrt

OMP=-qopenmp

SPEED=-O3

LARRAY=-mcmodel=medium

%.o: %.f90
	$(FC) $(OMP) $(SPEED) $(LARRAY) -c $<
%.o: %.f
	$(FC) $(OMP) $(SPEED) $(LARRAY) -c $<

OBJ=./machine.o ./mesh.o ./acceleration_profile_an.o ./convolution.o ./density_init_profile.o ./gravity_init_profile.o ./gravity_profile_an.o ./load.o ./main.o ./poisson_boundary.o ./projection.o ./save.o
# to create convolution.o please contact the authors (details in README.md)
OBJFILES=machine.o mesh.o acceleration_profile_an.o convolution.o density_init_profile.o gravity_init_profile.o gravity_profile_an.o load.o main.o poisson_boundary.o projection.o save.o

convolution.out: $(OBJ)
	$(FC) $(OBJ) $(OMP) $(SPEED) $(LARRAY) -o $@ $(ACML_LIB)

clean:
	rm -rf *.dat *.mod *.out *.o

run:
	./convolution.out