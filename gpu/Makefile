FC=nvfortran

ACML_DIR=.
ACML_LIB=$(ACML_DIR)/libacml_mp.a -lrt

INTEL_DIR=/opt/sw/vsc4/VSC/x86_64/glibc-2.17/skylake/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64 # specify the intel compiler library directory
INTEL_LIB=-L$(INTEL_DIR) -lifcore -lirc -limf -lsvml

OMP=-mp=multicore

SPEED=-fast -Mpreprocess

LARRAY=-mcmodel=medium

DEVICE=-gpu=cc80 -cuda -cudalib=cufft # specify destination gpu compute capability with the -gpu flag

%.o: %.cuf
	$(FC) $(OMP) $(SPEED) $(LARRAY) $(DEVICE) -c $<
%.o: %.f
	$(FC) $(OMP) $(SPEED) $(LARRAY) $(DEVICE) -c $<

OBJ=./precision_m.o ./cufft_m.o ./machine.o ./mesh.o ./density_init_profile.o ./gravity_init_profile.o ./load.o ./poisson_boundary.o ./gravity_profile_an.o ./acceleration_profile_an.o ./convolution.o ./save.o ./projection.o
# to create convolution.o please contact the authors (details in README.md)
OBJFILES=precision_m.o cufft_m.o machine.o mesh.o density_init_profile.o gravity_init_profile.o load.o poisson_boundary.o gravity_profile_an.o acceleration_profile_an.o convolution.o save.o projection.o

convolution.out: $(OBJ)
	$(FC) $(OBJ) $(SPEED) $(LARRAY) $(OMP) $(DEVICE) $(DEBUG) -o $@ $(ACML_LIB) $(INTEL_LIB) main.cuf

clean:
	rm -rf *.dat *.mod *.out *.o

run:
	./convolution.out