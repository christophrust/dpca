ARPACK_SOURCES = $(wildcard *.f *.F90)
ARPACK_OBJECTS = $(ARPACK_SOURCES:.f=.o)
ARPACK_OBJECTS := $(ARPACK_OBJECTS:.F90=.o)

# Pass FFLAGS and FPICFLAGS from parent Makevars
ARPACK_FFLAGS = $(FFLAGS) $(FPICFLAGS)

libarpack.a: $(ARPACK_OBJECTS)
	$(AR) rcs $@ $^

%.o: %.f
	$(FC) $(ARPACK_FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(ARPACK_FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a

.PHONY: clean
