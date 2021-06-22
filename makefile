PROG = $(wildcard prog_*.f90)
MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

FC      = gfortran
FCFLAGS = -fbacktrace -Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Wreal-q-constant -Wunused 
LPFLAGS = -llapack

PROGRAM = cmplx_lyap

default: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(PROG) $^ $(LPFLAGS) 

$(OBJS): %.o : %.f90 
	$(FC) $(FCFLAGS) -c $< $(LPFLAGS)

mod_proc.o mod_init.o: mod_prec.o
mod_lpck.o: mod_prec.o mod_proc.o


debug:
	@echo $(PROG)
	@echo $(MODS)
	@echo $(OBJS)

clean:
	rm $(PROGRAM) $(OBJS) $(patsubst %.o,%.mod,$(OBJS))


.PHONY = default debug clean