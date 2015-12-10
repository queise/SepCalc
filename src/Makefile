COMP   = gfortran
FLAGS  = -ffixed-line-length-150

all: mods exec

mods:
	@$(COMP) $(FLAGS) -c -o SepFunct.o SepFunct.f
	@$(COMP) $(FLAGS) -c -o SepOut.o SepOut.f

exec:
	@$(COMP) $(FLAGS) SepCalc.f SepFunct.f SepOut.f -o SepCalc
