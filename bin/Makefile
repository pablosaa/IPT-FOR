# Make file for ipt FORTRAN program
# Meteorologisches Institut - Univerität München
# Pablo Saavedra

FXX = ifort
MYFLAGS = -u -O3 -zero -132 -w95 -c

ipt : var_init.o ipt_knr_bbc2.o over_days.o absfunc.o deatte.o maths.o stp_fap.o make_K3_fap.o lwc_high_res.o
	$(FXX) -o $@ $^

lwc_high_res.o : lwc_high_res.f90 variables.mod
	$(FXX) $(MYFLAGS) $<

stp_fap.o : stp_fap.f90 variables.mod absorption.mod
	$(FXX) $(MYFLAGS) $<

make_K3_fap.o : make_K3_fap.f90 variables.mod absorption.mod
	$(FXX) $(MYFLAGS) $<

over_days.o : over_days.f90 variables.mod maths.mod
	$(FXX) $(MYFLAGS) $<

deatte.o : deatte.f90 variables.mod absorption.mod
	$(FXX) $(MYFLAGS) $<

maths.mod : maths.o
maths.o : maths.f90
	$(FXX) $(MYFLAGS) $<

ipt_knr_bbc2.o : ipt_knr_bbc2.f90 variables.mod
	$(FXX) $(MYFLAGS) $<

absorption.mod : absfunc.o absfunc.f90
absfunc.o : absfunc.f90 variables.mod
	$(FXX) $(MYFLAGS) $<

variables.mod : var_init.o var_init.f90
var_init.o : var_init.f90
	$(FXX) $(MYFLAGS) $<

# End of Make file
