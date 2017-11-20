#$Id: makefile,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $

include makefile.inc

.SUFFIXES: .o .f90 .C 

.f90.o:
	$(f90) $(f90appflags) -c $< -o $@

.C.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@


# Input files
C++FILES = ABE.o Ansorg.o Block.o misc.o monitor.o Parallel.o Patch.o var.o\
           cgh.o bssnEScalar_class.o surface_integral.o ShellPatch.o

F90FILES = enforce_algebra.o fmisc.o initial_puncture.o prolongrestrict.o\
	   rungekutta4_rout.o bssnEScalar_rhs.o diff_new.o kodiss.o kodiss_sh.o\
	   lopsidediff.o sommerfeld_rout.o getnp4EScalar.o diff_new_sh.o\
	   shellfunctions.o bssn_rhs_ss.o Set_Rho_ADM.o

$(C++FILES): Block.h enforce_algebra.h fmisc.h initial_puncture.h microdef.h\
	     misc.h monitor.h MyList.h Parallel.h Patch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_class.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h

$(C++FILES) $(F90FILES): microdef.fh

ABE: $(C++FILES) $(F90FILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES) $(F90FILES) $(LDLIBS)

ReadData: ReadData.C microdef.fh derivatives.h
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $< $(F90FILES) $(LDLIBS)

ReadCT: ReadCT.C microdef.fh
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

Converge: Converge.C microdef.fh microdef.h
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

CompareData: CompareData.C microdef.fh
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

construct_bbox: construct_bbox.C microdef.fh misc.h misc.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $< misc.o $(F90FILES) $(LDLIBS)

selectplane: selectplane.C
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

ri2aphi: ri2aphi.C
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

interpdata: interpdata.C
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $<

All: ABE ReadData ReadCT CompareData construct_bbox selectplane ri2aphi\
     interpdata

clear:
	rm *.o ABE make.log ReadData ReadCT CompareData construct_bbox\
	selectplane ri2aphi interpdata -f
