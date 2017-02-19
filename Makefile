# Makefile

EXE = run/heat_2d

OBJDIR = ./obj

default: $(OBJDIR) makefile.inc $(EXE)

makefile.inc:
	@echo "###############################################"
	@echo "# Before compiling mumps, you should have     #"
	@echo "# an appropriate file Makefile.inc available. #"
	@echo "###############################################"
	@exit 1

include makefile.inc


OBJS = obj/mod_precision.o \
       obj/mod_cartesian_grid.o \
       obj/mod_csr_matrix.o \
       obj/mod_utils.o \
       obj/mod_heat_2d.o \
       obj/radau5.o \
       obj/decsol.o \
       obj/dc_decsol.o \
       obj/mod_radau.o \
       obj/rock4.o \
       obj/rho.o \
       obj/mod_rock.o \
       obj/mod_forward_euler.o \
       obj/mod_backward_euler.o \
       obj/mod_integration.o \
       obj/heat_2d_main.o

$(OBJDIR):
	@echo "Creation du repertoire $@"
	@mkdir -p $@
	@echo

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

obj/%.o : src/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/%.o : src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@


#dependencies

obj/mod_utils.o: obj/mod_precision.o obj/mod_cartesian_grid.o

obj/mod_heat_2d.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_csr_matrix.o

obj/mod_integration : obj/mod_precision.o obj/mod_radau.o obj/mod_rock.o obj/mod_forward_euler.o obj/mod_backward_euler.o obj/mod_heat_2d

obj/mod_radau.o: obj/mod_precision.o 

obj/mod_backward_euler.o: obj/mod_precision.o obj/mod_csr_matrix.o

obj/mod_forward_euler.o: obj/mod_precision.o 

obj/heat_2d_main.o: obj/mod_precision.o obj/mod_cartesian_grid.o obj/mod_utils.o obj/mod_integration.o obj/mod_heat_2d.o

clean:
	rm -f $(OBJS) obj/*.mod $(EXE)
