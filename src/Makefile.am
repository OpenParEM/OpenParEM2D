bin_PROGRAMS = OpenParEM2D builder waveguide simplify process
OpenParEM2D_SOURCES = convergence.cpp fem2D.cpp fieldPoints.cpp findPoints.cpp modes.cpp OpenParEM2D.cpp results.cpp eigensolve.c Hsolve.c project.c 

builder_SOURCES = builder.cpp 
waveguide_SOURCES = waveguide.c
simplify_SOURCES = simplify.cpp project.c
process_SOURCES = process.cpp project.c

install:
	rm -f ../bin/OpenParEM2D
	cp OpenParEM2D ../bin/OpenParEM2D
	rm -f $(HOME)/bin/OpenParEM2D
	cp OpenParEM2D $(HOME)/bin/OpenParEM2D

	rm -f ../bin/process
	cp process ../bin/process
	rm -f $(HOME)/bin/process
	cp process $(HOME)/bin/process

	rm -f ../bin/simplify
	cp simplify ../bin/simplify
	rm -f $(HOME)/bin/simplify
	cp simplify $(HOME)/bin/simplify

	rm -f ../bin/builder
	cp builder ../bin/builder
	rm -f $(HOME)/bin/builder
	cp builder $(HOME)/bin/builder

	rm -f $(HOME)/bin/process.sh
	cp ../scripts/process.sh $(HOME)/bin/process.sh

	rm -f $(HOME)/bin/regression.sh
	cp ../scripts/regression.sh $(HOME)/bin/regression.sh

	rm -f $(HOME)/bin/checkBuilder.sh
	cp ../scripts/checkBuilder.sh $(HOME)/bin/checkBuilder.sh

	rm -f $(HOME)/.FreeCAD/Macro/OpenParEM2D_save.py
	cp ../scripts/OpenParEM2D_save.py $(HOME)/.FreeCAD/Macro/OpenParEM2D_save.py

	rm -f $(HOME)/.config/ParaView/Macros/field_plot.py
	cp ../scripts/field_plot.py $(HOME)/.config/ParaView/Macros/field_plot.py

