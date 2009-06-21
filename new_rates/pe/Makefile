mainprogram.exe: mainprogram.o pe.o
	ifc -cm -w90 -132 -o mainprogram.exe mainprogram.o pe.o

clean:
	rm -f *.o core* *~ 

mainprogram.o: mainprogram.f
	ifc -cm -w90 -132 -c   mainprogram.f

pe.o: pe.f
	ifc -cm -w90 -132 -c   pe.f

# end of makefile.


