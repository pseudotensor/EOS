mainprogram.exe: mainprogram.o ne.o
	ifc -cm -w90 -132 -o mainprogram.exe mainprogram.o ne.o

clean:
	rm -f *.o core* *~ 

mainprogram.o: mainprogram.f
	ifc -cm -w90 -132 -c   mainprogram.f

ne.o: ne.f
	ifc -cm -w90 -132 -c   ne.f

# end of makefile.


