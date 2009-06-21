mainprogram.exe: mainprogram.o neint.o qne2.o
	ifc -cm -w90 -132 -o mainprogram.exe mainprogram.o neint.o qne2.o

clean:
	rm -f *.o core* *~ 

mainprogram.o: mainprogram.f
	ifc -cm -w90 -132 -c   mainprogram.f

neint.o: /home/kohri/ndaf/code3/ne/neint/neint.f
	ifc -cm -w90 -132 -c   /home/kohri/ndaf/code3/ne/neint/neint.f

qne2.o: qne2.f
	ifc -cm -w90 -132 -c   qne2.f

# end of makefile.


