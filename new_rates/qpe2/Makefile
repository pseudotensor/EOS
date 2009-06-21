mainprogram.exe: mainprogram.o neint.o qpe2.o
	ifc -cm -w90 -132 -o mainprogram.exe mainprogram.o neint.o qpe2.o

clean:
	rm -f *.o core* *~ 

mainprogram.o: mainprogram.f
	ifc -cm -w90 -132 -c   mainprogram.f

neint.o: /home/kohri/ndaf/code3/ne/neint/neint.f
	ifc -cm -w90 -132 -c   /home/kohri/ndaf/code3/ne/neint/neint.f

qpe2.o: qpe2.f
	ifc -cm -w90 -132 -c   qpe2.f

# end of makefile.


