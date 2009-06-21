mainprogram.exe: mainprogram.o neint.o gammane.o subs.o
	ifc -cm -w90 -132 -o mainprogram.exe mainprogram.o neint.o gammane.o subs.o

clean:
	rm -f *.o core* *~ 

mainprogram.o: mainprogram.f
	ifc -cm -w90 -132 -c   mainprogram.f

neint.o: /home/kohri/ndaf/code3/ne/neint/neint.f
	ifc -cm -w90 -132 -c   /home/kohri/ndaf/code3/ne/neint/neint.f

gammane.o: gammane.f
	ifc -cm -w90 -132 -c   gammane.f

subs.o: subs.f
	ifc -cm -w90 -132 -c   subs.f

# end of makefile.


