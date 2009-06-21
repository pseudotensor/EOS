# copied off of HARM method

include maketail.kazeos.inc


all:	$(CMD) $(FINISH)



$(PREP):
	#
	( echo $(OBJSONEPOINT) > objsonepoint )
	# remove and touch (so exists) files to create
	( rm -rf objsonepoint.cmd.list ; rm -rf onepoint ; touch objsonepoint.cmd.list ; touch onepoint)
	# create list of prefixes
	( for fil in `cat objsonepoint` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> onepoint ; done )
	# create list of commands
	( for fil in `cat onepoint` ; do echo -e $$fil.o :\\n\\t$(FC) $(FFLAGS) -c $$fil.f $(LDFLAGS) >> objsonepoint.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.onepoint.list ; touch deps.onepoint.list )
	# create new dependencies file
	( for fil in `cat onepoint` ; do  echo $$fil.f ; $(FC) -MM $$fil.f >> deps.onepoint.list ; done )
	#
	
$(FINISH):
#	rm -rf objsonepoint objsmanypoints

$(CMD):	$(OBJSONEPOINT) Makefile
	$(FC)  $(FFLAGS) -o $(CMD) $(OBJSONEPOINT) $(LDFLAGS)


-include objsonepoint.cmd.list
-include deps.onepoint.list


clean:
	rm *.o

cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~


