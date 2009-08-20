c..
c..declaration for pipelining the eos routines

c..maximum length of the row vector
      integer   nrowmax
      parameter (nrowmax = 40000)

c     Extra storage when want to manipulate _row information
      integer   nrowextra
      parameter (nrowextra = 3)


c..maximum number of isotopes
      integer   irowmax
      parameter (irowmax = 30)


c..maximum number of ionization stages
      integer   jstagemax
      parameter (jstagemax = 30)
