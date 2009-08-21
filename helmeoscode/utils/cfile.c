#include <stdio.h>


// Compile:
// g77 -c fortranfile.f
// g++ -o myprogram cfile.c fortranfile.o -lg2c
// icc -c cfile.c

// f2c: http://astro.berkeley.edu/~wright/f2c.html
// http://wwwcompass.cern.ch/compass/software/offline/software/fandc/fandc.html
// http://www.chiralcomp.com/support/mixing_f77_c_cpp/defcall.html

// http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
// http://www.physiology.wisc.edu/comp/docs/notes/not017.html

// Prototype for Fortran

extern "C" {  extern void event_(int*,float*,float*); }



// main C function
int main(void)
{

  int kflavour=3;
  float energy=90.0;
  float result;
  

  event_(&kflavour,&energy,&result);


  fprintf(stdout,"result=%21.15g\n",result);

  return(1);


}
