

///////////////
//
// Program is C version of wrapper with MPI support so can call individual Fortran programs for each core or node
// Should be called by helmeoscode/script/chunkn.sh when using MPI mode.
//
// In MPI mode, use helmeoscode/scripts/chunkbunch.sh to have this code called
// Then do:
// 1) Compile helmeoscode as described in its Makefile
// 2) Ensure installdir is clean: rm -rf <installdir>/*
// 2) Copy required stuff: sh copyjonhelm.sh <installdir>
// 3) Ensure no running bad dead MPI binaries (e.g.): killall helmeosc
// 4) Enter <installdir> and do: sh chunkbunch.sh <systemtype> <datadir> <totalchunks> <chunklist>
// Where <chunklist> is not necessary (see chunkbunch.sh)
// Where first argument is 3 or 4 or 5 currently (see chunkbunch.sh).
// E.g. sh chunkbunch.sh 3 . 4
// 5) Then collate if did more than 1 chunk
// E.g. sh collatechunks.sh . 4
// 6) Link up to normal names if desired:
// ln -s eosother.final.dat eosother.dat
// ln -s eos.final.dat eos.dat
// ln -s eos.final.head eos.head
// ln -s eoscoulomb.final.dat eoscoulomb.dat
// ln -s eosazbar.final.dat eosazbar.dat
// 7) Process to eosnew and/or view in SM, etc.



// whether to allow use of MPI if desired
#define USEMPI (USINGMPI)
#define MAXCHUNKS 2000
#define MAXCHUNKSTRING (MAXCHUNKS*(1+3)) // roughly 3 digits with 1 space for MAXCHUNKS numbers.  Only valid for log10(MAXCHUNKS)<~3
#define MAXGENNAME 2000

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#include <sys/stat.h> // for system unix commands.  See also "info libc"
#include <sys/types.h> // for system unix commands.  See also "man -a mkdir"
#include <unistd.h> // see man -a chdir



#if(USEMPI)
#include <mpi.h>
#else
#define MPI_MAX_PROCESSOR_NAME 1000
#endif

/////////////////////
// Below are just ideas/comments while trying to figure out how to compile and link C and Fortran77 programs
//
// Compile:
// g77 -c fortranfile.f
// g++ -o myprogram cfile.c fortranfile.o -lg2c
// icc -c cfile.c

// f2c: http://astro.berkeley.edu/~wright/f2c.html
// http://wwwcompass.cern.ch/compass/software/offline/software/fandc/fandc.html
// http://www.chiralcomp.com/support/mixing_f77_c_cpp/defcall.html

// http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
// http://www.physiology.wisc.edu/comp/docs/notes/not017.html

// Prototype for Fortran function:
//extern "C" {  extern void teos_helm_(); }

extern void teos_helm_(void); // Fortran prototype


// C prototypes:
static int init_mpi(int *argc, char **argv[]);
static void myargs(int argc, char *argv[]);
static int get_chunklist(size_t strsize, char* chunkliststring, int *chunklist, int *numchunks);
static int print_chunklist(int numchunks,int *chunklist);
static int setup_teos_helm(int myid, int *chunklist, int totalchunks, char *jobprefix, char *cwdold, char *cwdnew);
static int finish_teos_helm(int myid, int *chunklist, int totalchunks, char *jobprefix, char *cwdold, char *cwdnew);

static void cpu0fprintf(FILE* fileptr, char *format, ...);
static void myffprintf(FILE* fileptr, char *format, ...);


int truenumprocs,myid;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int procnamelen;

int totalchunks;
char chunkliststring[MAXCHUNKSTRING];
int chunklist[MAXCHUNKSTRING];
int numchunks;

char DATADIR[MAXGENNAME];
char jobprefix[MAXGENNAME];

char cwdold[MAXGENNAME];
char cwdnew[MAXGENNAME];


/////////////////////
//
// main C function
//
/////////////////////
int main(int argc, char *argv[])
{

  int kflavour=3;
  float energy=90.0;
  float result;
  
  ///////////////////////////////////
  //
  myffprintf(stdout,"init_mpi Begin.\n");

  init_mpi(&argc,&argv);

  myffprintf(stdout,"init_mpi End: myid=%d.\n",myid);

  ///////////////////////////////////
  //
  myffprintf(stdout,"myargs Begin: myid=%d.\n",myid);

  myargs(argc,argv);

  myffprintf(stdout,"myargs End: myid=%d.\n",myid);



  ///////////////////////////////////
  //
  //  (i.e. runchunkn.sh is called as the mpi binary in Jon's scripts.)
  myffprintf(stdout,"C runchunkn.sh -like Begin: myid=%d.\n",myid);

  setup_teos_helm(myid,chunklist,totalchunks,jobprefix,cwdold,cwdnew);

  myffprintf(stdout,"C runchunkn.sh -like End: myid=%d.\n",myid);


  ///////////////////////////////////
  //
  // Note that only myid==0's Fortran call will necessarily show Fortran code output.  Rest of CPU's output may not be redirected (i.e. that's implementation dependent)
  myffprintf(stdout,"C EOS Begin: myid=%d.\n",myid);

  teos_helm_();

  myffprintf(stdout,"C EOS Done: myid=%d.\n",myid);

  finish_teos_helm(myid,chunklist,totalchunks,jobprefix,cwdold,cwdnew);


  myffprintf(stdout,"Done with jon_helm_programstart.c.\n");

  return(1);


}


// Note we pass pointers since MPI_init() modifies the values
static int init_mpi(int *argc, char **argv[])
{
  int ierr;

#if(USEMPI)
  ierr=MPI_Init(argc, argv);

  if(ierr!=0){
    myffprintf(stderr,"MPI Error during MPI_Init\n");
    exit(1);
  }
  
  MPI_Comm_size(MPI_COMM_WORLD, &truenumprocs); // WORLD total number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); // WORLD proc id
  MPI_Get_processor_name(processor_name, &procnamelen); // to ensure really on certain nodes

  myffprintf(stderr, "WORLD proc: %d of %d on %s\n", myid,truenumprocs,processor_name);

  myffprintf(stderr, "end: init_MPI\n");
  fflush(stderr);
#else
  truenumprocs=1;
  myid=0;
#endif
  
  return(0);

}



static void myargs(int argc, char *argv[])
{
  int argi,numargs,numextraargs;
  size_t strsize;

  numargs=4; // number of user arguments

  ////////////////
  //
  // Get arguments from command-line
  //
  ////////////////
  if(argc!=1+numargs){ // 1 is normal command-line argument
    myffprintf(stderr,"Not enough args given!\n");
    myffprintf(stderr,"Expected: <binaryname> \"string of chunk numbers separated by spaces.\" <totalchunks> <DATADIR> <jobprefix>\n");
    myffprintf(stderr,"E.g.: mpirun -np 4 ./helmeosc \"1 2 3 4\" 4 . eoschunk\n");
    exit(1);
  }
  else{
    // arguments should be same and in same order as for runchunkn.sh for easy calling in script that sets up the job
    // argci=0 would access command-line argument
    argi=1; // 1 is first true argument after command line and MPI extracted stuff
    // argi=1:
    myffprintf(stderr,"argv[%d]=%s\n",argi,argv[argi]);
    strcpy(chunkliststring,argv[argi]); argi++;
    strsize=strlen(chunkliststring);
    if(strsize>MAXCHUNKSTRING){
      myffprintf(stderr,"Increase MAXCHUNKSTRING or use malloc!\n");
      exit(1);
    }
    myffprintf(stderr,"strsize=%d chunkliststring=%s\n",(int)strsize,chunkliststring);
    // argi=2:
    myffprintf(stderr,"argv[%d]=%s\n",argi,argv[argi]);
    totalchunks = atoi(argv[argi]); argi++;
    // argi=3:
    myffprintf(stderr,"argv[%d]=%s\n",argi,argv[argi]);
    strcpy(DATADIR,argv[argi]); argi++;
    // argi=4:
    myffprintf(stderr,"argv[%d]=%s\n",argi,argv[argi]);
    strcpy(jobprefix,argv[argi]); argi++;



    ///////////////////////
    //
    // now get chunk list
    //
    ///////////////////////
    get_chunklist(strsize,chunkliststring,chunklist,&numchunks);
    print_chunklist(numchunks,chunklist);

    if(numchunks!=truenumprocs){
      myffprintf(stderr,"Must have numchunks=%d equal to truenumprocs=%d\n",numchunks,truenumprocs);
      myffprintf(stderr,"Required since cannot fork(), so each proc can only call 1 Fortran EOS call.\n");
      exit(1);
    }
    else{
      myffprintf(stderr,"Good chunk count: numchunks=%d\n",numchunks);
    }// end else if good numchunks
  }// end if good number of arguments



}

static int print_chunklist(int numchunks,int *chunklist)
{
  int i;

  myffprintf(stderr,"Got %d chunks\n",numchunks);
  myffprintf(stderr,"chunks are: ");
  for(i=0;i<numchunks;i++){
    myffprintf(stderr,"%d ",chunklist[i]);
  }
  myffprintf(stderr,"\n");

  return(0);
}


// fill chunklist[] with numbers from string
static int get_chunklist(size_t strsize, char* chunkliststring, int *chunklist, int *numchunks)
{
  char *nptr,*endptr;

  *numchunks=0; // initialize
  nptr=&chunkliststring[0]; // setup pointer to start reading string from

  while(1){
    chunklist[*numchunks]=(int)strtol(nptr,&endptr,10);

    if(strlen(nptr)==strlen(endptr)){
      // check if done with string
      // then not going anywhere anymore, so probably whitespace at end with no valid characters
      break;
    }
    else{

      if(chunklist[*numchunks]<=0){
	myffprintf(stderr,"Chunk number problem: %d\n",chunklist[*numchunks]);
	exit(1);
      }

      // iterate
      *numchunks = (*numchunks)+1;
      nptr=endptr;
      // DEBUG:
      //      myffprintf(stderr,"cl[%d]=%d %d %d : %s\n",*numchunks,chunklist[*numchunks-1],strlen(nptr),strlen(chunkliststring),nptr);
    }
  }

  return(0);
}




// do things like in runchunkn.sh script
static int setup_teos_helm(int myid, int *chunklist, int totalchunks, char *jobprefix, char *cwdold, char *cwdnew)
{
  int subchunk;
  int subjobnumber;
  char subjobname[MAXGENNAME];
  char subjobdir[MAXGENNAME];  
  int error;
  FILE *eoschunkfile;


  /////////////////////
  //
  // get original working directory
  //
  /////////////////////
  int sizecwdold=MAXGENNAME;
  if(getcwd(cwdold,sizecwdold)==NULL){
    myffprintf(stderr,"Could not get old working directory.\n");
    exit(1);
  }


  /////////////////////
  //
  // Change directory
  //
  /////////////////////
  // Chunk list should never be <1.  First chunk possible is 1 (i.e Fortran index type)
  subchunk=myid+1; // dose-out chunks by CPU id number.
  subjobnumber=chunklist[myid]; // access array with 0 as first element. (C index type.)
  sprintf(subjobname,"%sc%dtc%d",jobprefix,subjobnumber,totalchunks);
  sprintf(subjobdir,"%s/%s",DATADIR,subjobname);

  // assumes path exists.  Let fail if not, since then not setup properly using chunkbunch.sh script
  // new path for all future file accesses that don't give full path
  error=chdir(subjobdir);

  if(error!=0){
    myffprintf(stderr,"Failed to change to directory: %s\n",subjobdir);
    exit(1);
  }

  ////////////////////////
  //
  // setup eoschunk.dat
  //
  ////////////////////////
  eoschunkfile=fopen("eoschunk.dat","wt"); // new file, do not append
  if(eoschunkfile==NULL){
    myffprintf(stderr,"Failed to open eoschunk.dat file\n");
    exit(1);
  }

  myffprintf(eoschunkfile,"%d %d\n",subchunk,totalchunks);

  fclose(eoschunkfile);

  ////////
  //
  // Remove old EOS files
  //
  // remove() is just same as unlink() for files and rmdir() for directories, but directory must be empty to use remove() or rmdir() on directories
  //
  ///////
  remove("eos.head");
  remove("eosdetails.dat");
  remove("eosother.dat");
  remove("eos.dat");
  remove("eoscoulomb.dat");
  remove("eosazbar.dat");

  ////////////
  //
  // get new working directory
  //
  ///////////
  int sizecwdnew=MAXGENNAME;
  if(getcwd(cwdnew,sizecwdnew)==NULL){
    myffprintf(stderr,"Could not get new working directory.\n");
    exit(1);
  }


  // ready to call Fortran EOS code!


  return(0);

}



// do things like in runchunkn.sh script AFTER binary is called
static int finish_teos_helm(int myid, int *chunklist, int totalchunks, char *jobprefix, char *cwdold, char *cwdnew)
{

  // change back to old working directory (not necessary currently)
  chdir(cwdold);


  return(0);
}

// only myid==0 prints
static void cpu0fprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;

  if  (myid==0) {
    va_start (arglist, format);

    if(fileptr==NULL){
      fprintf(stderr,"tried to print to null file pointer: %s\n",format);
      fflush(stderr);
    }
    else{
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);
    }
    va_end(arglist);
  }
}

// fprintf but also flushes
static void myffprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;

  va_start (arglist, format);

  if(fileptr==NULL){
    fprintf(stderr,"tried to print to null file pointer: %s\n",format);
    fflush(stderr);
  }
  else{
    vfprintf (fileptr, format, arglist);
    fflush(fileptr);
  }
  va_end(arglist);
}


