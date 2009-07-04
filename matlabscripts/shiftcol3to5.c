#include <stdio.h>

// icc -o shiftcol3to5 shiftcol3to5.c

// cp shiftcol3to5.c ~/research/eos/ ; rm eos.dat.shifted  ; icc -o shiftcol3to5 shiftcol3to5.c ; ./shiftcol3to5 &> out.txt

int main(void)
{

  FILE* in;
  FILE* out;
  //  int numlines=2000000; // from wc
  //  int numchars=1924000000;
  int numlines=5000000; // from wc
  long i,numchars;
  int numcharsperline=962; // from head -1 | wc
  int numcol=31; // from head -1 | wc
  char ch;
  int col;
  char buffer[300];
  int sizecol3;
  int j;
  int didout;



  numchars = numlines*numcharsperline;

  fprintf(stderr,"numchars=%ld\n",numchars);



  //in = fopen("test.dat","rt"); //test
  in = fopen("eos.dat","rt");
  out = fopen("eos.dat.shifted","wt");

  didout=0;
  col=1;
  for(i=1;i<=numchars;i++){

    ch = fgetc(in);
    fputc(ch,out);
    if(ch==' ') col++;
    if(ch=='\n'){ didout=0; col=1;}
    //    fprintf(stderr,"col=%d\n",col); // test

    if(col==3*2){
      j=0;
      while(1){
	buffer[j]=fgetc(in);
	if(j>0 && buffer[j-1]==' ' && buffer[j]==' '){ sizecol3=j+1; break;}
	//else{ fprintf(stderr,"buffer[%d]=%c\n",j,buffer[j]);  j++; }
	else{ j++; }
      }
      col+=2; // so this if won't repeat
      i+=sizecol3;
    }
    else if(col==6*2 && didout==0){
      // then write 3rd
      for(j=0;j<sizecol3;j++){
	fputc(buffer[j],out);
      }
      didout=1; // so this if won't repeat
    }

  }

  fclose(in);
  fclose(out);


  return(0);

}
