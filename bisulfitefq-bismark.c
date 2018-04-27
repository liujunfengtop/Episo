/* bisulfitefq.c

   Junfeng Liu, 2015
   
   This is an program to bisulfite the fastq file in unit of each fragment.

   gcc -o bisulfitefq-fragment bisulfitefq-fragment.c -lm
   
   bisulfitefq-fragment <fastq_1 file> <fastq_2 file> <methylation_rate>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


/*The function is used to bisulfite the fastq file.*/
int BisulfiteFq (char *fastq1f, char *methylation_rate);

/*random number generator.*/
void rndu (double *r, double p[]);

int main (int argc, char* argv[])
{
   int l;
   
   char fastq1[1000], methylation_rate[1000];

   if(argc>2) {
   	strcpy(fastq1, argv[1]);
   	strcpy(methylation_rate, argv[2]);
   } else {
   	return -1;
   }

   l=BisulfiteFq(fastq1,methylation_rate);
     
   return (l);
}

int BisulfiteFq (char *fastq1f, char *methylation_rate)
{
  /*l for 'for loop'*/
	int l;
	/*n for counting the read number;
	  m1 for the number of fragment;
	*/
	long  m1;
	
	/*r is random number between 0 and 1;
	  s1,*r1 are the initial value for fastq;
	*/
	double r, s1, *r1;
	
	int i,j,read,k;
	
	/**/
	double p1[2];
	
	/*ls1[] for the row of the fastq1 file;
	  ls2[] for the row of the fastq2 file;
	*/
	char ls1[5000], ls2[5000], temp[5000];
	
	/*te is tag for bisulfiting*/
	char te;
	
	FILE *ffastq1,  *fout1;
	
	if((ffastq1=fopen(fastq1f, "r"))==NULL) {
   		return -1;
   }
   
   
	if((fout1=fopen("bisulfite_bismark_pe.txt", "w"))==NULL) {
   		return -1;
   }
   
   for(i=0;i<1;i++) {
   	fgets(ls1,5000,ffastq1);
   	fprintf(fout1,"%s",ls1);
   }
  
  m1=0;
  
  s1=abs(2*(int)time(NULL)+1);
	r1=&s1;
	
	while(fgets(ls1,5000,ffastq1) != NULL) {
   	  k=0;
   	  j=0;
   	  read=0;
   	  memset(temp,0,sizeof(temp));

   	  /*output the row of the read file(or the read1 file and the read2 file)*/
   	  while(read<15) {
   	  	temp[k]=ls1[j];
   	  	/*read each column of each row from the txt file*/
   	  	if(ls1[j]==' '||ls1[j]=='\t'||ls1[j]=='\n') {
   	  		/*read the first column of the row and record the seq-name for temp2*/
   	  		if((read<=6)||(read==8)||(read==9)||(read>10)) {
   	  			fprintf(fout1,"%s",temp);
   	  		}
   	  		if(read==7) {
   	  			te='0';
	 		      if(m1>0) {
	 			      s1=p1[1];
	 			      r1=&s1;
	 			      rndu(r1,p1);
	 		      } else {
	 			      rndu(r1,p1);
	 			      m1=m1+1;
	 		      }
	 		      if(p1[0]<atof(methylation_rate)) {
	 			      te='1';
	 		      }
	 		      for(l=0;l<=k;l++) {
	 		      	if(temp[l]=='.'||temp[l]=='\t') {
	 		      		fprintf(fout1,"%c",temp[l]);
	 		      	}
	 		      	if(temp[l]=='H') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'H');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'h');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='X') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'X');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'x');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='Z') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'Z');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'z');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='U') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'U');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'u');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='h') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'H');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'h');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='x') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'X');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'x');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='z') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'Z');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'z');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='u') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'U');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'u');
	 		      		}
	 		      	}
	 		      }
   	  		}
   	  		if(read==10) {
   	  			for(l=0;l<=k;l++) {
	 		      	if(temp[l]=='.'||temp[l]=='\t') {
	 		      		fprintf(fout1,"%c",temp[l]);
	 		      	}
	 		      	if(temp[l]=='H') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'H');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'h');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='X') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'X');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'x');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='Z') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'Z');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'z');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='U') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'U');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'u');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='h') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'H');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'h');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='x') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'X');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'x');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='z') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'Z');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'z');
	 		      		}
	 		      	}
	 		      	if(temp[l]=='u') {
	 		      		if(te=='1') {
	 		      			fprintf(fout1,"%c",'U');
	 		      		} else {
	 		      			fprintf(fout1,"%c",'u');
	 		      		}
	 		      	}
	 		      }
   	  		}
   	  		memset(temp,0,sizeof(temp));
   	  		read=read+1;
   	  		k=0;
   	  	} else {
   	  		k=k+1;
   	  	}
   	  	j=j+1;
   	  }
    }
	
  
	fclose(ffastq1);
	fclose(fout1);
	return 0;
}

void rndu (double *r, double p[])
{
	
	int i, m;
	double s, u, v;
	
	s=65536.0;
	u=2053.0;
	v=13849.0;
	for(i=0;i<1;i++) {
		*r=u*(*r)+v;
		m=(int)(*r/s);
		*r=*r-m*s;
		p[0]=*r/s;
		p[1]=*r;
	}
	return;
}