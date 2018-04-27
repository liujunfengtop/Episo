/* bisulfitefq.c

   Junfeng Liu, 2015
   
   This is an program to bisulfite the fastq file.

   gcc -o bisulfitefq bisulfitefq.c tools.c -lm
   
   bisulfitefq <fastq_1 file> <fastq_2 file> <methylation_rate>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


/*The function is used to bisulfite the fastq file.*/
int BisulfiteFq (char *fastq1f, char *fastq2f, char *methylation_rate);

/*random number generator.*/
void rndu (double *r, double p[]);

int main (int argc, char* argv[])
{
   int l;
   
   char fastq1[1000], fastq2[1000], methylation_rate[1000];

   if(argc>3) {
   	strcpy(fastq1, argv[1]);
   	strcpy(fastq2, argv[2]);
   	strcpy(methylation_rate, argv[3]);
   } else {
   	return -1;
   }

   l=BisulfiteFq(fastq1,fastq2,methylation_rate);
     
   return (l);
}

int BisulfiteFq (char *fastq1f, char *fastq2f, char *methylation_rate)
{
  /*l for 'for loop'*/
	int l;
	/*n for counting the read number;
	  m1 for the number of C;
	  m2 for the number of G
	*/
	long n, m1, m2;
	
	/*r is random number between 0 and 1;
	  s1,*r1 are the initial value for fastq1;
	  s2,*r2 are the initial value for fastq2;
	*/
	double r, s1, s2, *r1, *r2;
	
	/**/
	double p1[2], p2[2];
	
	/*ls1[] for the row of the fastq1 file;
	  ls2[] for the row of the fastq2 file;
	*/
	char ls1[5000], ls2[5000];
	
	FILE *ffastq1, *ffastq2, *fout1, *fout2;
	
	if((ffastq1=fopen(fastq1f, "r"))==NULL) {
   		return -1;
   }
   
  if((ffastq2=fopen(fastq2f, "r"))==NULL) {
   		return -1;
   }
   
	if((fout1=fopen("FluxSim_bisulfite_1.fastq", "w"))==NULL) {
   		return -1;
   }
   
  if((fout2=fopen("FluxSim_bisulfite_2.fastq", "w"))==NULL) {
   		return -1;
   }
  
  n=0;
  m1=0;
  m2=0;
  
  s1=abs(2*(int)time(NULL)+1);
  s2=s1+1.0;
	r1=&s1;
	r2=&s2;
  
	 /*read each row of the fastq file*/
	 while(fgets(ls1,5000,ffastq1) != NULL) {
	 	fgets(ls2,5000,ffastq2);
	 	if(fmod(n,4)==1) {
	 		for(l=0;l<strlen(ls1);l++) {
	 			if(ls1[l]=='C'||ls1[l]=='c') {
	 				if(m1>0) {
	 					s1=p1[1];
	 					r1=&s1;
	 					rndu(r1,p1);
	 				} else {
	 					rndu(r1,p1);
	 					m1=m1+1;
	 				}
	 				if(p1[0]<atof(methylation_rate)) {
	 					fprintf(fout1,"%c",ls1[l]);
	 				} else {
	 					fprintf(fout1,"%c",'T');
	 				}
	 			} else {
	 				fprintf(fout1,"%c",ls1[l]);
	 			}
	 		}
	 		for(l=0;l<strlen(ls2);l++) {
	 			if(ls2[l]=='G'||ls2[l]=='g') {
	 				if(m2>0) {
	 					s2=p2[1];
	 					r2=&s2;
	 					rndu(r2,p2);
	 				} else {
	 					rndu(r2,p2);
	 					m2=m2+1;
	 				}
	 				if(p2[0]<atof(methylation_rate)) {
	 					fprintf(fout2,"%c",ls2[l]);
	 				} else {
	 					fprintf(fout2,"%c",'A');
	 				}
	 			} else {
	 				fprintf(fout2,"%c",ls2[l]);
	 			}
	 		}
	 	} else {
	 		fprintf(fout1,"%s",ls1);
	 		fprintf(fout2,"%s",ls2);
	 	}
	 	n=n+1;
	}
	fclose(ffastq1);
	fclose(ffastq2);
	fclose(fout1);
	fclose(fout2);
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