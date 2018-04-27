/* bisulfitefq-fragment-multirate.c

   Junfeng Liu, 2017
   
   This is an program to bisulfite the fastq file in unit of each fragment according to multiple methylation rate.

   gcc -o bisulfitefq-fragment-multirate bisulfitefq-fragment-multirate.c -lm
   
   bisulfitefq-fragment-multirate <fastq_1 file> <fastq_2 file> <methylation_rate_1> <methylation_rate_2> <methylation_rate_3> <tag>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


/*The function is used to bisulfite the fastq file.*/
int BisulfiteFq (char *fastq1f, char *fastq2f, char *methylation_rate_1, char *methylation_rate_2, char *methylation_rate_3);
int BisulfiteFq_noanti (char *fastq1f, char *fastq2f, char *methylation_rate_1, char *methylation_rate_2, char *methylation_rate_3);

/*random number generator.*/
void rndu (double *r, double p[]);

int main (int argc, char* argv[])
{
   int l;
   
   char fastq1[1000], fastq2[1000], methylation_rate_1[1000], methylation_rate_2[1000], methylation_rate_3[1000], tag[2];

   if(argc>6) {
   	strcpy(fastq1, argv[1]);
   	strcpy(fastq2, argv[2]);
   	strcpy(methylation_rate_1, argv[3]);
   	strcpy(methylation_rate_2, argv[4]);
   	strcpy(methylation_rate_3, argv[5]);
   	strcpy(tag, argv[6]);
   } else {
   	return -1;
   }

   l=BisulfiteFq(fastq1,fastq2,methylation_rate_1,methylation_rate_2,methylation_rate_3);
   
   if(tag[0]=='1') {
   	l=BisulfiteFq(fastq1,fastq2,methylation_rate_1,methylation_rate_2,methylation_rate_3);
   } else {
   	l=BisulfiteFq_noanti(fastq1,fastq2,methylation_rate_1,methylation_rate_2,methylation_rate_3);
   }
     
   return (l);
}

int BisulfiteFq (char *fastq1f, char *fastq2f, char *methylation_rate_1, char *methylation_rate_2, char *methylation_rate_3)
{
  /*l for 'for loop'*/
	int l;
	/*n for counting the read number;
	  m1 for the number of fragment;
	*/
	long n, m1;
	
	/*r is random number between 0 and 1;
	  s1,*r1 are the initial value for fastq;
	*/
	double r, s1, *r1;
	
	/**/
	double p1[2];
	
	/*ls1[] for the row of the fastq1 file;
	  ls2[] for the row of the fastq2 file;
	*/
	char ls1[5000], ls2[5000];
	
	/*te is tag for bisulfiting*/
	char te;
	
	FILE *ffastq1, *ffastq2, *fout1, *fout2;
	
	if((ffastq1=fopen(fastq1f, "r"))==NULL) {
   		return -1;
   }
   
  if((ffastq2=fopen(fastq2f, "r"))==NULL) {
   		return -1;
   }
   
	if((fout1=fopen("multirate_1.fastq", "w"))==NULL) {
   		return -1;
   }
   
  if((fout2=fopen("multirate_2.fastq", "w"))==NULL) {
   		return -1;
   }
  
  n=0;
  m1=0;
  
  s1=abs(2*(int)time(NULL)+1);
	r1=&s1;
  
	 /*read each row of the fastq file*/
	 while(fgets(ls1,5000,ffastq1) != NULL) {
	 	fgets(ls2,5000,ffastq2);
	 	if(fmod(n,4)==0) {
	 		te='0';
	 		if(m1>0) {
	 			s1=p1[1];
	 			r1=&s1;
	 			rndu(r1,p1);
	 		} else {
	 			rndu(r1,p1);
	 			m1=m1+1;
	 		}
	 		if(strstr(ls1,".1:")!=NULL) {
	 			if(p1[0]<atof(methylation_rate_1)) {
	 				te='1';
	 			}
	 		}
	 		if(strstr(ls1,".2:")!=NULL) {
	 			if(p1[0]<atof(methylation_rate_2)) {
	 				te='1';
	 			}
	 		}
	 		if(strstr(ls1,".3:")!=NULL) {
	 			if(p1[0]<atof(methylation_rate_3)) {
	 				te='1';
	 			}
	 		}
	 		if(strstr(ls1,"methylated_liu")!=NULL) {
	 			if(te=='1') {
	 				fprintf(fout1,"%s",ls1);
	 			} else {
	 				fprintf(fout1,"%c",'@');
	 				for(l=15;l<strlen(ls1);l++) {
	 					fprintf(fout1,"%c",ls1[l]);
	 				}
	 			}
	 		} else {
	 			if(te=='1') {
	 				fprintf(fout1,"%s","@methylated_liu");
	 				for(l=1;l<strlen(ls1);l++) {
	 					fprintf(fout1,"%c",ls1[l]);
	 				}
	 			} else {
	 				fprintf(fout1,"%s",ls1);
	 			}
	 		}
	 		if(strstr(ls2,"methylated_liu")!=NULL) {
	 			if(te=='1') {
	 				fprintf(fout2,"%s",ls2);
	 			} else {
	 				fprintf(fout2,"%c",'@');
	 				for(l=15;l<strlen(ls2);l++) {
	 					fprintf(fout2,"%c",ls2[l]);
	 				}
	 			}
	 		} else {
	 			if(te=='1') {
	 				fprintf(fout2,"%s","@methylated_liu");
	 				for(l=1;l<strlen(ls2);l++) {
	 					fprintf(fout2,"%c",ls2[l]);
	 				}
	 			} else {
	 				fprintf(fout2,"%s",ls2);
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

int BisulfiteFq_noanti (char *fastq1f, char *fastq2f, char *methylation_rate_1, char *methylation_rate_2, char *methylation_rate_3)
{
  /*l for 'for loop'*/
	int l;
	/*n for counting the read number;
	  m1 for the number of fragment;
	*/
	long n, m1;
	
	/*r is random number between 0 and 1;
	  s1,*r1 are the initial value for fastq;
	*/
	double r, s1, *r1;
	
	/**/
	double p1[2];
	
	/*ls1[] for the row of the fastq1 file;
	  ls2[] for the row of the fastq2 file;
	*/
	char ls1[5000], ls2[5000], methylation_rate[1000];
	
	/*te is tag for bisulfiting*/
	char te;
	
	FILE *ffastq1, *ffastq2, *fout1, *fout2;
	
	if((ffastq1=fopen(fastq1f, "r"))==NULL) {
   		return -1;
   }
   
  if((ffastq2=fopen(fastq2f, "r"))==NULL) {
   		return -1;
   }
   
	if((fout1=fopen("multirate_1.fastq", "w"))==NULL) {
   		return -1;
   }
   
  if((fout2=fopen("multirate_2.fastq", "w"))==NULL) {
   		return -1;
   }
  
  n=0;
  m1=0;
  
  s1=abs(2*(int)time(NULL)+1);
	r1=&s1;
  
	 /*read each row of the fastq file*/
	 while(fgets(ls1,5000,ffastq1) != NULL) {
	 	fgets(ls2,5000,ffastq2);
	 	if(fmod(n,4)==0) {
	 		memset(methylation_rate,0,sizeof(methylation_rate));
	 		strcpy(methylation_rate,"0.2");
	 		if(strstr(ls1,".1:")!=NULL) {
	 			memset(methylation_rate,0,sizeof(methylation_rate));
	 			strcpy(methylation_rate,methylation_rate_1);
	 		}	 		
	 		if(strstr(ls1,".2:")!=NULL) {
	 			memset(methylation_rate,0,sizeof(methylation_rate));
	 			strcpy(methylation_rate,methylation_rate_2);
	 		}	 
	 		if(strstr(ls1,".3:")!=NULL) {
	 			memset(methylation_rate,0,sizeof(methylation_rate));
	 			strcpy(methylation_rate,methylation_rate_3);
	 		}	 				
	 	}
	 	if(fmod(n,4)==1) {
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
	 		for(l=0;l<strlen(ls1);l++) {
	 			if(ls1[l]=='C'||ls1[l]=='c') {
	 				if(te=='1') {
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
	 				if(te=='1') {
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