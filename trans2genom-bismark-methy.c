/* trans2genom-bismark-methy.c

   Junfeng Liu, 2016
   
   This is an program to compute methylation rate at each site according to methylation_genom file.

   gcc -o trans2genom-bismark-methy trans2genom-bismark-methy.c -lm
   
   trans2genom-bismark-methy <methylation_genom>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>

/*The function is used to compute methylation rate at each site according to methylation_genom file*/
int convert (char *summary_complex);


int main (int argc, char* argv[])
{
	 int l;
   
   char out_trans[1000], summary_complex[1000];

   if(argc>1) {
   	strcpy(summary_complex, argv[1]);
   } else {
   	return -1;
   }

   l=convert(summary_complex);
     
   return (l);
}

int convert (char *summary_complex)
{
	/* j for chracter loop about each row of the methylation_summary-complex file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;
	
	long n;
	
	double m;
	
	
	/*ls[] for each row of the methylation_summary-complex file;
	  temp[] for each column of the row;
	*/
	char ls[5000], temp[1000];
	
  /*
    pre_site_position[] is for the previous site position;
  */
  char pre_site_position[1000];
  
  char tag;
	
	FILE *fcomplex, *fout;
	
	if((fcomplex=fopen(summary_complex, "r"))==NULL) {
		return -1;
  }
     
	if((fout=fopen("methylation_genom-rate", "w"))==NULL) {
   	return -1;
  }
  
  
  
  /*pre_site_position[] initial value*/
  strcpy(pre_site_position,"presiteposition");
  
  n=0.0;
  
  m=0.0;
  
  tag='I'; //initial by ljf 2017-06-08;
  
  /*read each row of the methylation_summary-complex file*/
	while(fgets(ls,5000,fcomplex) != NULL) {
		k=0;
	 	j=0;
	 	read=0;
	 	memset(temp,0,sizeof(temp));
	 	while(read<5) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==3) {
	 				if(strncmp(temp,pre_site_position,k)!=0) {
	 					if(strcmp(pre_site_position,"presiteposition")==0) {
	 						memset(pre_site_position,0,sizeof(pre_site_position));
	 						for(l=0;l<k;l++) {
	 							pre_site_position[l]=temp[l];
	 						}
	 					} else {
	 						fprintf(fout,"%s\t%ld\t%f\t%c\n",pre_site_position,n,m/n,tag); //add tag by ljf 2017-06-07; change %.2f to %f by ljf 2017-11-21;
	 						n=0;
	 						m=0.0;
	 						tag='I'; //add tag by ljf 2017-06-07;
	 						memset(pre_site_position,0,sizeof(pre_site_position));
	 						for(l=0;l<k;l++) {
	 							pre_site_position[l]=temp[l];
	 						}
	 					}
	 				}
	 				n=n+1;
	 			}
	 			if(read==4) {
	 				if(temp[0]=='X'||temp[0]=='H'||temp[0]=='Z') {
	 					m=m+1;
	 					tag=temp[0]; //add tag by ljf 2017-06-07;
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
	fprintf(fout,"%s\t%ld\t%f\t%c\n",pre_site_position,n,m/n,tag); //add tag by ljf 2017-06-08; change %.2f to %f by ljf 2017-11-21;
	fclose(fcomplex);
	fclose(fout);
	return 0;
}

