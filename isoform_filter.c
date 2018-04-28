/* isoform_filter.c

   Junfeng Liu, 2015
   
   This is an program to filter the isoform file in order to drop the isofrom whose FPKM is less than the assign number.

   gcc -o isoform_filter isoform_filter.c
   
   isoform_filter <isforms.fpkm_tracking> <isforms_methylation.fpkm_tracking> <FPKM_number>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


/*The function is used to filter the isoform file in order to drop the isofrom whose FPKM is less than the assign number*/
int filter (char *isoform, char *isoform_methylation, char *FPKM_number_1, char *FPKM_number_2);

int main (int argc, char* argv[])
{
   int l;
   
   char isoform[1000], isoform_methylation[1000], FPKM_number_1[1000], FPKM_number_2[1000];

   if(argc>4) {
   	strcpy(isoform, argv[1]);
   	strcpy(isoform_methylation, argv[2]);
   	strcpy(FPKM_number_1, argv[3]);
   	strcpy(FPKM_number_2, argv[4]);
   } else {
   	return -1;
   }

   
   l=filter(isoform,isoform_methylation,FPKM_number_1,FPKM_number_2);
   
   
     
   return (l);
}

int filter (char *isoform, char *isoform_methylation, char *FPKM_number_1, char *FPKM_number_2)
{
	/* j for chracter loop about each row of the isoform file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;
	
	/* j1 for chracter loop about each row of the isoform_methylation file;
     read1 for column loop about each row;
     k1 for character loop about each column;
  */
	int k1, j1, read1;
	
	
	/* n for fpkm from the isoform file;
	*/
	double n;
	
	/* m for fpkm from the isoform_methylation file;
	*/
	double m;
	
	/*ls[] for each row of the isoform file;
	  ls1[] for each row of the isoform_methylation file;
	  temp[] for each column of the row from the isoform file;
	  temp1[] for each column of the row from the isoform_methylation file;
	*/
	char ls[5000], ls1[5000], temp[1000], temp1[1000];
	
	FILE *fisoform, *fisoform_methylation, *fout, *fout1;
	
	if((fisoform=fopen(isoform, "r"))==NULL) {
   		return -1;
   }
   
  if((fisoform_methylation=fopen(isoform_methylation, "r"))==NULL) {
   		return -1;
   }
  
	
	if((fout=fopen("isoform_filter_out", "w"))==NULL) {
   		return -1;
   }

   
   if((fout1=fopen("isoform_methylation_filter_out", "w"))==NULL) {
   		return -1;
   }
   
  
	/*skipping the first row*/
   for(l=0;l<1;l++) {
   	fgets(ls,5000,fisoform);
   	fgets(ls1,5000,fisoform_methylation);
   	fprintf(fout,"%s",ls);
   	fprintf(fout1,"%s",ls1);
   }
   
   /*filter the isoform file*/
   while(fgets(ls,5000,fisoform) != NULL) {
   	k=0;
	 	j=0;
	 	read=0;
	 	n=0.0;
	 	k1=0;
	 	j1=0;
	 	read1=0;
	 	m=0.0;
	 	memset(temp,0,sizeof(temp));
	 	memset(temp1,0,sizeof(temp1));
	 	while(read<10) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==9) {
	 				n=atof(temp);
	 			}
	 			memset(temp,0,sizeof(temp));
	 			read=read+1;
	 			k=0;
	 		} else {
	 			k=k+1;
	 		}
	 		j=j+1;
	 	}
	  fgets(ls1,5000,fisoform_methylation);
	 	while(read1<10) {
	 		temp1[k1]=ls1[j1];
	 		if(ls1[j1]==' '||ls1[j1]=='\t'||ls1[j1]=='\n') {
	 			if(read1==9) {
	 				m=atof(temp1);
	 				if((n>=atof(FPKM_number_1))&&(m>=atof(FPKM_number_2))) {
	 					fprintf(fout,"%s",ls);
   	        fprintf(fout1,"%s",ls1);
	 				}
	 			}
	 			memset(temp1,0,sizeof(temp1));
	 			read1=read1+1;
	 			k1=0;
	 		} else {
	 			k1=k1+1;
	 		}
	 		j1=j1+1;
	 	}
   }
	 

	 

	fclose(fisoform);
	fclose(fisoform_methylation);
	fclose(fout);
	fclose(fout1);
	return 0;
}
