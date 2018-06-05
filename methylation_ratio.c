/* methylation_ratio.c

   Junfeng Liu, 2018
   
   This is an program to compute methylation ratio on each transcript. 
   The program is the update verison on base of the version from 2015.
   The version of 2015 did not consider the effective length of transcript when computing the m5c level.

   gcc -o methylation_ratio methylation_ratio.c -lm
   
   methylation_ratio <isoform_filter> <isoform_methylation_filter> <total_number> <methylation_number>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


/*The function is used to compute the methylation ratio according to the file isoform and the file isoform_methylation*/
int ratio (char *isoform_filter, char *isoform_methylation_filter, char *total_number, char *methylation_number);

int main (int argc, char* argv[])
{
   int l;
   
   char isoform_filter[1000], isoform_methylation_filter[1000], total_number[1000], methylation_number[1000];

   if(argc>4) {
   	strcpy(isoform_filter, argv[1]);
   	strcpy(isoform_methylation_filter, argv[2]);
   	strcpy(total_number, argv[3]);
   	strcpy(methylation_number, argv[4]);
   } else {
   	return -1;
   }

   l=ratio(isoform_filter,isoform_methylation_filter,total_number,methylation_number);
     
   return (l);
}

int ratio (char *isoform_filter, char *isoform_methylation_filter, char *total_number, char *methylation_number)
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
	
	
	/* n1 for the fpkm value of the isoform from the isoform file;
	   n2 for the low fpkm of the isoform from the isoform file;
	   n3 for the high fpkm of the isoform from the isoform file;
	*/
	double n1, n2, n3;
	
	/* m1 for the fpkm value of the isoform from the isoform_methylation file;
	   m2 for the low fpkm of the isoform from the isoform_methylation file;
	   m3 for the high fpkm of the isoform from the isoform_methylation file;
	*/
	double m1, m2, m3;
	
	/*
	  ratio_estimated for the estimated value about methylation ratio;
	  ratio_estimated_lo for confidence interval low value about methylation ratio;
	  ratio_estimated_hi for confidence interval hi value about methylation ratio;
	*/
	double ratio_estimated, ratio_estimated_lo, ratio_estimated_hi;
	/*ls[] for each row of the isoform file;
	  ls1[] for each row of the isoform_methylation file;
	  temp[] for each column of the row from the isoform file;
	  temp1[] for each column of the row from the isoform_methylation file;
	  temp2[] for the transcript name;
	  temp3[] for the class_code;
	  temp4[] for nearest_ref_id;
	*/
	char ls[5000], ls1[5000], temp[1000], temp1[1000], temp2[1000], temp3[1000], temp4[1000];
	
	FILE *fisoform_filter, *fisoform_methylation_filter, *fout;
	
   
  if((fisoform_filter=fopen(isoform_filter, "r"))==NULL) {
   		return -1;
   }
   
  if((fisoform_methylation_filter=fopen(isoform_methylation_filter, "r"))==NULL) {
   		return -1;
   }
  
	if((fout=fopen("ratio_out", "w"))==NULL) {
   		return -1;
   }
   
  fprintf(fout,"%s\t%s\t%s\t%s\t%s\t%s\n","tracking_id","class_code","nearest_ref_id","ratio_estimated","ratio_estimated_lo","ratio_estimated_hi");
  
   
	 
	 /*skipping the first row*/
   for(l=0;l<1;l++) {
   	fgets(ls,5000,fisoform_filter);
   	fgets(ls1,5000,fisoform_methylation_filter);
   }
   
   /*compute ratio_estimated, ratio_estimated_lo and ratio_estimated_hi*/
   while(fgets(ls,5000,fisoform_filter) != NULL) {
   	k=0;
	 	j=0;
	 	read=0;
	 	n1=0.0;
	 	n2=0.0;
	 	n3=0.0;
	 	k1=0;
	 	j1=0;
	 	read1=0;
	 	m1=0.0;
	 	m2=0.0;
	 	m3=0.0;
	 	ratio_estimated=0.0;
	 	ratio_estimated_lo=0.0;
	 	ratio_estimated_hi=0.0;
	 	memset(temp,0,sizeof(temp));
	 	memset(temp1,0,sizeof(temp1));
	 	memset(temp2,0,sizeof(temp2));
	 	memset(temp3,0,sizeof(temp3));
	 	memset(temp4,0,sizeof(temp4));
	 	while(read<12) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==0) {
	 				for(l=0;l<k;l++) {
	 					temp2[l]=temp[l];
	 				}
	 			}
	 			if(read==3) {
	 				for(l=0;l<k;l++) {
	 					temp3[l]=temp[l];
	 				}
	 			}
	 			if(read==6) {
	 				for(l=0;l<k;l++) {
	 					temp4[l]=temp[l];
	 				}
	 			}
	 			if(read==9) {
	 				n1=atof(temp);
	 			}
	 			if(read==10) {
	 				n2=atof(temp);
	 			}
	 			if(read==11) {
	 				n3=atof(temp);
	 			}
	 			memset(temp,0,sizeof(temp));
	 			read=read+1;
	 			k=0;
	 		} else {
	 			k=k+1;
	 		}
	 		j=j+1;
	 	}
	  fgets(ls1,5000,fisoform_methylation_filter);
	 	while(read1<12) {
	 		temp1[k1]=ls1[j1];
	 		if(ls1[j1]==' '||ls1[j1]=='\t'||ls1[j1]=='\n') {
	 			if(read1==9) {
	 				m1=atof(temp1);
	 			}
	 			if(read1==10) {
	 				m2=atof(temp1);
	 			}
	 			if(read1==11) {
	 				m3=atof(temp1);
	 				if(n1==0.0||pow(n1,2)==0.0||pow(n1,4)==0.0) {
	 					fprintf(fout,"%s\t%s\t%s\t%f\t%f\t%f\n",temp2,temp3,temp4,0.0,0.0,0.0);
	 				} else {
	 					ratio_estimated=(atof(methylation_number)/atof(total_number))*(m1/n1+(m1*pow((n3-n1)/1.96,2))/pow(n1,3));
	 					if(ratio_estimated>1.0) {
	 						ratio_estimated=1.0;
	 					}
	 					ratio_estimated_lo=ratio_estimated-1.96*(atof(methylation_number)/atof(total_number))*sqrt(pow((m3-m1)/1.96,2)/pow(n1,2)+pow(m1,2)*pow((n3-n1)/1.96,2)/pow(n1,4));
	 					if(ratio_estimated_lo<0.0) {
	 						ratio_estimated_lo=0.0;
	 					}
	 					ratio_estimated_hi=ratio_estimated+1.96*(atof(methylation_number)/atof(total_number))*sqrt(pow((m3-m1)/1.96,2)/pow(n1,2)+pow(m1,2)*pow((n3-n1)/1.96,2)/pow(n1,4));
	 					if(ratio_estimated_hi>1.0) {
	 						ratio_estimated_hi=1.0;
	 					}
	 					fprintf(fout,"%s\t%s\t%s\t%f\t%f\t%f\n",temp2,temp3,temp4,ratio_estimated,ratio_estimated_lo,ratio_estimated_hi);
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
	 
	fclose(fisoform_filter);
	fclose(fisoform_methylation_filter);
	fclose(fout);
	return 0;
}
