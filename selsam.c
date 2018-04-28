/* selsam.c

   Junfeng Liu, 2015
   
   This is an program to select the rows which include 'methylated_liu' from the file accepted_hits.sam.

   gcc -o selsam selsam.c -lm
   
   selsam <accpted_hits.sam> <skipped_number>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


/*The function is used to select the rows which include 'methylated_liu' from the file accepted_hits.sam and to output accepted_hits_methylation file*/
int SelSam (char *hitsf, char *skipped_number);

int main (int argc, char* argv[])
{
   int l;
   
   char hits[1000], skipped_number[1000];

   if(argc>2) {
   	strcpy(hits, argv[1]);
   	strcpy(skipped_number, argv[2]);
   } else {
   	return -1;
   }

   l=SelSam(hits,skipped_number);
     
   return (l);
}

int SelSam (char *hitsf, char *skipped_number)
{
	/* j for chracter loop about each row of the hits file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;
	
	
	/*ls[] for each row of the hits file;
	  temp[] for each column of the row;
	*/
	char ls[5000], temp[1000];
	
	
	FILE *fhits, *fout;
	
	if((fhits=fopen(hitsf, "r"))==NULL) {
   		return -1;
   }
   
  
	if((fout=fopen("accepted_hits_methylation.sam", "w"))==NULL) {
   		return -1;
   }
  
	/*skipping the first skipped_number rows of the hits file*/
   for(l=0;l<atol(skipped_number);l++) {
   	fgets(ls,5000,fhits);
   	fprintf(fout,"%s",ls);
   }
	
	 /*read each row of the hits file*/
	 while(fgets(ls,5000,fhits) != NULL) {
	 	k=0;
	 	j=0;
	 	read=0;
	 	memset(temp,0,sizeof(temp));
	 	while(read<1) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==0) {
	 				if(strncmp(temp,"methylated_liu",14)==0) {
	 					fprintf(fout,"%s",ls);
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
	fclose(fhits);
	fclose(fout);
	return 0;
}

