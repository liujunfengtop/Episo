/* selsam-single.c

   Junfeng Liu, 2015
   
   This is an program to select the rows which include assigned site from the file accepted_hits.sam and methylation_summary_sam.

   gcc -o selsam-single selsam-single.c
   
   selsam-single <accpted_hits.sam> <methylation_summary_sam> <skipped_number> <tag>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


/*The function is used to select the rows which include assigned site from the file accepted_hits.sam and methylation_summary_sam*/
int SelSam (char *hitsf, char *summf, char *skipped_number, char *tag, char *read_number, char *tag1, char *out_name);

int main (int argc, char* argv[])
{
   int l;
   
   char hits[1000], summ[1000], skipped_number[1000], tag[1000], read_number[1000], tag1[1000], out_name[1000];

   if(argc>7) {
   	strcpy(hits, argv[1]);
   	strcpy(summ, argv[2]);
   	strcpy(skipped_number, argv[3]);
   	strcpy(tag, argv[4]);
   	strcpy(read_number, argv[5]);
   	strcpy(tag1, argv[6]);
   	strcpy(out_name, argv[7]);
   } else {
   	return -1;
   }

   l=SelSam(hits,summ,skipped_number,tag, read_number, tag1, out_name);
     
   return (l);
}

int SelSam (char *hitsf, char *summf, char *skipped_number, char *tag, char *read_number, char *tag1, char *out_name)
{
	/* j for chracter loop about each row of the hits file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read, n, m;
	
	
	/*ls[] for each row of the hits file;
	  ls1[] for each row of the summ file;
	  temp[] for each column of the row from the hits file;
	  temp1[] for each column of the row from the summ file;
	*/
	char ls[5000], ls1[5000], temp[1000], temp1[1000];
	
	char **summ;
	
	
	
	FILE *fhits, *fout, *fsumm;
	
	if((fhits=fopen(hitsf, "r"))==NULL) {
   		return -1;
   }
   
  if((fsumm=fopen(summf, "r"))==NULL) {
   		return -1;
   }
   
  
	if((fout=fopen(out_name, "w"))==NULL) {
   		return -1;
   }
  
  
	/*skipping the first skipped_number rows of the hits file*/
   for(l=0;l<atol(skipped_number);l++) {
   	fgets(ls,5000,fhits);
   	if(tag1[0]=='1') {
   		fprintf(fout,"%s",ls);
   	}
   }
	
	 n=0;
	 
	 while(fgets(ls1,5000,fsumm) != NULL) {
	 	n=n+1;
	 }
	 
	 rewind(fsumm);
	 
	 
	 summ=(char**)malloc(n*sizeof(char*));
	 
	 
	 for(l=0;l<n;l++) {
	 	summ[l]=(char*)malloc(1000*sizeof(char));
	}
	
	n=0;
	
	 while(fgets(ls1,5000,fsumm) != NULL) {
	 	for(l=0;l<strlen(ls1)-2;l++) {
	 		summ[n][l]=ls1[l];
	 	}
	 	n=n+1;
	}
	
   m=0;
	
	 /*read each row of the hits file*/
	 while(fgets(ls,5000,fhits) != NULL) {
	 	k=0;
	 	j=0;
	 	read=0;;
	 	memset(temp,0,sizeof(temp));
	 	while(read<1) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==0) {
	 				//rewind(fsumm);
	 				/*
	 				while(fgets(ls1,5000,fsumm) != NULL) {
	 					memset(temp1,0,sizeof(temp1));
	 					for(l=0;l<strlen(ls1)-2;l++) {
	 						temp1[l]=ls1[l];
	 					}
	 					if(strcmp(temp,temp1)==0||strstr(temp,temp1)!=NULL) {
	 						if(tag[0]=='1') {
	 							if(strncmp(temp,"methylated_liu",14)==0) {
	 								fprintf(fout,"%s",ls);
	 							}
	 						} else {
	 							fprintf(fout,"%s",ls);
	 						}
	 						break;
	 					}
	 				}
	 				*/
	 				for(l=0;l<n;l++) {
	 					//te='1';
	 					//for(m=0;m<strlen(temp);m++) {
	 						
	 					//}
	 					if(strcmp(temp,summ[l])==0||strstr(temp,summ[l])!=NULL) {
	 						if(tag[0]=='1') {
	 							if(strncmp(summ[l],"methylated_liu",14)==0) {
	 								fprintf(fout,"%s",ls);
	 							}
	 						} else {
	 							fprintf(fout,"%s",ls);
	 						}
	 						break;
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
	 	m=m+1;
	 	if(m==atol(read_number)) {
	 		break;
	 	}
	}
	
	free(summ);
	
	fclose(fhits);
	fclose(fout);
	fclose(fsumm);
	return 0;
}

