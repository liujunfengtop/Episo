/* compare.c

   Junfeng Liu, 2015
   
   This is an program to see if the locations from two .

   gcc -o compare compare.c -lm
   
   compare <transfile> <pre_transcript_id> <pre_position> <transcript_id> <position>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>

/*The function is used to convert the location in the chromsome to the location in the transcript and to compare*/
int GetNumber (char *transf, char *pre_transcript_id, char *pre_position, char *transcript_id, char *position, char *length);

int main (int argc, char* argv[])
{
   int l;
   
   char trans[1000], pre_transcript_id[1000], pre_position[1000], transcript_id[1000], position[1000], length[1000];

   if(argc>6) {
   	strcpy(trans, argv[1]);
   	strcpy(pre_transcript_id, argv[2]);
   	strcpy(pre_position, argv[3]);
    strcpy(transcript_id, argv[4]);
    strcpy(position, argv[5]);
    strcpy(length, argv[6]);
   } else {
   	return -1;
   }

   l=GetNumber(trans,pre_transcript_id,pre_position,transcript_id,position,length);
     
   return (l);
}

int GetNumber (char *transf, char *pre_transcript_id, char *pre_position, char *transcript_id, char *position, char *length)
{
	/* j for chracter loop about each row of the trans file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
        */
	int k, j, l, read;
	
	/*n1 for the start position of the exon;
	  n2 for the end position of the exon;
	  n3 for the location in the transcript;
	  n4 for the number of the transcript which include the assigned methylation site;
	*/
	long n1, n2, n3, n4, n5, n6, length_remain_1, length_remain_2;
	
	/*ls[] for each row of the trans file;
	  temp[] for each column of the row;
	  temp2[] for chromosome name and the transcript name;
	*/
	char ls[5000], temp[1000], temp2[1000];
	
	char te1, te2, te3, te4;
	
	char pre_chrom_name[1000], chrom_name[1000];
	
	char record_all_1[1000], record_1[1000], record_all_2[1000], record_2[1000];
	
	
	FILE *ftrans=fopen (transf, "r");
	
	n4=0;
	
	n5=0;
	
	n6=0;
	
	length_remain_1=0;
	
	length_remain_2=0;
	
	te3='G';
	
	te4='G';
	
	 /*read each row of the trans file*/
	 while(fgets(ls,5000,ftrans) != NULL) {
	 	k=0;
	 	j=0;
	 	read=0;
	 	n1=0;
	 	n2=0;
	 	n3=0;
	 	te1='0';
	 	te2='0';
	 	memset(temp,0,sizeof(temp));
   	memset(temp2,0,sizeof(temp2));
	 	while(read<1000) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==0) {
	 				for(l=0;l<k;l++) {
	 					temp2[l]=temp[l];
	 				}
	 			}
	 			if(read==1) {
	 				if(strncmp(temp,pre_transcript_id,k)==0) {
	 					te1='1';
	 					strcpy(pre_chrom_name,temp2);
	 				}
	 				if(strncmp(temp,transcript_id,k)==0) {
	 					te2='1';
	 					strcpy(chrom_name,temp2);
	 				}
	 			}
	 			if(te1=='0'&&te2=='0'&&read>1) {
	 				break;
	 			}
	 			if((read>=2)&&(fmod(read,2)==0)) {
	 				n1=atol(temp);
	 			}
	 			if((read>=2)&&(fmod(read,2)==1)) {
	 				n2=atol(temp);
	 			}
	 			/*get the location in the transcript and output*/
	 			if(n2>n1) {
	 				n3=n3+n2-n1+1;
	 				if(te1=='1') {
	 					if(n4>0) {
	 						if(length_remain_1>(n2-n1+1)) {
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n1);
	 							strcat(record_all_1,record_1);
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n2);
	 							strcat(record_all_1,record_1);
	 							length_remain_1=length_remain_1-(n2-n1+1);
	 						} else {
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n1);
	 							strcat(record_all_1,record_1);
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",(n1+length_remain_1-1));
	 							strcat(record_all_1,record_1);
	 							te3='1';
	 							break;
	 						}
	 					}
	 					if((atol(pre_position)<=n3)&&(n4==0)) {
	 						n4=n2+atol(pre_position)-n3;
	 						if((n4+atol(length)-1)<=n2) {
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n4);
	 							strcat(record_all_1,record_1);
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",(n4+atol(length)-1));
	 							strcat(record_all_1,record_1);
	 							te3='1';
	 							break;
	 						} else {
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n4);
	 							strcat(record_all_1,record_1);
	 							memset(record_1,0,sizeof(record_1));
	 							sprintf(record_1,"%ld",n2);
	 							strcat(record_all_1,record_1);
	 							length_remain_1=atol(length)-(n2-n4+1);
	 						}
	 					}
	 				}
	 				if(te2=='1') {
	 					if(n5>0) {
	 						if(length_remain_2>(n2-n1+1)) {
	 						memset(record_2,0,sizeof(record_2));
	 						sprintf(record_2,"%ld",n1);
	 						strcat(record_all_2,record_2);
	 						memset(record_2,0,sizeof(record_2));
	 						sprintf(record_2,"%ld",n2);
	 						strcat(record_all_2,record_2);
	 						length_remain_2=length_remain_2-(n2-n1+1);
	 						} else {
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",n1);
	 							strcat(record_all_2,record_2);
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",(n1+length_remain_2-1));
	 							strcat(record_all_2,record_2);
	 							te4='1';
	 							break;
	 					  }
	 					}
	 					if((atol(position)<=n3)&&(n5==0)) {
	 						n5=n2+atol(position)-n3;
	 						if((n5+atol(length)-1)<=n2) {
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",n5);
	 							strcat(record_all_2,record_2);
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",(n5+atol(length)-1));
	 							strcat(record_all_2,record_2);
	 							te4='1';
	 							break;
	 						} else {
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",n5);
	 							strcat(record_all_2,record_2);
	 							memset(record_2,0,sizeof(record_2));
	 							sprintf(record_2,"%ld",n2);
	 							strcat(record_all_2,record_2);
	 							length_remain_2=atol(length)-(n2-n5+1);
	 						}
	 					}
	 				}
	 			}
	 			if(ls[j]=='\n') {
	 				break;
	 			}
	 			memset(temp,0,sizeof(temp));
	 			read=read+1;
	 			k=0;
	 		} else {
	 			k=k+1;
	 		}
	 		j=j+1;
	 	}
	 	if(te3=='1'&&te4=='1') {
	 		if(strcmp(pre_chrom_name,chrom_name)==0&&strcmp(record_all_1,record_all_2)==0) {
	 			n6=1;
	 		} else {
	 			n6=0;
	 		}
	 		break;
	 	}
	}
	fclose(ftrans);
	return(n6);
}

