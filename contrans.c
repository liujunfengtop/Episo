/* contrans.c

   Junfeng Liu, 2015
   
   This is an program to convert the gtf file given by cufflinks program to the fa file for bismark program.

   gcc -o contrans contrans.c
   
   contrans <ctlfile>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<ctype.h> /*add head file for compile warning about isalnum function. by liujf 2017-05-18*/

struct CommonInfo {
   char *z[3], *spname[3], outf[128], seqf[128], gtff[128], faf[128], transf[128], ratef[128], ctlf[128], seqlen[128], fix_locusrate;
   int model, ncode, cleandata, seed, npoints, ncatBeta, UseMedianBeta, getSE;
   int ndata, ngene, seqtype, ns, ls, posG[1+1], lgene[1], *pose, npatt, readpattern;
   int *Nij, nGtree;
   double *fpatt, kappa, alpha, rho, rgene[1], pi[4], piG[1][4];
   double *lnLmax, *locusrate;
   double *pDclass, *tau1beta, *bp0124[5], *wwprior[5];
}  com;




int GetOptions (char *ctlf);


FILE *fout, *fgtf, *ffa, *ftrans, *fseq;

int main (int argc, char* argv[])
{
   char VerStr[32] = "Version 1.0, March 2015";
   
   /*i for row loop about the gtf file; j for chracter loop about each row of the gtf file; read for column loop about each row;
     k for character loop about each column; l for 'for loop'; m for counter to print '\n'
   */
   int i,j,read,k,l,m;
   /*read1 for the length of the chromosome name*/
   int read1;
   /*read2 for the position that character is in fafile*/
   long read2=1;
   /*n1 for the start position of the transcript or exon on each row of the gtf file;
     n2 for the end position of the transcript or exon on each row of the gtf file;
     n3 for the current position of ffa;
   */
   long n1,n2,n3;
   /*ls[] for each row of the gtf file; ls1[] for each row of ffa file*/
   char ls[1000],ls1[1000];
   /*temp[] for each column of each row of the gtf file;
   temp1[] for chromosome name of fa file;
   temp2[] for transcript or exon;temp3[] for transcrip_id;
   temp4[] for chromosome name of gtf file,temp4_pre[] for the previous chromosome name of gtf file
   */
   char temp[1000],temp1[1000],temp2[20],temp3[100],temp4[1000],temp4_pre[1000];
   /*te for the character*/
   char te;

   /*read the control file*/
   strcpy(com.ctlf, "contrans.ctl");
   if(argc>1) strcpy(com.ctlf, argv[1]);
   GetOptions (com.ctlf);

   /*open the file*/
   if((fout=fopen(com.outf, "w"))==NULL) {
   		printf("Open %s file Error: %s\n",com.outf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.outf,strerror(errno));
   		return -1;
   }
   
   if((ftrans=fopen(com.transf, "w"))==NULL) {
   	  printf("Open %s file Error: %s\n",com.transf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.transf,strerror(errno));
   		return -1;
   }
   
   if((fseq=fopen(com.seqf, "w"))==NULL) {
   	  printf("Open %s file Error: %s\n",com.seqf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.seqf,strerror(errno));
   		return -1;
   }
   
   if((fgtf=fopen(com.gtff,"r"))==NULL) {
   	  printf("Open %s file Error: %s\n",com.gtff,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.gtff,strerror(errno));
   		return -1;
   }
   
   if((ffa=fopen(com.faf,"rb"))==NULL) {
   	  printf("Open %s file Error: %s\n",com.faf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.faf,strerror(errno));
   		return -1;
   }
   
   printf("contrans (%s)\n", VerStr);
   fprintf(fout, "contrans (%s)\n", VerStr);
   
   i=0;
   m=0;
   
   /*read each row of gtf file to generate trans file which contains the information about the transcript and seq file which contains DNA sequence about the transcript*/
   while(fgets(ls,1000,fgtf) != NULL) {
   	  k=0;
   	  j=0;
   	  read=0;
   	  /*generate the information about the transcript for the trans file and the DNA sequence for the seq file according ls[j]*/
   	  while(read<12) {
   	  	temp[k]=ls[j];
   	  	/*read each column of each row*/
   	  	if(ls[j]==' '||ls[j]==','||ls[j]=='\t') {
   	  		/*output the chromosome name for trans file and find the same chromosome from the fa file when got the first column*/
   	  		if(read==0) {
   	  			for(l=0;l<k;l++) {
   	  				temp4[l]=temp[l];
   	  			}
   	  			read1=k;
   	  		}
   	  		if(read==2) {
   	  			for(l=0;l<k;l++) {
   	  				temp2[l]=temp[l];
   	  			}
   	  		}
   	  		if(read==3) {
   	  			n1=atol(temp);
   	  		}
   	  		if(read==4) {
   	  			n2=atol(temp);
   	  		}
   	  		/*generate the information about the transcript for the trans file and the DNA sequence for the seq file when got the twelfth column*/
   	  		if(read==11) {
   	  			/*output the information when the row of gtf file is transcript*/
   	  			if(strncmp(temp2,"transcript",10)==0) {
   	  				/*remove the pointer to the beginning position of the ffa file when temp4_pre!=temp4 
   	  				and remove the pointer to the position of the given chromosome in the ffa file when temp4_pre=temp4
   	  				and reset the variable read2
   	  				*/
   	  				if(strlen(temp4_pre)==0||strcmp(temp4_pre,temp4)!=0) {
   	  					rewind(ffa);
   	  				  /*locate the position of the given chromosome in the fa file*/
   	  				  while(fgets(ls1,1000,ffa)!=NULL) {
   	  					  for(l=0;l<read1;l++) {
   	  						  temp1[l]=ls1[l+1];
   	  					  }
   	  					  //if(strncmp(temp1,temp4,read1)==0&&read1==(strlen(ls1)-2)) { //by ljf 2017-5-18
   	  					  if(strncmp(temp1,temp4,read1)==0&&(ls1[read1+1]==' '||ls1[read1+1]=='\t'||ls1[read1+1]=='\n')) { //add the case mouse.fa by ljf 2017-05-18
   	  						  break;
   	  					  }
   	  				  }
   	  				  n3=ftell(ffa);
   	  				  if(n3==1L) {
   	  				  	printf("re-location error\n");
   	  				  	fprintf(fout,"chrosome %s:re-location error\n",temp4);
   	  				  }
   	  					strcpy(temp4_pre,temp4);
   	  				} else {
   	  					rewind(ffa);
   	  					fseek(ffa,n3,0);
   	  				}
   	  				read2=1;
   	  				/*output the information for trans file and seq file*/
   	  				if(i==0) {
   	  					for(l=0;l<read1;l++) {
   	  						fprintf(ftrans,"%c",temp4[l]);
   	  					}
   	  					fprintf(ftrans,"\t");
   	  					fprintf(fseq,">");
   	  					for(l=0;l<k-3;l++) {
   	  						fprintf(ftrans,"%c",temp[l+1]);
   	  						fprintf(fseq,"%c",temp[l+1]);
   	  					}
   	  					fprintf(fseq,"\n");
   	  				} else {
   	  					fprintf(ftrans,"\n");
   	  					for(l=0;l<read1;l++) {
   	  						fprintf(ftrans,"%c",temp4[l]);
   	  					}
   	  					fprintf(ftrans,"\t");
   	  					if(m>0&&m<atol(com.seqlen)) {
   	  						fprintf(fseq,"\n>");
   	  						m=0;
   	  					} else {
   	  						fprintf(fseq,">");
   	  					}
   	  					for(l=0;l<k-3;l++) {
   	  						fprintf(ftrans,"%c",temp[l+1]);
   	  						fprintf(fseq,"%c",temp[l+1]);
   	  					}
   	  					fprintf(fseq,"\n");
   	  				}
   	  			}
   	  			/*output the information when the row of gtf file is exon*/
   	  			if(strncmp(temp2,"exon",4)==0) {
   	  				/*output the information for trans file*/
   	  				fprintf(ftrans,"\t%ld\t%ld",n1,n2);
   	  				/*output the information for seq file*/
   	  				te=fgetc(ffa);
   	  				if(read2!=1&&te!='\n') {
   	  					read2=read2+1;
   	  				}
   	  				while(te!=EOF) {
   	  					if(read2==n2) {
   	  						fprintf(fseq,"%c",te);
   	  						m=m+1;
   	  						if(m==atol(com.seqlen)) {
   	  							fprintf(fseq,"\n");
   	  							m=0;
   	  						}
   	  						break;
   	  					} else {
   	  						if(read2>=n1&&read2<n2) {
   	  							if(te!='\n') {
   	  								fprintf(fseq,"%c",te);
   	  								m=m+1;
   	  								if(m==atol(com.seqlen)) {
   	  									fprintf(fseq,"\n");
   	  									m=0;
   	  								}
   	  							}
   	  						}
   	  						te=fgetc(ffa);
   	  						if(te!='\n') {
   	  							read2=read2+1;
   	  						}
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
   	  i=i+1;
    }
 

   fprintf(ftrans,"\n");   //add \n at the last line by ljf 2017-05-19
   fprintf(fseq,"\n");     //add \n at the last line by ljf 2017-05-19

   printf("Completing! Please check the file (%s) and (%s)\n", com.transf,com.seqf);
   fprintf(fout,"Completing! Please check the file (%s) and (%s)\n", com.transf,com.seqf);
   


   fclose(ftrans);
   fclose(fseq);
   fclose(fout);
   fclose(fgtf);
   fclose(ffa);
   return 0;
}


int GetOptions (char *ctlf)
{
   int iopt,i, nopt=9, lline=4096;
   char line[4096],*pline, opt[32], *comment="*#", *seqerrstr="0EF";
   char *optstr[] = {"outfile","gtffile","fafile","transfile","seqfile","seqlength"};
   double t=1;
   FILE  *fctl=fopen (ctlf, "r");
   

   if (fctl) {
      for (;;) {
         if(fgets(line, lline, fctl) == NULL) break;
         if(line[0]=='/' && line[1]=='/') 
            break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (strchr(comment,line[i])) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "="))==NULL)
            continue;

         for (iopt=0; iopt<nopt; iopt++) {
            if (strncmp(opt, optstr[iopt], 9)==0)  {
               switch (iopt) {
                  case ( 0): sscanf(pline+1, "%s", com.outf);    break;
                  case ( 1): sscanf(pline+1, "%s", com.gtff);    break;
                  case ( 2): sscanf(pline+1, "%s", com.faf);     break;
                  case ( 3): sscanf(pline+1, "%s", com.transf);  break;
                  case ( 4): sscanf(pline+1, "%s", com.seqf);    break;
                  case ( 5): sscanf(pline+1, "%s", com.seqlen);  break;
               }
               break;
            }
         }
         if (iopt==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose(fctl);
   }

   return (0);
}


