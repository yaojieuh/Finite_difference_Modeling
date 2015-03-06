#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"


void params_read( char *file,FILE *file1, int *nx, int*nz,double *dx,double *dz, int*ne,int*nt,double *dt,double *f0,int*veloread,char *file_v0,int *nbsou, int *nbrec, char *file_sou){


  FILE *file_id;
  char *line, *string, *head;
  int  done = 0, i;

  // open file
  if ( (file_id = fopen(file, "r")) == NULL )
  {                                         
    fprintf(stderr, "\n Problem opening file %s", file);
    exit(EXIT_FAILURE);                                
  }
  int MAX_LINE_LEN=140;
  string = (char *) malloc(MAX_LINE_LEN);
  line = (char *) malloc(MAX_LINE_LEN);
  head = (char *) malloc(MAX_LINE_LEN);


  fgets(line, MAX_LINE_LEN, file_id);

  
  strcpy(string, "read_param");
  while ( fgets(line, MAX_LINE_LEN, file_id) && !done )
  {
    if ( strstr(line, string) != NULL )
    {
      while ( fgets(line, MAX_LINE_LEN, file_id) && !done )
      {
        if ( strstr(line, "dx") != NULL )
        {
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *dx = (double)atof(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *dz = (double)atof(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nx = atoi(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nz = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *ne = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nt = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *dt = (double)atof(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *f0 = (double)atof(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *veloread = atoi(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(file_v0, strtok(NULL, "'")); 
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nbsou = atoi(strtok(head, " "));      
	  fgets(line, MAX_LINE_LEN, file_id);
	  strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nbrec = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(file_sou, strtok(NULL, "'")); 
	  done = 1;
	}
      }
      done = 0;
    }
  }
  

  free(line);
  free(head);
  free(string);
  

  fprintf(file1, "--!\tRead parameter file for 2D FD Modeling...\n");
  fprintf(file1, "--!\t\tFull Model dimension : nz %d nx %d \n", *nz, *nx);
  fprintf(file1, "--!\t\tFull Model discretization [m]: dz %lf dx %lf \n", *dz, *dx);
  fprintf(file1, "--!\t\tTime step number %d\n",*nt);
  fprintf(file1, "--!\t\tTime step %lf\n", *dt);
  fprintf(file1, "--!\t\tSource domininat frenquency %lf\n", *f0);
  fprintf(file1, "--!\t\tTaper dimension: ne %d\n", *ne);
  if(*veloread==1) fprintf(file1, "--!\t\tRead P-wave velocity from %s\n", file_v0);
  if(*veloread==0) fprintf(file1, "--!\t\tGenerate P-wave velocity file in the code\n");
  fprintf(file1, "--!\t\tNumber of source %d\n", *nbsou);
  fprintf(file1, "--!\t\tNumber of reciver %d\n", *nbrec);
  fprintf(file1, "--!\t\tRead  source coordinate from %s\n", file_sou);

}
