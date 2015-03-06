#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"


void params_read( char *file, int *nx, int*nz,double *dx,double *dz, int*ne,int*nt,double *dt,double *f0,int*veloread,char *file_v0,int *nbsou, int *nbrec, char *file_sou){


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


}
