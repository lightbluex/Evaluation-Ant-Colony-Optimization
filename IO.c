#include<stdio.h>
#include<stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include"AS.h"
#include"IO.h"

long int n;          /* number of cities in the instance to be solved */
long int tn;
double  (*distance)(long int, long int);  /* pointer to function returning distance */
FILE *report;
FILE *best_report;
FILE *try_result;

double round_distance (long int i, long int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
	/*if integer,insert (long int), if double, delete (long int) */
 
	double r  = sqrt(xd*xd + yd*yd) + 0.5;
    return (long int) r;
	
	/*if integer,insert (long int), if double, delete (long int) */
/*
	double r  = sqrt(xd*xd + yd*yd);
    return r;
*/

}


struct point * read_etsp(const char *tsp_file_name) 
/*    
      FUNCTION: read instance file
      INPUT:    instance name
      OUTPUT:   list of coordinates for all nodes
      COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
*/
{
    FILE         *tsp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    struct point *nodeptr;

    tsp_file = fopen(tsp_file_name, "r");
    if ( tsp_file == NULL ) {
	fprintf(stderr,"No instance file specified, abort\n");
	exit(1);
    }
    assert(tsp_file != NULL);
    printf("\nreading tsp-file %s ... \n\n", tsp_file_name);

    fscanf(tsp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
	if ( strcmp("NAME", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); )
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance.name, buf);
	    TRACE ( printf("%s \n", instance.name); )
	    buf[0]=0;
	}
	else if ( strcmp("NAME:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance.name, buf);
	    TRACE ( printf("%s \n", instance.name); )
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); )
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT:", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); )
	    buf[0]=0;
	}
	else if ( strcmp("TYPE", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE( printf("%s ", buf); );
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if ( strcmp("TYPE:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if( strcmp("DIMENSION", buf) == 0 ){
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); );
	    fscanf(tsp_file, "%ld", &n);
	    instance.n = n;
	    TRACE ( printf("%ld\n", n); );
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if ( strcmp("DIMENSION:", buf) == 0 ) {
	    fscanf(tsp_file, "%ld", &n);
	    instance.n = n;
	    TRACE ( printf("%ld\n", n); );
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); );
	    buf[0]=0;
	}
	else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); );
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); );
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); );
	   if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
 //	    if ( strcmp("CEIL_2D", buf) == 0 ) {
//		distance = ceil_distance;
//	    }
	    //else if ( strcmp("GEO", buf) == 0 ) {
		//distance = geo_distance;
	    //}
	    //else if ( strcmp("ATT", buf) == 0 ) {
		//distance = att_distance;
	    //}
	    else
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
	    //strcpy(instance.edge_weight_type, buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
	       EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
		printf("%s\n", buf);
	    printf("%s\n", buf);
	   if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	   }
 //	    if ( strcmp("CEIL_2D", buf) == 0 ) {
//		distance = ceil_distance;
// 	    }
	    //else if ( strcmp("GEO", buf) == 0 ) {
		//distance = geo_distance;
	    //}
	    //else if ( strcmp("ATT", buf) == 0 ) {
		//distance = att_distance;
	    //}
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
		exit(1);
	    }
	    //strcpy(instance.edge_weight_type, buf);
	    buf[0]=0;
	}
	buf[0]=0;
	fscanf(tsp_file,"%s", buf);
    }


    if( strcmp("NODE_COORD_SECTION", buf) == 0 ){
	TRACE ( printf("found section contaning the node coordinates\n"); )
	    }
    else{
	fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
	exit(1);
    }

    if( (nodeptr = malloc(sizeof(struct point) * n)) == NULL )
	exit(EXIT_FAILURE);
    else {
	for ( i = 0 ; i < n ; i++ ) {
	    fscanf(tsp_file,"%ld %lf %lf", &nodeptr[i].citynumber, &nodeptr[i].x, &nodeptr[i].y);
//		printf("\n%ld %lf %lf", nodeptr[i].citynumber, nodeptr[i].x, nodeptr[i].y);
	}
    }
    TRACE ( printf("number of cities is %ld\n",n); )
    TRACE ( printf("\n... done\n"); )
	return (nodeptr);
}
void write_params( void ) 
/*    
      FUNCTION:       writes chosen parameter settings in standard output and in 
                      report files 
      INPUT:          none
      OUTPUT:         none
*/
{
  printf("\nParameter-settings of:  %s\n\n", instance.name);

 
  printf("neighbor city number %ld\n", instance.n_near);

  printf("max-tries %ld\n", max_tries);
  printf("optimum %ld\n", optimal);
  printf("time %f\n", max_time);
  printf("num-ants %ld\n", n_ants);

  printf("alpha %f\n", alpha);
  printf("beta %f\n", beta);
  printf("rho %f\n", rho);
  printf("Q %f\n", Q);

  printf("x_num %ld\n",x_num);
  printf("y_num %ld\n",y_num);

  printf("divided map %ld\n",n_map-1);



  printf("\n");
  fprintf(report,"\nParameter-settings of: %s\n\n", instance.name);
  fprintf(report,"neighbor city number %ld\n", instance.n_near);
  fprintf(report,"max-tries %ld\n", max_tries);
  fprintf(report,"optimum %ld\n", optimal);
  fprintf(report,"time %f\n", max_time);
  fprintf(report,"num-ants %ld\n", n_ants);
  fprintf(report,"alpha %f\n", alpha);
  fprintf(report,"beta %f\n", beta);
  fprintf(report,"rho %f\n", rho);
  fprintf(report,"Q %f\n", Q); 

  fprintf(report,"x_num %ld\n",x_num);
  fprintf(report,"y_num %ld\n",y_num);

  fprintf(report,"divided map %ld\n",n_map-1);

  fprintf(report,"\n");
  fflush(report);

    fprintf(report,"\nParameter-settings of: %s\n\n", instance.name);
  fprintf(best_report,"neighbor city number %ld\n", instance.n_near);
  fprintf(best_report,"max-tries %ld\n", max_tries);
  fprintf(best_report,"optimum %ld\n", optimal);
  fprintf(best_report,"time %f\n", max_time);
  fprintf(best_report,"num-ants %ld\n", n_ants);
  fprintf(best_report,"alpha %f\n", alpha);
  fprintf(best_report,"beta %f\n", beta);
  fprintf(best_report,"rho %f\n", rho);
  fprintf(best_report,"Q %f\n", Q); 
  fprintf(best_report,"\n");
  fflush(best_report);

}
