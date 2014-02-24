#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"utilities.h"
#include"AS.h"
#include"IO.h"
#include"Map.h"

clock_t start_time;
long int seed = 12345678;

void sort2(double v[], long int v2[], long int left, long int right)
/*    
      FUNCTION:       recursive routine (quicksort) for sorting one array; second 
                      arrays does the same sequence of swaps  
      INPUT:          two arrays, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  long int k, last;

  if (left >= right) 
    return;
  swap2(v, v2, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap2(v, v2, ++last, k);
  swap2(v, v2, left, last);
  sort2(v, v2, left, last);
  sort2(v, v2, last+1, right);
}


void swap2(double v[], long int v2[], long int i, long int j)
/*    
      FUNCTION:       auxiliary routine for sorting an integer array  
      INPUT:          two arraya, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  double tmp;
  long int tmp2;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
  tmp2 = v2[i];
  v2[i] = v2[j];
  v2[j] = tmp2;
}


double ** generate_double_matrix( long int n, long int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{

  long int i;
  double **matrix;

  if((matrix = malloc(sizeof(double) * n * m +
		      sizeof(double *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*m;
  }
  return matrix;
}

long int ** generate_int_matrix( long int n, long int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{

  long int i;
  long int **matrix;

  if((matrix = malloc(sizeof(long int) * n * m +
		      sizeof(long int *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (long int *)(matrix + n) + i*m;
  }
  return matrix;
}
void map_start_timers(long int i_map)
/*    
      FUNCTION:		  start timer
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:   
*/
{
    map[i_map].start_time = clock();
}
void start_timers()
/*    
      FUNCTION:		  start timer
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:   
*/
{
    start_time = clock();
}
double map_elapsed_time(long int i_map)
	
/*    
      FUNCTION:       return the time used in seconds (virtual or real, depending on type) 
      INPUT:          TIMER_TYPE (virtual or real time)
      OUTPUT:         seconds since last call to start_timers (virtual or real)
      (SIDE)EFFECTS:  none
*/
{
    double interval;

	interval = clock()- map[i_map].start_time;
    return interval / CLOCKS_PER_SEC;
}
double elapsed_time()
	
/*    
      FUNCTION:       return the time used in seconds (virtual or real, depending on type) 
      INPUT:          TIMER_TYPE (virtual or real time)
      OUTPUT:         seconds since last call to start_timers (virtual or real)
      (SIDE)EFFECTS:  none
*/
{
    double interval;

	interval = clock()- start_time;
    return interval / CLOCKS_PER_SEC;
}




double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

double compute_tour_length( long int *t ) 
/*    
      FUNCTION: compute the tour length of tour t
      INPUT:    pointer to tour t
      OUTPUT:   tour length of tour t
*/
{
    long int      i;
    double tour_length = 0.0;
  
    for ( i = 0 ; i < n ; i++ )
	{
		tour_length += instance.distance[t[i]][t[i+1]];
    }
    return tour_length;
}

long int map_find_best( long int i_map ) 
/*    
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   min=0.0;
    long int   k, k_min;

    min = map[i_map].ant[0].tour_length;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
	if( map[i_map].ant[k].tour_length < min ) {
	    min = map[i_map].ant[k].tour_length;
	    k_min = k;
	}
    }
    return k_min;
}
long int find_best( void ) 
/*    
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   min;
    long int   k, k_min;

    min = ant[0].tour_length;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
	if( ant[k].tour_length < min ) {
	    min = ant[k].tour_length;
	    k_min = k;
	}
    }
    return k_min;
}


void output_solution( void ) 
/*    
      FUNCTION:       output a solution together with node coordinates
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       not used in the default implementation but may be useful anyway
*/
{

    long int i;
    
    for ( i = 0 ; i < n ; i++ ) {
	fprintf(report," %ld %f %f\n",(*best_so_far_ant).tour[i],instance.nodeptr[(*best_so_far_ant).tour[i]].x, instance.nodeptr[(*best_so_far_ant).tour[i]].y);
    }
    printf("\n"); 
}

void checkTour( long int *t ) 
/*    
      FUNCTION:       make a simple check whether tour *t can be feasible
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i, sum=0;

    for( i = 0 ; i < n ; i++ ) {
	sum += t[i];
    }
    if ( sum != (n-1) * n / 2 ) {
	fprintf(stderr,"Next tour must be flawed !!\n");
	printTour( t );
	exit(1);
    }
}
void printTour( long int *t ) 
/*    
      FUNCTION:       print the tour *t
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i;

    printf("\n");
    for( i = 0 ; i <= n ; i++ ) {
	if (!i%25)
	    printf("\n");
	printf("%ld ", t[i]);
    }
    printf("\n");
    printf("Tour length = %ld\n\n",compute_tour_length( t ));
}


double mean( long int *values, long int max ) 
/*    
      FUNCTION:       compute the average value of an integer array of length max 
      INPUT:          pointer to array, length of array
      OUTPUT:         average 
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += (double)values[j];
  }
  m = m / (double)max;
  return m;
}



double meanr( double *values, long int max ) 
/*    
      FUNCTION:       compute the average value of a floating number array of length max 
      INPUT:          pointer to array, length of array
      OUTPUT:         average 
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += values[j];
  }
  m = m / (double)max;
  return m;
}



double std_deviation( long int *values, long int max, double mean ) 
/*    
      FUNCTION:       compute the standard deviation of an integer array  
      INPUT:          pointer to array, length of array, mean 
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   dev = 0.;

  if (max <= 1)
    return 0.;
  for ( j = 0 ; j < max; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}



double std_deviationr( double *values, long int max, double mean ) 
/*    
      FUNCTION:       compute the standard deviation of a floating number array  
      INPUT:          pointer to array, length of array, mean 
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   dev;

  if (max <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}


double best_of_vector( double *values, long int l ) 
/*    
      FUNCTION:       return the minimum value in an integer value  
      INPUT:          pointer to array, length of array
      OUTPUT:         smallest number in the array
      (SIDE)EFFECTS:  none
*/
{
  double min;
  long int k;

  k = 0;
  min = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] < min ) {
      min = values[k];
    }
  }
  return min;
}



double worst_of_vector( double *values, long int l ) 
/*    
      FUNCTION:       return the maximum value in an integer value  
      INPUT:          pointer to array, length of array
      OUTPUT:         largest number in the array
      (SIDE)EFFECTS:  none
*/
{
  double max;
  long int k;

  k = 0;
  max = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] > max ){
      max = values[k];
    }
  }
  return max;
}

void empty_memory( long int*  a, long int size ) 
/*    
      FUNCTION:       empty the ants's memory regarding visited cities
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited cities is reinitialized to FALSE
*/
{
    long int   i;

    for( i = 0 ; i < size ; i++ ) {
	a[i]=FALSE;
    }
}

void empty_matrix_memory( long int**  a, long int n, long int m) 
/*    
      FUNCTION:       empty the ants's memory regarding visited cities
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited cities is reinitialized to FALSE
*/
{
    long int   i,j;

    for( i = 0 ; i < n ; i++ )
	{
		for( j = 0 ; j < m ; j++ )
			a[i][j]=FALSE;
    }
}

void swap(long int v[], long int i, long int j)
/*    
      FUNCTION:       auxiliary routine for sorting an integer array  
      INPUT:          array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of array are swapped
*/
{
  long int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}

long int roulette(double* elments, long element_size)
{
	int i=0;
	double sum=0.0;
	double rnd=0.0;
	double* prob=malloc( element_size*sizeof(double) );
	rnd = ran01( &seed );
	for(i=0;i<element_size; i++)
	{
		sum += elments[i];
	}
	for(i=0;i<element_size; i++)
	{
		prob[i]=elments[i]/sum;
	}
	i=0;
	sum=prob[i];
	while(sum<=rnd)
	{
		i++;
		sum+=prob[i];
	}
	free(prob);
	return i;
}
