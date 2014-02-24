#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876


extern clock_t start_time;
extern long int seed;


void sort2(double v[], long int v2[], long int left, long int right);
void swap2(double v[], long int v2[], long int i, long int j);
double ** generate_double_matrix( long int n, long int m);
long int ** generate_int_matrix( long int n, long int m);
void map_start_timers(long int i_map);
void start_timers();
double map_elapsed_time(long int i_map);
double elapsed_time();


double ran01( long *idum );
double compute_tour_length( long int *t );
long int map_find_best( long int i_map );
long int find_best( void );
void output_solution( void );
void checkTour( long int *t);
void printTour( long int *t );

double mean( long int *values, long int max );
double meanr( double *values, long int max );
double std_deviation( long int *values, long int max, double mean );
double std_deviationr( double *values, long int max, double mean );
double best_of_vector( double *values, long int l );
double worst_of_vector( double *values, long int l );
void empty_memory( long int*  a, long int size );

void swap(long int v[], long int i, long int j);
void empty_matrix_memory( long int**  a, long int n, long int m);


long int roulette(double* elments, long element_size);














