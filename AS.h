#define LINE_BUF_LEN     100
#define INFTY            LONG_MAX
#define TRUE  1
#define FALSE 0

struct point {
  long int citynumber;
  double x;
  double y;
  long int group;
  long int subcitynumber;
};



long int **relationship;


/*struct map {
	double fitness;
	double ave_best;
	double **pheromone;
};*/

struct problem
{
  char          name[LINE_BUF_LEN];      	 /* instance name */
  double        optimum;					/* optimal tour length if known, otherwise a bound */
  long int      n;                      /* number of cities */
  long int      n_near;                 /* number of nearest neighbors */
  struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
  double        **distance;			    /* distance matrix: distance[i][j] gives distance between city i und j */
  long int      **nn_list;              /* nearest neighbor list; contains for each node i a
                                           sorted list of n_near nearest neighbors */
};


typedef struct 
{
  long int  *tour;
  char      *visited;
  double  tour_length;
} ant_struct;

//extern clock_t start_time;
extern double elapsed;

extern long int n_map;

extern double   *best_in_try;		/* best length of each try */
extern long int *best_found_at;		/* iteration when best length is found */
extern double   *time_best_found;	/* time when best length is found  */
extern double   *time_total_run;	/* total time of a try */
extern long int iter_to_best;

extern ant_struct *ant;
extern ant_struct *Mant;			/*memory ants*/
extern ant_struct *best_so_far_ant;
extern double   *prob_of_selection;

extern long int n_ants;			/* number of ants */

extern long int iteration;         /* iteration counter */
extern long int max_iteration;
extern long int max_tries;         /* maximum number of independent tries */
extern double   max_time;          /* maximal allowed run time of a try  */
extern long int optimal;           /* optimal solution or bound to find */

extern double alpha;         /* importance of trail */
extern double beta;          /* importance of heuristic evaluate */
extern double rho;           /* parameter for evaporation */
extern double Q;			  /* pheromone deposition weight */
extern double q_0;           /* probability of best choice in tour construction */

extern long int acs_flag;

extern char name_buf[LINE_BUF_LEN]; /* directory (include file name) of data file */

extern struct problem instance;
extern struct problem cinstance;
extern long int n;          /* number of cities in the instance to be solved */
extern long int tn;

long int x_num;				/* the individ number of x */ 
long int y_num;				/* the individ number of y */

long int *rx_num;
long int *ry_num;

long int redivide_number;

long int cross_flag;

extern	double   trail_0;	/*pheromone initialization */

extern long int**  tmp_tour;
extern long int**  pos;

extern double **pheromone;
extern double **heuristic;
extern double **total;

extern double **phex;


extern long int ras_ranks;       /* additional parameter for rank-based version */
extern double epsilon;			  /* EACO parameter */

extern long int *sort_ant_no;
extern long int *sort_Mant_no;
extern double *save_ant_length;
extern double *save_Mant_length;

extern long int dlb_flag;
extern long int ls_flag;
long int evapor_nn_flag;

extern double  (*distance)(long int, long int);  /* pointer to function returning distance */


void init_program( long int argc, char *argv[] );
void init_try( long int ntry );
void init_pheromone_trails( double initial_trail );
void map_init_pheromone_trails(long int i_map, double initial_trail );
void set_default_parameters();
double ** compute_distances(void);
long int ** compute_nn_lists( void );
void allocate_ants ( void );
double nn_tour( void );
void compute_total_information( void );
void map_compute_total_information( long int i_map );
void compute_heuristic();
void place_ant( ant_struct *a , long int step );
void choose_closest_next( ant_struct *a, long int phase );
void ant_empty_memory( ant_struct *a );
long int map_termination_condition( long int i_map );
long int termination_condition( void );
void construct_solutions( void );
void construct_solutions_revised( void );
void map_construct_solutions_revised( long int i_map );
void memory_construct_solutions( void );
void move_to_next( ant_struct *a , long int phase );
void move_to_next_normalized( ant_struct *a , long int phase );
void colony_ACS_move_to_next( long int n_colony, ant_struct *a , long int phase );
void map_ACS_move_to_next( long int i_map, ant_struct *a , long int phase );
void ACS_move_to_next( ant_struct *a , long int phase );
void map_choose_best_next( long int i_map, ant_struct *a, long int phase );
void choose_best_next( ant_struct *a, long int phase );
void map_neighbour_choose_best_next( long int i_map, ant_struct *a, long int phase );
void neighbour_choose_best_next( ant_struct *a, long int phase );
void map_update_statistics(long int i_map, long int n_try );
void update_statistics( long int );
void copy_from_to(ant_struct *a1, ant_struct *a2);
void pheromone_trail_update( void );
void map_evaporation( long int i_map );
void evaporation( void );
void map_evaporation_nn_list( long int i_map );
void evaporation_nn_list( void );
void as_update( void );
void ras_update( void );
void EACO_ras_update( void );
void global_update_pheromone( ant_struct *a );
void map_global_acs_pheromone_update( long int i_map, ant_struct *a );
void global_acs_pheromone_update( ant_struct *a );
void global_update_pheromone_weighted( ant_struct *a, double weight);
void map_report( long int i_map, long int ntry ) ;
void exit_try( long int ntry );
void exit_program( void );
void copy_ant_to_tour(ant_struct *ant, long int *tour);
void copy_tour_to_ant(long int *tour, ant_struct *ant);
void copy_from_tour_to_tour(long int *src, long int *des);
void add_tmp_phe( ant_struct *a );
void minus_tmp_phe( ant_struct *a, double tmp_length );
void minus_tmp_phe_tour( long int *tour, double tmp_length );
void map_ACS_construct_solutions(long int i_map);
void ACS_construct_solutions();
void EACS_construct_solutions();
void M_ACS_construct_solutions();
void map_local_acs_pheromone_update(long int i_map, ant_struct *a, long int phase );
void local_acs_pheromone_update( ant_struct *a, long int phase );
void EACO_construct_solutions( void );
void M_EACO_construct_solutions( void );
void EACO_pheromone_trail_update( void );
void map_ACS_pheromone_trail_update(long int i_map);
void ACS_pheromone_trail_update(void);
void EACS_pheromone_trail_update(void);
void map_2opt_local_search(long int i_map);
void map_local_search(long int i_map);
void local_search( void );
double individual_local_search(long int* tour );
void two_opt_first( long int *tour );
void three_opt_first( long int *tour );
long int * generate_random_permutation( long int n );
long int check_convergency(long int i_map);
void swap_pheromone(double**a, double**b);
void average_pheromone(double**a, double**b);
void island_crossover(long int island_num);

void island_crossover_c(long int island_num);