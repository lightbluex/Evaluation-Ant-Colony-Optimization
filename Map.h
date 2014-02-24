typedef struct 
{
	double   **pheromone;
	double   **total;

	double   *best_in_try;		/* best length of each try */
	long int *best_found_at;	/* iteration when best length is found */
	double   *time_best_found;	/* time when best length is found  */
	double   *time_total_run;	/* total time of a try */

	long int iteration;			/*count iteration*/
	long int iter_to_best;		/*record best so far iteration*/
	long int convergance_num;   /*count convergance number*/
	long int mn;
	clock_t  start_time;
	double	 elapsed;

	ant_struct *ant;
	ant_struct *Mant;			/*memory ants*/
	ant_struct *best_so_far_ant;/*ont ant, best so far */
	double     *prob_of_selection;/*the probability point for the roullete*/



} map_struct;

extern map_struct * map;

void generate_map( long int n_map );
void init_map(long int i_map, long int n_map );
void map_allocate_ants(long int i_map);
void end_map(long int s_map);