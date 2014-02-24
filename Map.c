#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#include"AS.h"
#include"Map.h"
#include"utilities.h"


map_struct *map;

void generate_map( long int n_map )
/*
    COMMENT: Dont forget to delete memory
*/
{
	long int i;

	if((map = malloc(sizeof( map_struct ) * n_map)) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
    }


/*
	for ( i = 0 ; i < n_map ; i++ ) 
	{
		map[i].pheromone = generate_double_matrix( n, n );
		map[i].total = generate_double_matrix( n, n );

		map[i].best_in_try = calloc(max_tries, sizeof(double));
		map[i].best_found_at = calloc(max_tries, sizeof(long int));
		map[i].time_best_found = calloc(max_tries, sizeof(double));
		map[i].time_total_run = calloc(max_tries, sizeof(double));

		map_allocate_ants(i);
*/
		/*
		if(acs_flag)
			trail_0 = 1. / ( (double) n * (double) nn_tour( ) ) ;
		else
			trail_0 = 1. / ( (rho) * nn_tour() );
		map_init_pheromone_trails(i, trail_0 );//initiate phe of each map
		map_compute_total_information(i);

		map[i].best_so_far_ant->tour_length = INFTY;
		*/
/*
   }
////////////////////change by 2014.1.23 to create memory////////////
*/



}

void init_map(long int i_map, long int n_map )
/*
    COMMENT: Dont forget to delete memory
*/
{		
		if(acs_flag)
			trail_0 = 1. / ( (double) n * (double) nn_tour( ) ) ;
		else
			trail_0 = 1. / ( (rho) * nn_tour() );
		map_init_pheromone_trails(i_map, trail_0 );//initiate phe of each map
		map_compute_total_information(i_map);

		map[i_map].iteration =1;
		map[i_map].convergance_num=1;
		map[i_map].best_so_far_ant->tour_length = INFTY;

		map_start_timers(i_map);

}

void map_allocate_ants(long int i_map)
{
	long int i;
	if((map[i_map].ant = malloc(sizeof( ant_struct ) * n_ants +
			 sizeof(ant_struct *) * n_ants	 )) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
	}

	for ( i = 0 ; i < n_ants ; i++ )
	{
		map[i_map].ant[i].tour        = calloc(n+1, sizeof(long int));
		map[i_map].ant[i].visited     = calloc(n, sizeof(char));
	}

/**********************************************************************/
	if((map[i_map].Mant = malloc(sizeof( ant_struct ) * n_ants +
			 sizeof(ant_struct *) * n_ants	 )) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
	}
	for ( i = 0 ; i < n_ants ; i++ )
	{
		map[i_map].Mant[i].tour        = calloc(n+1, sizeof(long int));
		map[i_map].Mant[i].visited     = calloc(n, sizeof(char));
	}

/**********************************************************************/

	if((map[i_map].best_so_far_ant = malloc(sizeof( ant_struct ) )) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
	}
	map[i_map].best_so_far_ant->tour    = calloc(n+1, sizeof(long int));
	map[i_map].best_so_far_ant->visited = calloc(n, sizeof(char));
//	(*best_so_far_ant).tour        = calloc(n+1, sizeof(long int));
//  (*best_so_far_ant).visited     = calloc(n, sizeof(char));

	if((map[i_map].prob_of_selection = malloc(sizeof( double ) * instance.n_near )) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
	}

}




void end_map(long int s_map)
{
	long int i_map;
	long int i;
	for(i_map=0; i_map<s_map; i_map++)
	{
		free(map[i_map].pheromone);
		free(map[i_map].total);
		free(map[i_map].best_in_try);
		free(map[i_map].time_best_found);
		free(map[i_map].time_total_run);
	    for ( i = 0 ; i < n_ants ; i++ )
		{
		free( map[i_map].ant[i].tour );
		free( map[i_map].ant[i].visited );
		free( map[i_map].Mant[i].tour );
		free( map[i_map].Mant[i].visited );
		}
		free( map[i_map].ant );
        free( map[i_map].Mant );

		free( map[i_map].best_so_far_ant->tour );
		free( map[i_map].best_so_far_ant->visited );
		free( map[i_map].best_so_far_ant );
		free( map[i_map].prob_of_selection );
	}

}