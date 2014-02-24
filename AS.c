#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include"AS.h"
#include"IO.h"
#include"utilities.h"
#include"pheromone_process.h"
#include"Map.h"
#include"GA.h"


double elapsed;

long int n_map;


double *best_in_try;		/* best length of each try */
long int *best_found_at;	/* iteration when best length is found */
double   *time_best_found;	/* time when best length is found  */
double   *time_total_run;	/* total time of a try */
long int iter_to_best;

ant_struct *ant;
ant_struct *Mant;			/*memory ants*/
ant_struct *best_so_far_ant;/*one ant, best so far */
double   *prob_of_selection;/*the probability point for the roullete*/

long int n_ants;			/* number of ants */

long int iteration;         /* iteration counter */
long int max_tries;         /* maximum number of independent tries */
long int max_iteration;
double   max_time;          /* maximal allowed run time of a try  */
long int optimal;           /* optimal solution or bound to find */
long int end_flag;

long int x_num;				/* the individ number of x */ 
long int y_num;				/* the individ number of y */

long int redivide_number;

long int cross_flag;

double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
double rho;           /* parameter for evaporation */
double Q;			  /* pheromone deposition weight */
double q_0;

long int acs_flag=1;

double   trail_0;	/*pheromone initialization */

long int**  tmp_tour;
long int**  pos;

double **pheromone;
double **heuristic;/*簡単化できる！β乗して自分に与える*/
double **total;

double **phex;

double total_elapsed;

long int ras_ranks;				/* additional parameter for rank-based version */
double epsilon;					/* EACO parameter */
char name_buf[LINE_BUF_LEN];	/* directory (include file name) of data file */
struct problem instance;
struct problem cinstance;

long int *sort_ant_no;
long int *sort_Mant_no;
double *save_ant_length;
double *save_Mant_length;

long int dlb_flag = TRUE;
long int ls_flag = 1;
long int evapor_nn_flag = 1;


int main(int argc, char *argv[]) //new main since 2013.9
{
	long int i_map=0;
	long int i=0;
	long int j,m,l;
	long int n_try=0;
	long int swap_flag=0;
	long int cross_max=1;/*very large*/
	long int commu_i=0;
	long int converge_max=200;
	long int converge_max_divided=200;

	start_timers();
	srand( (unsigned)time(NULL) ); /* 乱数の種を設定 */
	seed = (long int) time(NULL);

	init_program(argc, argv);/*inlclude iniate heristic*/



	generate_map(n_map);
	generate_cinstance(n_map);

	indi_map(x_num,y_num);

	set_relation();

	instance.nn_list = compute_nn_lists();

//	printf("commu_max=%ld\n converge_max=%ld\n", commu_max, converge_max);
	fprintf(report,"island_num=%ld\nconverge_max=%ld\ncross_max=%ld\n",n_map,converge_max,cross_max);
	fprintf(best_report,"island_num=%ld\nconverge_max=%ld\ncross_max=%ld\n",n_map,converge_max,cross_max);
	elapsed = elapsed_time();
    printf("Initialization took %.10f seconds\n",elapsed);

	n_map=x_num*y_num+1;


	for ( n_try = 0 ; n_try < max_tries ; n_try++ )
	{
		end_flag=0;
		init_try(n_try);
		/*
		for(i_map=0; i_map<n_map; i_map++)
		{
			set_map(i_map);
			init_map(i_map,n_map );
		}
		*/

		for(i=0; i<cross_max; i++)
		{
			for(i_map=0; i_map<(n_map-1); i_map++)
			{

				set_map(i_map);
				init_map(i_map,n_map);

				
				map[i_map].convergance_num=0;
				while ( !termination_condition(i_map ) ) 
				{
					map_ACS_construct_solutions(i_map);
					map_local_search(i_map);
					map_update_statistics(i_map,n_try);
					map_ACS_pheromone_trail_update(i_map);				
					map[i_map].iteration++;
					if(map[i_map].convergance_num==converge_max_divided)
					{
						printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
						fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
						break;//see as converged, so pause
					}
				}//end while
//				map_report(i_map,n_try);

				for(i=0;i<=n;i++)
				{
					fprintf(report,"%d\t%f\t%f\n",relationship[i_map][map[i_map].best_so_far_ant->tour[i]],instance.nodeptr[map[i_map].best_so_far_ant->tour[i]].x,instance.nodeptr[map[i_map].best_so_far_ant->tour[i]].y);
				}


				if(	map[i_map].best_so_far_ant->tour_length==optimal || elapsed_time() >= max_time)
				{
					end_flag=1;
					fprintf(best_report,"best so far %lf\n",(*best_so_far_ant).tour_length);
					fprintf(best_report,"\n");

					fprintf(try_result,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
							n_try,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
					fflush(report); 

					break;
				}
			}//end i_map
			if(end_flag==1)
			{
				printf("break from cross");
				fprintf(report,"break from cross");
			break;
			}
			fprintf(report,"island cross %ld\n",i);
			printf("island cross %ld\n",i);
//			if(i<1){
//			island_crossover_c(n_map);
//			}
		}//end cross
//////////////////////////////new for 2013 09 19////////////////////////////////


//				n=tn;
//				instance=cinstance;

//				printf("tn: %d n:%d\n",tn, cinstance.n);

				i_map=(n_map-1);

//				instance.nodeptr = read_etsp(name_buf);

				for(i=0;i<tn;i++)
				{
				instance.nodeptr[i]=cinstance.nodeptr[i];
				}

//				printf("%d\n",tn);

				for ( i = 0 ; i < tn ; i++ )
				{
//					printf("\n%ld %ld %lf %lf %ld %ld",i, instance.nodeptr[i].citynumber, instance.nodeptr[i].x, instance.nodeptr[i].y, instance.nodeptr[i].group, instance.nodeptr[i].subcitynumber );
				}



				n=tn;
				map[n_map-1].mn=tn;

				instance.distance = compute_distances();
				instance.nn_list  = compute_nn_lists();

				compute_heuristic();

				///////////////////here starts2014.1.24//////////////////////////////////////

				i=n_map-1;

				map[i].pheromone = generate_double_matrix( n, n );
				map[i].total = generate_double_matrix( n, n );
		
				map[i].best_in_try = calloc(max_tries, sizeof(double));
				map[i].best_found_at = calloc(max_tries, sizeof(long int));
				map[i].time_best_found = calloc(max_tries, sizeof(double));
				map[i].time_total_run = calloc(max_tries, sizeof(double));

				map_allocate_ants(i);

				////////////////here ends////////////////////////////////////////////////

				init_map(i_map,n_map);
				island_crossover_c((n_map-1));
				

				map[i_map].convergance_num=0;
				while ( !map_termination_condition(i_map ) ) 
				{
					map_ACS_construct_solutions(i_map);
					if(cross_flag==1)
					{
					map_local_search(i_map);
					}
					map_update_statistics(i_map,n_try);
					map_ACS_pheromone_trail_update(i_map);				
					map[i_map].iteration++;
					if(map[i_map].convergance_num==converge_max)
					{
						printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
						fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
						break;//see as converged, so pause
					}
				}//end while

				map_report(i_map,n_try);

				for(m=0;m<tn;m++)
				{
					for(l=0;l<tn;l++)
					{
						phex[m][l]=map[i_map].pheromone[m][l];
					}
				}

				end_map(n_map);


				/************redivide starts**********/
				for(j=0 ; j<redivide_number ; j++)
				{

				n_map=rx_num[j]*ry_num[j]+1;

				redivide(rx_num[j],ry_num[j]);

				instance.nn_list = compute_nn_lists();

				for(i=0; i<cross_max; i++)
				{
					for(i_map=0; i_map<(n_map-1); i_map++)
					{

						set_map(i_map);

						/***************initial map**************************/

						if(acs_flag)
							trail_0 = 1. / ( (double) n * (double) nn_tour( ) ) ;
						else
							trail_0 = 1. / ( (rho) * nn_tour() );

						map_init_pheromone_trails(i_map, trail_0 );//initiate phe of each map
						pheromone_divide(i_map);
						map_compute_total_information(i_map);

						map[i_map].iteration =1;
						map[i_map].convergance_num=1;
						map[i_map].best_so_far_ant->tour_length = INFTY;

						map_start_timers(i_map);
						/***************initial map**************************/

				
						map[i_map].convergance_num=0;
						while ( !termination_condition(i_map ) ) 
						{
							map_ACS_construct_solutions(i_map);
//							map_local_search(i_map);
							map_update_statistics(i_map,n_try);
							map_ACS_pheromone_trail_update(i_map);				
							map[i_map].iteration++;
							if(map[i_map].convergance_num==converge_max_divided)
							{
								printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
								fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
								break;//see as converged, so pause
							}
						}//end while
//						map_report(i_map,n_try);
						if(	map[i_map].best_so_far_ant->tour_length==optimal || elapsed_time() >= max_time)
						{
							end_flag=1;
							fprintf(best_report,"best so far %lf\n",(*best_so_far_ant).tour_length);
							fprintf(best_report,"\n");
		
							fprintf(try_result,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
									n_try,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
							fflush(report); 

							break;
						}
					}//end i_map
					if(end_flag==1)
					{
						printf("break from cross");
						fprintf(report,"break from cross");
					break;
					}
					fprintf(report,"island cross %ld\n",i);
					printf("island cross %ld\n",i);
//					if(i<1){
//					island_crossover_c(n_map);
//					}
				}//end cross


				i_map=(n_map-1);

//				instance.nodeptr = read_etsp(name_buf);

				for(i=0;i<tn;i++)
				{
				instance.nodeptr[i]=cinstance.nodeptr[i];
				}

//				printf("%d\n",tn);

				for ( i = 0 ; i < tn ; i++ )
				{
//					printf("\n%ld %ld %lf %lf %ld %ld",i, instance.nodeptr[i].citynumber, instance.nodeptr[i].x, instance.nodeptr[i].y, instance.nodeptr[i].group, instance.nodeptr[i].subcitynumber );
				}

				

				n=tn;
				map[n_map-1].mn=tn;

				instance.distance = compute_distances();
				instance.nn_list  = compute_nn_lists();

				compute_heuristic();

				///////////////////here starts2014.1.24//////////////////////////////////////

				i=(n_map-1);

				map[i].pheromone = generate_double_matrix( n, n );
				map[i].total = generate_double_matrix( n, n );
		
				map[i].best_in_try = calloc(max_tries, sizeof(double));
				map[i].best_found_at = calloc(max_tries, sizeof(long int));
				map[i].time_best_found = calloc(max_tries, sizeof(double));
				map[i].time_total_run = calloc(max_tries, sizeof(double));

				map_allocate_ants(i);

				////////////////here ends////////////////////////////////////////////////

				i_map=(n_map-1);

				init_map(i_map,n_map);
				island_crossover_c((n_map-1));
				

				map[i_map].convergance_num=0;
				while ( !map_termination_condition(i_map ) ) 
				{
					map_ACS_construct_solutions(i_map);
					if(cross_flag==1)
					{
					map_local_search(i_map);
					}
					map_update_statistics(i_map,n_try);
					map_ACS_pheromone_trail_update(i_map);				
					map[i_map].iteration++;
					if(map[i_map].convergance_num==converge_max)
					{
						printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
						fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
						break;//see as converged, so pause
					}
				}//end while

				map_report(i_map,n_try);

				for(m=0;m<tn;m++)
				{
					for(l=0;l<tn;l++)
					{
						phex[m][l]=map[i_map].pheromone[m][l];
					}
				}

				end_map(n_map);

				}
				/************redivide ends************/


				if(	map[i_map].best_so_far_ant->tour_length==optimal || elapsed_time() >= max_time)
				{
					end_flag=1;
					fprintf(best_report,"best so far %lf\n",(*best_so_far_ant).tour_length);
					fprintf(best_report,"\n");

					fprintf(try_result,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
							n_try,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
					fflush(report); 

					break;
				}
//////////////////////////////////end new/////////////////////////////////////

//////////////////////////////new for 2013 09 20////////////////////////////////
	/*

				n=tn;
				instance=cinstance;
				
				island_crossover_c(n_map);

				i_map=0;
				
				map[i_map].convergance_num=0;
				while ( !termination_condition ) 
				{
					ACS_construct_solutions(i_map);
				//	local_search(i_map);
					update_statistics(i_map,n_try);
					ACS_pheromone_trail_update(i_map);				
					iteration++;
				if(map[i_map].convergance_num==converge_max)
					{
						printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
						fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
						break;//see as converged, so pause
					}
				}//end while
				map_report(i_map,n_try);
				if(best_so_far_ant->tour_length==optimal || elapsed_time() >= max_time)
				{
					end_flag=1;
					fprintf(best_report,"best so far %lf\n",(*best_so_far_ant).tour_length);
					fprintf(best_report,"\n");

					fprintf(try_result,"IN TRY %ld:\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
							n_try, best_so_far_ant->tour_length, iter_to_best,total_elapsed);
					fflush(report); 

					break;
				}
				*/
////////////////////////////////////end new/////////////////////////////////////


	}//end try for
//	exit_program();
	end_map(n_map);
	free( sort_ant_no );
	free( sort_Mant_no );
	free( save_ant_length );
	free( save_Mant_length );
    free( instance.distance );
    free( instance.nn_list );
    free( pheromone );
	free( phex );
	free( relationship );
    free( heuristic );
    free( total );
	free(tmp_tour);
	free(pos);
    free( time_best_found );
    free( best_found_at );
    free( best_in_try );
    free( time_total_run );
    for ( i = 0 ; i < n_ants ; i++ )
		{
		free( ant[i].tour );
		free( ant[i].visited );
		free( Mant[i].tour );
		free( Mant[i].visited );
    }
    free( ant );
    free( Mant );
	free( (*best_so_far_ant).tour );
    free( (*best_so_far_ant).visited );
	free( best_so_far_ant );
    free( prob_of_selection );
	return(0);
}




/********************************************************************/
/******************here starts the old main**************************/
/********************************************************************/
/*

int main(int argc, char *argv[]) 
{
	long int i_map=0;
	long int i=0;
	long int n_try=0;
	long int swap_flag=0;
	long int cross_max=200000;//very large
	long int commu_i=0;
	long int converge_max=700;
//	double   total_elapsed=0.0;  //2012.1.4
	start_timers();
	srand( (unsigned)time(NULL) ); // 乱数の種を設定
	seed = (long int) time(NULL);

	init_program(argc, argv);//inlclude iniate heristic
	generate_map(n_map);

	instance.nn_list = compute_nn_lists();

//	printf("commu_max=%ld\n converge_max=%ld\n", commu_max, converge_max);
	fprintf(report,"island_num=%ld\nconverge_max=%ld\ncross_max=%ld\n",n_map,converge_max,cross_max);
	fprintf(best_report,"island_num=%ld\nconverge_max=%ld\ncross_max=%ld\n",n_map,converge_max,cross_max);
	elapsed = elapsed_time();
    printf("Initialization took %.10f seconds\n",elapsed);


	for ( n_try = 0 ; n_try < max_tries ; n_try++ )
	{
		end_flag=0;
		init_try(n_try);
		for(i_map=0; i_map<n_map; i_map++)
		{
			init_map(i_map,n_map );
		}

		for(i=0; i<cross_max; i++)
		{
			for(i_map=0; i_map<n_map; i_map++)
			{
				map[i_map].convergance_num=0;
				while ( !map_termination_condition(i_map ) ) 
				{
					map_ACS_construct_solutions(i_map);
					map_local_search(i_map);
					map_update_statistics(i_map,n_try);
					map_ACS_pheromone_trail_update(i_map);				
					map[i_map].iteration++;
					if(map[i_map].convergance_num==converge_max)
					{
						printf("see as converged!(iteration=%ld)\n",map[i_map].iteration);
						fprintf(report,"see as converged!(iteration=%ld)\n",map[i_map].iteration);
						break;//see as converged, so pause
					}
				}//end while
				map_report(i_map,n_try);
				if(	map[i_map].best_so_far_ant->tour_length==optimal || elapsed_time() >= max_time)
				{
					end_flag=1;
					fprintf(best_report,"best so far %lf\n",(*best_so_far_ant).tour_length);
					fprintf(best_report,"\n");

					fprintf(try_result,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
							n_try,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
					fflush(report); 

					break;
				}
			}//end i_map
			if(end_flag==1)
			{
				printf("break from cross");
				fprintf(report,"break from cross");
			break;
			}
			fprintf(report,"island cross %ld\n",i);
			printf("island cross %ld\n",i);
			//crossover:1.randomly crossover 2.ruellete crossover
			island_crossover(n_map);
		}//end cross
		for(swap_flag=0; swap_flag<commu_max; swap_flag++)
		{
			i_map=swap_flag%2;
			map[i_map].convergance_num=0;
			map_construct_solutions_revised(i_map);
			while ( !map_termination_condition(i_map ) ) 
			{
				map_ACS_construct_solutions(i_map);
				switch(i_map)//i_map
				{
				case 0:
				//if(ls_flag==1)		
					map_local_search(i_map);
					break;
				default:
					//GA
					//map_GA(i_map);
					//map_local_search(i_map);
					map_local_search(i_map);
					//map_2opt_local_search(i_map);
					break;
				}
				map_update_statistics(i_map,n_try);
				map_ACS_pheromone_trail_update(i_map);				
				map[i_map].iteration++;
				if(map[i_map].convergance_num==converge_max)
				{
					printf("see as converged!\n");
					break;//see as converged, so pause
				}
//				iteration++; 
			}//end while
			map_report(i_map,n_try);
			if(
				(map[0].best_so_far_ant->tour_length==optimal
				||map[1].best_so_far_ant->tour_length==optimal)==1  )
				break;
			if(i_map==1)
			{
				printf("swap %ld\n", swap_flag/2+1);
				fprintf(report,"swap %ld\n", swap_flag/2+1);
				swap_pheromone(map[0].pheromone,map[1].pheromone);
//				average_pheromone(map[0].pheromone, map[1].pheromone);
				//average_pheromone(map[0].pheromone,map[1].pheromone)
			}
		}//end swap_flag for

	}//end try for
//	exit_program();
	end_map();
	free( sort_ant_no );
	free( sort_Mant_no );
	free( save_ant_length );
	free( save_Mant_length );
    free( instance.distance );
    free( instance.nn_list );
    free( pheromone );
    free( heuristic );
    free( total );
	free(tmp_tour);
	free(pos);
    free( time_best_found );
    free( best_found_at );
    free( best_in_try );
    free( time_total_run );
    for ( i = 0 ; i < n_ants ; i++ )
		{
		free( ant[i].tour );
		free( ant[i].visited );
		free( Mant[i].tour );
		free( Mant[i].visited );
    }
    free( ant );
    free( Mant );
	free( (*best_so_far_ant).tour );
    free( (*best_so_far_ant).visited );
	free( best_so_far_ant );
    free( prob_of_selection );
	return(0);
}

*/
/********************************************************************/
/******************here starts the old main**************************/
/********************************************************************/




//		construct_solutions_revised();//when EACO, commenting it will be better?
//		while ( !termination_condition() ) 
//		{
//			memory_construct_solutions();
//			construct_solutions();
//			construct_solutions_revised();
//			ACS_construct_solutions();
//			EACS_construct_solutions();
//			M_ACS_construct_solutions();
//			EACO_construct_solutions();
//			M_EACO_construct_solutions();
/*******************************************************/
/*			for ( i = 0 ; i < n_ants ; i++ )
			{
				sort_ant_no[ i ]= i ;
				save_ant_length[ i ]=ant[ i ].tour_length;
				sort_Mant_no[ i ]= i ;
				save_Mant_length[ i ]=Mant[ i ].tour_length;
			}
			sort2(save_ant_length, sort_ant_no, 0, n_ants-1);
			sort2(save_Mant_length, sort_Mant_no, 0, n_ants-1);
			for(i=0; i<10; i++)
			{
				rnd_ant_no=(rand()%(n_ants-1)+1);// random number from 1~n-1 
				copy_from_to( &Mant[ sort_Mant_no[i] ], &ant[ sort_ant_no[rnd_ant_no] ] );
//				copy_from_to(&Mant[ sort_Mant_no[i] ], &ant[ sort_ant_no[ (n_ants-1)-i ] ]);
			}
*/
/*******************************************************/

			/***deposit before applying local search****/
			
//			update_statistics(n_try);
//			ACS_pheromone_trail_update();  
//			global_acs_pheromone_update( best_so_far_ant );

			/*******************************************/


/*
			if(ls_flag==1)		
				local_search();
			update_statistics(n_try);
			ACS_pheromone_trail_update();  
*/
//			EACO_pheromone_trail_update();

//			global_acs_pheromone_update( best_so_far_ant );

//			pheromone_trail_update();  


			
//			EACS_pheromone_trail_update();  
//			iteration++; 
//		}//end while





void init_program( long int argc, char *argv[] ) 
{

	set_default_parameters();

	best_in_try = calloc(max_tries, sizeof(double));
	best_found_at = calloc(max_tries, sizeof(long int));
	time_best_found = calloc(max_tries, sizeof(double));
	time_total_run = calloc(max_tries, sizeof(double));

	sort_ant_no = calloc(n_ants, sizeof(long int));
	sort_Mant_no = calloc(n_ants, sizeof(long int));
	save_ant_length = calloc(n_ants, sizeof(double));
	save_Mant_length = calloc(n_ants, sizeof(double));


	instance.nodeptr = read_etsp(name_buf);

	cinstance.nodeptr = read_etsp(name_buf);

	report = fopen(instance.name, "w");
	best_report = fopen("best", "w");
	try_result = fopen("try_result", "w");

	instance.distance = compute_distances();

	cinstance.distance = compute_distances();

	tmp_tour   = generate_int_matrix( n_ants, n+1);/*n comes from read_etsp()*/
	pos		   = generate_int_matrix( n_ants, n);

	allocate_ants();	

  	pheromone = generate_double_matrix( n, n );

	phex=generate_double_matrix( n, n );

	relationship = generate_int_matrix( n, n );

	heuristic = generate_double_matrix( n, n );
	
	total = generate_double_matrix( n, n );
	compute_heuristic();
	write_params( ); 
}


void set_default_parameters() 
/*    
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{

    n_ants         = 30;    /* number of ants */
    max_tries      = 20;

	max_iteration  = 800;

    max_time       = 3600;

    alpha          = 1.0;
    beta           = 5.0;
    rho            = 0.675;
    Q			   = 1; 	/* pheromone deposition weight */
	q_0			   = 0;

	epsilon		   = 1;		/*set it to 1 if ras_randks*/
	
	ras_ranks	   = 6;
	instance.n_near= 20;

	ls_flag        = 1; /*if 1, use local search, and only nn deposition */
	evapor_nn_flag = 1;
	acs_flag	   = 1; /* traio_0 changes when acs_flag */
	cross_flag		=1;

	x_num			=3;
	y_num			=3;

	n_map			=11;

	redivide_number	=3;

	rx_num = calloc(redivide_number, sizeof(long int));
	ry_num = calloc(redivide_number, sizeof(long int));

	rx_num[0]		=3;
	ry_num[0]		=3;

	rx_num[1]		=3;
	ry_num[1]		=3;

	rx_num[2]		=3;
	ry_num[2]		=3;

//	rx_num[3]		=12;
//	ry_num[3]		=8;

//	rx_num[4]		=8;
//	ry_num[4]		=12;
	




//#define rat783
//#define pcb442
//#define pcb1173
#define tsp225
//#define pr2392
//#define eil51
//#define berlin52
//#define oliver30
//#define st70
//#define eil76
//#define pr76
//#define rat99
//#define kroA100
//#define kroB100
//#define kroC100
//#define kroD100
//#define kroE100
//#define rd100
//#define eil101
//#define lin105
//#define pr107
//#define pr124
//#define bier127
//#define ch130
//#define pr136
//#define pr144
//#define ch150
//#define kroA150
//#define kroB150
//#define pr152
//#define u159
//#define rat195
//#define d198
//#define kroA200
//#define kroB200
//#define ts225
//
//#define pr226
//#define gil262
//#define pr264
//#define a280
//#define pr299
//#define lin318
//#define linhp318
//#define rd400
//#define fl417
//#define pr439
//#define pcb442
//#define d493
//#define u574
//#define rat575
//#define p654
//#define d657
//#define u724
//#define rat783
//#define pr1002
//#define u1060
//#define vm1084
//#define pcb1173
//#define d1291
//#define rl1323
//#define rl1304
//#define nrw1379
//#define fl1400
//#define u1432
//#define fl1577 
//#define d1655
//#define vm1748
//#define u1817
//#define rl1889
//#define d2103
//#define u2152
//#define u2319
//#define pr2392
//#define pcb3038


#ifdef eil51
optimal = 426;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\eil51.tsp", LINE_BUF_LEN);
#endif

#ifdef berlin52
optimal = 7542;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\berlin52.tsp", LINE_BUF_LEN);
#endif

#ifdef	oliver30
optimal = 420;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\oliver30.tsp", LINE_BUF_LEN);
#endif

#ifdef	st70
optimal = 675;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\st70.tsp", LINE_BUF_LEN);
#endif

#ifdef	eil76
optimal = 538;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\eil76.tsp", LINE_BUF_LEN);
#endif

#ifdef	pr76
optimal = 108159;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr76.tsp", LINE_BUF_LEN);
#endif

#ifdef	rat99
optimal = 1211;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rat99.tsp", LINE_BUF_LEN);
#endif

#ifdef kroA100
optimal = 21282;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroA100.tsp", LINE_BUF_LEN);
#endif

#ifdef kroB100
optimal = 22141;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroB100.tsp", LINE_BUF_LEN);
#endif

#ifdef kroC100
optimal = 20749;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroC100.tsp", LINE_BUF_LEN);
#endif

#ifdef kroD100
optimal = 21294;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroD100.tsp", LINE_BUF_LEN);
#endif

#ifdef kroE100
optimal = 22068;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroE100.tsp", LINE_BUF_LEN);
#endif

#ifdef	rd100
optimal = 7910;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rd100.tsp", LINE_BUF_LEN);
#endif

#ifdef	eil101
optimal = 629;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\eil101.tsp", LINE_BUF_LEN);
#endif

#ifdef	lin105
optimal = 14379;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\lin105.tsp", LINE_BUF_LEN);
#endif

#ifdef	pr107
optimal = 44303;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr107.tsp", LINE_BUF_LEN);
#endif

#ifdef	pr124
optimal = 59030;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr124.tsp", LINE_BUF_LEN);
#endif

#ifdef bier127
optimal = 118282;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\bier127.tsp", LINE_BUF_LEN);
#endif

#ifdef ch130
optimal = 6110;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\ch130.tsp", LINE_BUF_LEN);
#endif

#ifdef pr136
optimal = 96772;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr136.tsp", LINE_BUF_LEN);
#endif

#ifdef pr144
optimal = 58537;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr144.tsp", LINE_BUF_LEN);
#endif

#ifdef ch150
optimal = 6528;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\ch150.tsp", LINE_BUF_LEN);
#endif

#ifdef kroA150
optimal = 26524;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroA150.tsp", LINE_BUF_LEN);
#endif

#ifdef kroB150
optimal = 26130;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroB150.tsp", LINE_BUF_LEN);
#endif

#ifdef pr152
optimal = 73682;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr152.tsp", LINE_BUF_LEN);
#endif

#ifdef u159
optimal = 42080;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u159.tsp", LINE_BUF_LEN);
#endif

#ifdef rat195
optimal = 2323;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rat195.tsp", LINE_BUF_LEN);
#endif

#ifdef d198
optimal = 15780;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d198.tsp", LINE_BUF_LEN);
#endif

#ifdef kroA200
optimal = 29368;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroA200.tsp", LINE_BUF_LEN);
#endif

#ifdef kroB200
optimal = 29437;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\kroB200.tsp", LINE_BUF_LEN);
#endif

#ifdef ts225
optimal = 126643;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\ts225.tsp", LINE_BUF_LEN);
#endif

#ifdef tsp225
//optimal = 3916;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\tsp225.tsp", LINE_BUF_LEN);
#endif

#ifdef pr226
optimal = 80369;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr226.tsp", LINE_BUF_LEN);
#endif

#ifdef gil262
optimal = 2378;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\gil262.tsp", LINE_BUF_LEN);
#endif

#ifdef pr264
optimal = 49135;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr264.tsp", LINE_BUF_LEN);
#endif

#ifdef a280
//optimal = 2579;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\a280.tsp", LINE_BUF_LEN);
#endif

#ifdef pr299
optimal = 48191;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr299.tsp", LINE_BUF_LEN);
#endif

#ifdef lin318
optimal = 42029;	//lin318
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\lin318.tsp", LINE_BUF_LEN);
#endif

#ifdef linhp318
optimal = 41345;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\linhp318.tsp", LINE_BUF_LEN);
#endif

#ifdef rd400
optimal = 15281;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rd400.tsp", LINE_BUF_LEN);
#endif

#ifdef fl417
optimal = 11861;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\fl417.tsp", LINE_BUF_LEN);
#endif

#ifdef pr439
optimal = 107217;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr439.tsp", LINE_BUF_LEN);
#endif

#ifdef pcb442
//optimal = 50778;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pcb442.tsp", LINE_BUF_LEN);
#endif

#ifdef d493
optimal = 35002;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d493.tsp", LINE_BUF_LEN);
#endif

#ifdef u574
optimal = 36905;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u574.tsp", LINE_BUF_LEN);
#endif

#ifdef rat575
optimal = 6773;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rat575.tsp", LINE_BUF_LEN);
#endif

#ifdef p654
optimal = 34643;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\p654.tsp", LINE_BUF_LEN);
#endif

#ifdef d657
optimal = 48912;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d657.tsp", LINE_BUF_LEN);
#endif

#ifdef u724
optimal = 41910;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u724.tsp", LINE_BUF_LEN);
#endif

#ifdef rat783
//optimal = 8806;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rat783.tsp", LINE_BUF_LEN);
#endif

#ifdef pr1002
optimal = 259045;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr1002.tsp", LINE_BUF_LEN);
#endif

#ifdef u1060
optimal = 224094;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u1060.tsp", LINE_BUF_LEN);
#endif

#ifdef vm1084
optimal = 239297;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\vm1084.tsp", LINE_BUF_LEN);
#endif


#ifdef pcb1173
//optimal = 56892;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pcb1173.tsp", LINE_BUF_LEN);
#endif

#ifdef d1291
optimal = 50801;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d1291.tsp", LINE_BUF_LEN);
#endif

#ifdef rl1304
optimal = 252948;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rl1304.tsp", LINE_BUF_LEN);
#endif

#ifdef rl1323
optimal = 270199;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rl1323.tsp", LINE_BUF_LEN);
#endif

#ifdef nrw1379
optimal = 56638;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\nrw1379.tsp", LINE_BUF_LEN);
#endif

#ifdef fl1400
optimal = 20127;
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\fl1400.tsp", LINE_BUF_LEN);
#endif

#ifdef u1432
optimal = 152970;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u1432.tsp", LINE_BUF_LEN);
#endif

#ifdef fl1577
optimal = 22249;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\fl1577.tsp", LINE_BUF_LEN);
#endif

#ifdef d1655
optimal = 62128;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d1655.tsp", LINE_BUF_LEN);
#endif

#ifdef vm1748 
optimal = 336556;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\vm1748.tsp", LINE_BUF_LEN);
#endif

#ifdef u1817 
optimal = 57201;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u1817.tsp", LINE_BUF_LEN);
#endif

#ifdef rl1889 
optimal = 316536;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\rl1889.tsp", LINE_BUF_LEN);
#endif

#ifdef d2103 
optimal = 80450;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\d2103.tsp", LINE_BUF_LEN);
#endif

#ifdef u2152 
optimal = 64253;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u2152.tsp", LINE_BUF_LEN);
#endif

#ifdef u2319
optimal = 234256;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\u2319.tsp", LINE_BUF_LEN);
#endif

#ifdef pr2392
//optimal = 378032;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pr2392.tsp", LINE_BUF_LEN);
#endif

#ifdef pcb3038
optimal = 137694;	
strncpy( name_buf,"C:\\islandsTSP\\tsp_data\\pcb3038.tsp", LINE_BUF_LEN);
#endif


}


long int ** compute_nn_lists( void )	/*problem maybe exist...*/
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    long int i, node, nn;
    double *distance_vector;
    long int *help_vector;
    long int **m_nnear;
 
    TRACE ( printf("\n computing nearest neighbor lists, "); )

	nn = instance.n_near;
    if((m_nnear = malloc(sizeof(long int) * n * nn
			     + n * sizeof(long int *))) == NULL){
	exit(EXIT_FAILURE);
    }
    distance_vector = calloc(n, sizeof(double));
    help_vector = calloc(n, sizeof(long int));
 
    for ( node = 0 ; node < n ; node++ )
	{/* compute cnd-sets for all node */
		m_nnear[node] = (long int *)(m_nnear + n) + node * nn;

		for ( i = 0 ; i < n ; i++ )
		{  /* Copy distances from nodes to the others */
			distance_vector[i] = instance.distance[node][i];
			help_vector[i] = i;
		}
		distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
		sort2(distance_vector, help_vector, 0, n-1);
		for ( i = 0 ; i < nn ; i++ ) 
		{
			m_nnear[node][i] = help_vector[i];
		}
    }
    free(distance_vector);
    free(help_vector);
    TRACE ( printf("\n    .. done\n"); )
    return m_nnear;
}



double ** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
    long int     i, j;
    double     **matrix;

    if((matrix = malloc(sizeof(double) * n * n +
			sizeof(double *) * n	 )) == NULL){
	fprintf(stderr,"Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n ; i++ ) {
	matrix[i] = (double *)(matrix + n) + i*n;
	for ( j = 0  ; j < n ; j++ ) {
	    matrix[i][j] = distance(i, j);
	}
    }
    return matrix;
}

void allocate_ants ( void )
/*    
      FUNCTION:       allocate the memory for the ant colony, the best-so-far and 
                      the iteration best ant
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that 
                      store intermediate tours

*/
{
    long int i;
  
    if((ant = malloc(sizeof( ant_struct ) * n_ants +
		     sizeof(ant_struct *) * n_ants	 )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
	ant[i].tour        = calloc(n+1, sizeof(long int));
	ant[i].visited     = calloc(n, sizeof(char));
    }

/**********************************************************************/
    if((Mant = malloc(sizeof( ant_struct ) * n_ants +
		     sizeof(ant_struct *) * n_ants	 )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
	Mant[i].tour        = calloc(n+1, sizeof(long int));
	Mant[i].visited     = calloc(n, sizeof(char));
    }
/**********************************************************************/

    if((best_so_far_ant = malloc(sizeof( ant_struct ) )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
    (*best_so_far_ant).tour        = calloc(n+1, sizeof(long int));
    (*best_so_far_ant).visited     = calloc(n, sizeof(char));

    if((prob_of_selection = malloc(sizeof( double ) * instance.n_near )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
}


void init_try( long int ntry ) 
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
{ 
	start_timers();
 //   elapsed = elapsed_time();

    iteration    = 1;         
    (*best_so_far_ant).tour_length = INFTY;

	if(acs_flag)
		trail_0 = 1. / ( (double) n * (double) nn_tour( ) ) ;
	else
		trail_0 = 1. / ( (rho) * nn_tour() );

	init_pheromone_trails( trail_0 );
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();    
//    fprintf(report,"begin try %li \n",ntry);


}


double nn_tour( void )
/*    
      FUNCTION:       generate some nearest neighbor tour and compute tour length
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  needs ant colony and one statistic ants
*/
{
    long int phase;
	double   help;

    ant_empty_memory( &ant[0] );

    phase = 0; /* counter of the construction steps */
    place_ant( &ant[0], phase);

    while ( phase < n-1 ) {
	phase++;
	choose_closest_next( &ant[0],phase);
    }
    phase = n;
    ant[0].tour[n] = ant[0].tour[0];
/*   copy_from_to( &ant[0], best_so_far_ant ); */
    ant[0].tour_length = compute_tour_length( ant[0].tour );
    help = ant[0].tour_length;
    ant_empty_memory( &ant[0] );
    return help;
}

void init_pheromone_trails( double initial_trail )
/*    
      FUNCTION:      initialize pheromone trails
      INPUT:         initial value of pheromone trails "initial_trail"
      OUTPUT:        none
      (SIDE)EFFECTS: pheromone matrix is reinitialized
*/
{
    long int i, j;

    TRACE ( printf(" init trails with %.15f\n",initial_trail); );

    /* Initialize pheromone trails */
    for ( i = 0 ; i < n ; i++ ) {
	for ( j =0 ; j <= i ; j++ ) {
	    pheromone[i][j] = initial_trail;
	    pheromone[j][i] = initial_trail;
	    total[i][j] = initial_trail;
	    total[j][i] = initial_trail;
	}
    }
}

void map_init_pheromone_trails(long int i_map, double initial_trail )
/*    
      FUNCTION:      initialize pheromone trails
      INPUT:         initial value of pheromone trails "initial_trail"
      OUTPUT:        none
      (SIDE)EFFECTS: pheromone matrix is reinitialized
*/
{
    long int i, j;
    /* Initialize pheromone trails */
    for ( i = 0 ; i < n ; i++ )
	{
		for ( j =0 ; j <= i ; j++ )
		{
			map[i_map].pheromone[i][j] = initial_trail;
			map[i_map].pheromone[j][i] = initial_trail;
			map[i_map].total[i][j] = initial_trail;//set arbitartily
			map[i_map].total[j][i] = initial_trail;
		}
    }
}

void compute_total_information( void )
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{
    long int     i, j;

    TRACE ( printf("compute total information\n"); );

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
	    total[i][j] = pow(pheromone[i][j],alpha) * pow(heuristic[i][j],beta);
	    total[j][i] = total[i][j];
	}
    }
}


void map_compute_total_information( long int i_map )
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{
    long int     i, j;

    TRACE ( printf("compute total information\n"); );

    for ( i = 0 ; i <  n ; i++ ) {
	for ( j = 0 ; j <= i ; j++ ) {
	    map[i_map].total[i][j] = pow(map[i_map].pheromone[i][j],alpha) * pow(heuristic[i][j],beta);
	    map[i_map].total[j][i] = map[i_map].total[i][j];
//		printf("total: %d,%d  %0.15f\n",i,j,map[i_map].total[i][j]);
	}
    }
}

void compute_heuristic()
{
	long int i, j;

    for ( i = 0 ; i < n ; i++ )
	{
		for ( j =0 ; j <= i ; j++ )
		{
			heuristic[i][j] = (1.0 / ((double) instance.distance[i][j] + 0.1));
			heuristic[j][i] = heuristic[i][j];
		}
    }
}

void place_ant( ant_struct *a , long int step )
/*    
      FUNCTION:      place an ant on a randomly chosen initial city
      INPUT:         pointer to ant and the number of construction steps 
      OUTPUT:        none
      (SIDE)EFFECT:  ant is put on the chosen city
*/
{
    long int     rnd;

    rnd = (long int) (ran01( &seed ) * (double) n); /* random number between 0 .. n-1 */
    (*a).tour[step] = rnd; 
    (*a).visited[rnd] = TRUE;
}

void choose_closest_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      Choose-
	  for an ant the closest city as the next one 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
	double min_distance;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    min_distance = INFTY;             /* Search shortest edge */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( (*a).visited[city] ) 
	    ; /* city already visited */
	else {
	    if ( instance.distance[current_city][city] < min_distance) {
		next_city = city;
		min_distance = instance.distance[current_city][city];
	    }
	} 
    }
    DEBUG( assert ( 0 <= next_city && next_city < n); );
    (*a).tour[phase] = next_city;
    (*a).visited[next_city] = TRUE;
}

void ant_empty_memory( ant_struct *a ) 
/*    
      FUNCTION:       empty the ants's memory regarding visited cities
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited cities is reinitialized to FALSE
*/
{
    long int   i;

    for( i = 0 ; i < n ; i++ ) {
	(*a).visited[i]=FALSE;
//	(*a).tour[i]=FALSE;
    }
}

long int map_termination_condition( long int i_map )
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return ( ((iteration >= max_iteration) && (elapsed_time() >= max_time)) || 
	  (map[i_map].best_so_far_ant->tour_length <= optimal)); 

/******************************stable state of the tour************************/


/******************************stable state of the tour************************/

}

long int termination_condition( void )
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return ( ((iteration >= max_iteration) && (elapsed_time() >= max_time)) || 
	  ((*best_so_far_ant).tour_length <= optimal)); 

/******************************stable state of the tour************************/


/******************************stable state of the tour************************/

}
void map_local_acs_pheromone_update(long int i_map, ant_struct *a, long int phase )
/*    
      FUNCTION:      removes some pheromone on edge just passed by the ant
      INPUT:         pointer to ant and number of constr. phase
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
      COMMENTS:      I did not do experiments with with different values of the parameter 
                     xi for the local pheromone update; therefore, here xi is fixed to 0.1 
		     as suggested by Gambardella and Dorigo for the TSP. If you wish to run 
		     experiments with that parameter it may be reasonable to use it as a 
		     commandline parameter
*/
{  
    long int  h, j;
 
	j = (*a).tour[phase];
    h = (*a).tour[phase-1];
	/* still additional parameter has to be introduced */
	map[i_map].pheromone[h][j] = (1. - 0.1) * map[i_map].pheromone[h][j] + 0.1 * trail_0;
	map[i_map].pheromone[j][h] = map[i_map].pheromone[h][j];
    map[i_map].total[h][j] = pow(map[i_map].pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
    map[i_map].total[j][h] = map[i_map].total[h][j];
}

void local_acs_pheromone_update( ant_struct *a, long int phase )
/*    
      FUNCTION:      removes some pheromone on edge just passed by the ant
      INPUT:         pointer to ant and number of constr. phase
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
      COMMENTS:      I did not do experiments with with different values of the parameter 
                     xi for the local pheromone update; therefore, here xi is fixed to 0.1 
		     as suggested by Gambardella and Dorigo for the TSP. If you wish to run 
		     experiments with that parameter it may be reasonable to use it as a 
		     commandline parameter
*/
{  
    long int  h, j;
 
	j = (*a).tour[phase];
    h = (*a).tour[phase-1];
	/* still additional parameter has to be introduced */
	pheromone[h][j] = (1. - 0.1) * pheromone[h][j] + 0.1 * trail_0;
    pheromone[j][h] = pheromone[h][j];
    total[h][j] = pow(pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
    total[j][h] = total[h][j];
}


void construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */
    for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
    }
    
    step = 0; 
    /* Place the ants on same initial city */
    for ( k = 0 ; k < n_ants ; k++ )
		place_ant( &ant[k], step);

    while ( step < n-1 )
	{
		step++;
		for ( k = 0 ; k < n_ants ; k++ )
		{
		   move_to_next_normalized( &ant[k], step);
		}
    }

    step = n;
    for ( k = 0 ; k < n_ants ; k++ )
	{
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );
	}
}



void construct_solutions_revised( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */

 for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
	    step = 0; 
		place_ant( &ant[k], step);
		while ( step < n-1 )
			{
				step++;
				move_to_next_normalized( &ant[k], step);
			}

        step = n;
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );
    }
}

void map_construct_solutions_revised( long int i_map )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */

 for ( k = 0 ; k < n_ants ; k++)
	{
//		ant_empty_memory( &ant[k] );
		ant_empty_memory( &map[i_map].ant[k] );
		step = 0; 
//		place_ant( &ant[k], step);
		place_ant( &map[i_map].ant[k], step);
		while ( step < n-1 )
			{
				step++;
				move_to_next_normalized( &map[i_map].ant[k], step);
			}

        step = n;
		map[i_map].ant[k].tour[n] = map[i_map].ant[k].tour[0];
		map[i_map].ant[k].tour_length = compute_tour_length( map[i_map].ant[k].tour );
    }
}


void map_ACS_construct_solutions(long int i_map)
{
	long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */

	for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &map[i_map].ant[k] );
	    step = 0; 
		place_ant( &map[i_map].ant[k], step);
//		printf("imap:%d ant k: %d step: %d",i_map,k,step);
		while ( step < n-1 )
			{
				step++;
				map_ACS_move_to_next(i_map, &map[i_map].ant[k], step);//
				map_local_acs_pheromone_update(i_map, &map[i_map].ant[k], step );//local evaporation
			}

        step = n;
		map[i_map].ant[k].tour[n] = map[i_map].ant[k].tour[0];
		map[i_map].ant[k].tour_length = compute_tour_length( map[i_map].ant[k].tour );
    }
}

void ACS_construct_solutions()
{
	long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */

	for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
	    step = 0; 
		place_ant( &ant[k], step);
		while ( step < n-1 )
			{
				step++;
				ACS_move_to_next( &ant[k], step);//
				local_acs_pheromone_update( &ant[k], step );//
			}

        step = n;
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );
    }
}
void EACS_construct_solutions()
{
	long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */
	double	tmp_length;
	long int *tem_tour;

	tem_tour=calloc(n+1,sizeof(long int));

    /* Mark all cities as unvisited */

	for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
		/***************** add tmp phe, BEGIN *****************/
		tmp_length=ant[k].tour_length;
		copy_from_tour_to_tour(ant[k].tour, tem_tour);
		add_tmp_phe( &ant[k] );
		/***************** add tmp phe, END *******************/

	    step = 0; 
		place_ant( &ant[k], step);
		while ( step < n-1 )
			{
				step++;
				ACS_move_to_next( &ant[k], step);//
//				local_acs_pheromone_update( &ant[k], step );//
			}

        step = n;
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );
		/***************** minus tmp phe, BEGIN *****************/
		minus_tmp_phe_tour( tem_tour, tmp_length);
		/***************** minus tmp phe, END *******************/

    }
}
void M_ACS_construct_solutions()
{
	long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

	long int *tem_tour;

	long int i; 
	long int pos_current;
	long int pos_new;
	double length;
    empty_matrix_memory( tmp_tour, n_ants, n+1); 	
    empty_matrix_memory( pos, n_ants, n);
	tem_tour=calloc(n+1,sizeof(long int));

	/* Mark all cities as unvisited */
	for ( k = 0 ; k < n_ants ; k++)
	{
		copy_ant_to_tour( &ant[k], tmp_tour[k]); /*last tour, temparary tour*/
		for ( i = 0 ; i < n ; i++ )
		{
			pos[k][ tmp_tour[k][i] ] = i;/* position */
		}
		copy_from_tour_to_tour(ant[k].tour, tem_tour);
		ant_empty_memory( &Mant[k] );
	    step = 0; 
		place_ant( &Mant[k], step);
		while ( step < n-1 )
		{
				step++;
				ACS_move_to_next( &Mant[k], step);//
				local_acs_pheromone_update( &Mant[k], step );//
				pos_current = pos[k][ Mant[k].tour[step-1] ] + 1;
				if(Mant[k].tour[step] != tmp_tour[k][ pos_current ] )
				{
					pos_new = pos[k][ Mant[k].tour[step] ];
					//printTour(ant[k].tour);
					//printTour(tmp_tour);
					//printNchar(ant[k].visited);
					if(pos_current%n==0)
					{
						pos_current=0;
						swap( tmp_tour[k], pos_current, pos_new );
						swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
						tmp_tour[k][n]=tmp_tour[k][0];
					}
					else
					{
						swap( tmp_tour[k], pos_current, pos_new );
						swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
					}
					//visited_change(&ant[k], step, tmp_tour[pos_current], tmp_tour[pos_new] );
					//printTour(ant[k].tour);
					//printTour(tmp_tour);
					//printNchar(ant[k].visited);
					length = compute_tour_length( tmp_tour[k] );
					if( length < ant[k].tour_length )
					{
						copy_tour_to_ant(tmp_tour[k], &Mant[k]);/*firstly, save it into Mant*/
						Mant[k].tour_length=length;
						break;
					}
					else
						continue;
				}
				else
					continue;	
			}/*end while*/

        step = n;
		Mant[k].tour[n] = Mant[k].tour[0];
		Mant[k].tour_length = compute_tour_length( Mant[k].tour );

		if(Mant[k].tour_length<ant[k].tour_length)
			copy_from_to( &Mant[ k ], &ant[ k ] );
    }/*end for*/
}
void EACO_construct_solutions( void )
/*    
      FUNCTION:       
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */
	double	tmp_length;
	long int *tem_tour;
    /* Mark all cities as unvisited */
	tem_tour=calloc(n+1,sizeof(long int));
	for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
		/***************** add tmp phe, BEGIN *****************/
		tmp_length=ant[k].tour_length;
		copy_from_tour_to_tour(ant[k].tour, tem_tour);
		add_tmp_phe( &ant[k] );
		/***************** add tmp phe, END *******************/
	    step = 0; 
		place_ant( &ant[k], step);
		while ( step < n-1 )
			{
				step++;
				move_to_next_normalized( &ant[k], step);
			}

        step = n;
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );
		/***************** minus tmp phe, BEGIN *****************/
		minus_tmp_phe( &ant[k], tmp_length);/*this is better, but not reasonalbe?*/
//		minus_tmp_phe_tour(tem_tour , tmp_length);
//		add_tmp_phe( &ant[k] );
		/***************** minus tmp phe, END *******************/
    }
}

void M_EACO_construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */
	
	double	tmp_length;
	long int *tem_tour;

	long int i; 
	long int pos_current;
	long int pos_new;
	double length;
    empty_matrix_memory( tmp_tour, n_ants, n+1); 	
    empty_matrix_memory( pos, n_ants, n);

	tem_tour=calloc(n+1,sizeof(long int));

    /* Mark all cities as unvisited */

 for ( k = 0 ; k < n_ants ; k++)
	{
		copy_ant_to_tour( &ant[k], tmp_tour[k]); /*last tour, temparary tour*/
		for ( i = 0 ; i < n ; i++ )
		{
			pos[k][ tmp_tour[k][i] ] = i;/* position */
		}
		/***************** add tmp phe, BEGIN *****************/
		tmp_length=ant[k].tour_length;
		copy_from_tour_to_tour(ant[k].tour, tem_tour);
		add_tmp_phe( &ant[k] );
		/***************** add tmp phe, END *******************/
		ant_empty_memory( &Mant[k] );

		step = 0; 
		place_ant( &Mant[k], step);
		while ( step < n-1 )
			{
				step++;
				move_to_next_normalized( &Mant[k], step);
				pos_current = pos[k][ Mant[k].tour[step-1] ] + 1;

				if(Mant[k].tour[step] != tmp_tour[k][ pos_current ] )
				{
					pos_new = pos[k][ Mant[k].tour[step] ];
					//printTour(ant[k].tour);
					//printTour(tmp_tour);
					//printNchar(ant[k].visited);
					if(pos_current%n==0)
					{
						pos_current=0;
						swap( tmp_tour[k], pos_current, pos_new );
						swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
						tmp_tour[k][n]=tmp_tour[k][0];
					}
					else
					{
						swap( tmp_tour[k], pos_current, pos_new );
						swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
					}
					//visited_change(&ant[k], step, tmp_tour[pos_current], tmp_tour[pos_new] );
					//printTour(ant[k].tour);
					//printTour(tmp_tour);
					//printNchar(ant[k].visited);
					length = compute_tour_length( tmp_tour[k] );
					if( length < ant[k].tour_length )
					{
						copy_tour_to_ant(tmp_tour[k], &Mant[k]);/*firstly, save it into Mant*/
						Mant[k].tour_length=length;
						break;
					}
					else
						continue;
				}
				else
					continue;	
			}

        step = n;
		Mant[k].tour[n] = Mant[k].tour[0];
		Mant[k].tour_length = compute_tour_length( Mant[k].tour );

		if(Mant[k].tour_length<ant[k].tour_length)
			copy_from_to( &Mant[ k ], &ant[ k ] );

		/***************** minus tmp phe, BEGIN *****************/
		minus_tmp_phe( &ant[k], tmp_length);
//		minus_tmp_phe_tour(tem_tour , tmp_length);
//		add_tmp_phe( &ant[k] );
		/***************** minus tmp phe, END *******************/
    }
	free(tem_tour);
}
void memory_construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

	long int i; 
	long int pos_current;
	long int pos_new;
	double length;
    empty_matrix_memory( tmp_tour, n_ants, n+1); 	
    empty_matrix_memory( pos, n_ants, n);
	/* Mark all cities as unvisited */
    for ( k = 0 ; k < n_ants ; k++)
	{
//		copy_from_to( &ant[ k ], &Mant[ k ] );
//		Mant[ k ].tour_length = ant[ k ].tour_length;
		copy_ant_to_tour( &ant[k], tmp_tour[k]); /*last tour, temparary tour*/
		for ( i = 0 ; i < n ; i++ )
		{
			pos[k][ tmp_tour[k][i] ] = i;/* position */
		}
		ant_empty_memory( &Mant[k] );
		step = 0;
		place_ant( &Mant[k], step);
		while ( step < n-1 )
		{
			step++;
			move_to_next_normalized( &Mant[k], step);
			pos_current = pos[k][ Mant[k].tour[step-1] ] + 1;
/***************** ↓　BEGIN　↓ ********************/
//   		printf("step%ld: %ld\n",step, ant[k].tour[step]);
			if(Mant[k].tour[step] != tmp_tour[k][ pos_current ] )
			{
				pos_new = pos[k][ Mant[k].tour[step] ];
				//printTour(ant[k].tour);
				//printTour(tmp_tour);
				//printNchar(ant[k].visited);
				if(pos_current%n==0)
				{
					pos_current=0;
					swap( tmp_tour[k], pos_current, pos_new );
					swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
					tmp_tour[k][n]=tmp_tour[k][0];
				}
				else
				{
					swap( tmp_tour[k], pos_current, pos_new );
					swap( pos[k], tmp_tour[k][pos_current], tmp_tour[k][pos_new] );
				}
				//visited_change(&ant[k], step, tmp_tour[pos_current], tmp_tour[pos_new] );
				//printTour(ant[k].tour);
				//printTour(tmp_tour);
				//printNchar(ant[k].visited);
				length = compute_tour_length( tmp_tour[k] );
				if( length < ant[k].tour_length )
				{
					copy_tour_to_ant(tmp_tour[k], &Mant[k]);/*firstly, save it into Mant*/
					Mant[k].tour_length=length;
					break;
				}
				else
					continue;
			}
			else
				continue;	
		}/*end while*/
		step = n;
		Mant[k].tour[n] = Mant[k].tour[0];
		Mant[k].tour_length = compute_tour_length( Mant[k].tour );

		if(Mant[k].tour_length<ant[k].tour_length)
			copy_from_to( &Mant[ k ], &ant[ k ] );
	
//		checkTour( Mant[ k ].tour );
    }/*end for*/
}  

void move_to_next( ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
	double tmp_rnd;/*look rnd*/

    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0.0, sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else
	{

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		tmp_rnd = rnd;
		rnd *= sum_prob;
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum < rnd )//
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}

void map_move_to_next_normalized(long int i_map, ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 

    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0.0, sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = map[i_map].total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];/* prepare to calculate prob */
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else
	{
		/* calculate probability */
		for ( i = 0 ; i < instance.n_near ; i++ )
		{
			prob_ptr[i]=prob_ptr[i]/sum_prob;	
		}

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum <= rnd )//
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}
void move_to_next_normalized( ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 

    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0.0, sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];/* prepare to calculate prob */
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else
	{
		/* calculate probability */
		for ( i = 0 ; i < instance.n_near ; i++ )
		{
			prob_ptr[i]=prob_ptr[i]/sum_prob;
		}

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum <= rnd )//
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}

void map_ACS_move_to_next( long int i_map, ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 

    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0.0, sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

	if ( (q_0 > 0.0) && (ran01( &seed ) < q_0)  )
	{
		/* with a probability q_0 make the best possible choice
		   according to pheromone trails and heuristic information */
		/* we first check whether q_0 > 0.0, to avoid the very common case
		   of q_0 = 0.0 to have to compute a random number, which is
		   expensive computationally */
		map_neighbour_choose_best_next(i_map,a, phase);
		return;
    }

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = map[i_map].total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];/* prepare to calculate prob */
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	map_choose_best_next( i_map,a, phase );
    }     
    else
	{
		/* calculate probability */
		for ( i = 0 ; i < instance.n_near ; i++ )
		{
			prob_ptr[i]=prob_ptr[i]/sum_prob;
		}

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum <= rnd )//
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}



void ACS_move_to_next( ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 

    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0.0, sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

	if ( (q_0 > 0.0) && (ran01( &seed ) < q_0)  )
	{
		/* with a probability q_0 make the best possible choice
		   according to pheromone trails and heuristic information */
		/* we first check whether q_0 > 0.0, to avoid the very common case
		   of q_0 = 0.0 to have to compute a random number, which is
		   expensive computationally */
		neighbour_choose_best_next(a, phase);
		return;
    }

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];/* prepare to calculate prob */
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else
	{
		/* calculate probability */
		for ( i = 0 ; i < instance.n_near ; i++ )
		{
			prob_ptr[i]=prob_ptr[i]/sum_prob;
		}

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum <= rnd )//
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}
void map_choose_best_next( long int i_map, ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
    double   value_best;

    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( (*a).visited[city] ) 
	    ; /* city already visited, do nothing */
	else {
	    if ( map[i_map].total[current_city][city] > value_best ) {
		next_city = city;
		value_best = map[i_map].total[current_city][city];
	    }
	} 
    }
    (*a).tour[phase] = next_city;
    (*a).visited[next_city] = TRUE;
}
void choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
    double   value_best;

    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( (*a).visited[city] ) 
	    ; /* city already visited, do nothing */
	else {
	    if ( total[current_city][city] > value_best ) {
		next_city = city;
		value_best = total[current_city][city];
	    }
	} 
    }
    (*a).tour[phase] = next_city;
    (*a).visited[next_city] = TRUE;
}

void map_neighbour_choose_best_next( long int i_map, ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int i, current_city, next_city, help_city;
    double   value_best, help;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    DEBUG ( assert ( 0 <= current_city && current_city < n ); )
    value_best = -1.;             /* values in total matix are always >= 0.0 */    
	for ( i = 0 ; i < instance.n_near; i++ ) {
	help_city = instance.nn_list[current_city][i];
	if ( (*a).visited[help_city] ) 
	    ;   /* city already visited, do nothing */
	else {
	    help = map[i_map].total[current_city][help_city];
	    if ( help > value_best ) {
		value_best = help;
		next_city = help_city;
	    }
	}
    }
    if ( next_city == n )
	/* all cities in nearest neighbor list were already visited */
	choose_best_next( a, phase );
    else {
	(*a).tour[phase] = next_city;
	(*a).visited[next_city] = TRUE;
    }
}

//some problem happened here
//the space for the ant became strange

void neighbour_choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int i, current_city, next_city, help_city;
    double   value_best, help;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    DEBUG ( assert ( 0 <= current_city && current_city < n ); )
    value_best = -1.;             /* values in total matix are always >= 0.0 */    
	for ( i = 0 ; i < instance.n_near; i++ ) {
	help_city = instance.nn_list[current_city][i];
	if ( (*a).visited[help_city] ) 
	    ;   /* city already visited, do nothing */
	else {
	    help = total[current_city][help_city];
	    if ( help > value_best ) {
		value_best = help;
		next_city = help_city;
	    }
	}
    }
    if ( next_city == n )
	/* all cities in nearest neighbor list were already visited */
	choose_best_next( a, phase );
    else {
	DEBUG( assert ( 0 <= next_city && next_city < n); )
	DEBUG( assert ( value_best > 0.0 ); )
	DEBUG( assert ( (*a).visited[next_city] == FALSE ); )
	(*a).tour[phase] = next_city;
	(*a).visited[next_city] = TRUE;
    }
}
void map_update_statistics(long int i_map, long int n_try )
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          
      OUTPUT:         
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    long int iteration_best_ant;

    iteration_best_ant = map_find_best(i_map); /* iteration_best_ant is a global variable */	
	map[i_map].elapsed = map_elapsed_time(i_map);
	check_convergency(i_map);
    if ( map[i_map].ant[iteration_best_ant].tour_length < map[i_map].best_so_far_ant->tour_length )
	{
		map[i_map].convergance_num=1;
		copy_from_to( &map[i_map].ant[iteration_best_ant], map[i_map].best_so_far_ant );
		map[i_map].iter_to_best = map[i_map].iteration;
		printf("try %ld: Map %ld, %lf in iteration\t%ld\tafter seconds\t%lf\n",
			   n_try,i_map,map[i_map].ant[iteration_best_ant].tour_length, map[i_map].iter_to_best, map[i_map].elapsed);
		fprintf(report,"try %ld: Map %ld, %lf in iteration\t%ld\tafter seconds\t%lf\n",
			   n_try,i_map,map[i_map].ant[iteration_best_ant].tour_length, map[i_map].iter_to_best, map[i_map].elapsed);

	}
	else
		map[i_map].convergance_num++;
	if(map[i_map].iteration%100==0)
	{
		printf("***** try %ld: Map %ld, %lf  iteration %ld  atfer seconds\t%lf *****\n",
		       n_try,i_map,map[i_map].best_so_far_ant->tour_length,map[i_map].iteration, map[i_map].elapsed);
	}
}

void update_statistics( long int n_try )
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    long int iteration_best_ant;

    iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */

	elapsed = elapsed_time(); /* best sol found after time_used */
    if ( ant[iteration_best_ant].tour_length < (*best_so_far_ant).tour_length )
	{

		copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
		iter_to_best = iteration;

		printf("%lf in iteration\t%ld\t\tafter seconds\t%lf\n",ant[iteration_best_ant].tour_length, iter_to_best, elapsed);

	}
	if(iteration%100==0)
		printf("***** %lf  iteration %ld  atfer seconds\t%lf *****\n",(*best_so_far_ant).tour_length,iteration, elapsed);
}



void copy_from_to(ant_struct *a1, ant_struct *a2) 
{
/*    
      FUNCTION:       copy solution from ant a1 into ant a2
      INPUT:          pointers to the two ants a1 and a2 
      OUTPUT:         none
      (SIDE)EFFECTS:  a2 is copy of a1
*/
    long int   i;
  
    (*a2).tour_length = (*a1).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	(*a2).tour[i] = (*a1).tour[i];
    }
    (*a2).tour[n] = (*a2).tour[0];
}


void pheromone_trail_update( void )  
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
                      according to the rules defined by the various ACO algorithms.
*/
{
	if(ls_flag==1)
		evaporation_nn_list();
	else
	    evaporation();

	ras_update(); 
//	as_update(); 

    compute_total_information();
}
void map_ACS_pheromone_trail_update(long int i_map)
{
/*
	if(evapor_nn_flag==1)
		map_evaporation_nn_list(i_map);
	else
*/
		map_evaporation(i_map);

	map_global_acs_pheromone_update( i_map, map[i_map].best_so_far_ant );
}

void ACS_pheromone_trail_update(void)
{

	if(evapor_nn_flag==1)
		evaporation_nn_list();
	else
		evaporation();

	global_acs_pheromone_update( best_so_far_ant );
}

void EACS_pheromone_trail_update(void)
{
	if(evapor_nn_flag==1)
		evaporation_nn_list();
	else
		evaporation();

	global_acs_pheromone_update( best_so_far_ant );
}

void EACO_pheromone_trail_update( void )  
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
                      according to the rules defined by the various ACO algorithms.
*/
{
/*
	if(evapor_nn_flag==1)
		evaporation_nn_list();
	else
	    evaporation();
*/
		EACO_ras_update(); 

	    compute_total_information();
}

void map_evaporation( long int i_map )
/*    
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{ 
    long int    i, j;

    TRACE ( printf("pheromone evaporation\n"); );

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j <= i ; j++ ) {
	    map[i_map].pheromone[i][j] = (1 - rho) * map[i_map].pheromone[i][j];
	    map[i_map].pheromone[j][i] = map[i_map].pheromone[i][j];
	}
    }
}
void evaporation( void )
/*    
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{ 
    long int    i, j;

    TRACE ( printf("pheromone evaporation\n"); );

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j <= i ; j++ ) {
	    pheromone[i][j] = (1 - rho) * pheromone[i][j];
	    pheromone[j][i] = pheromone[i][j];
	}
    }
}

void map_evaporation_nn_list( long int i_map )
/*    
      FUNCTION:      simulation of the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure 
                     only considers links between a city and those cities
		     of its candidate list
*/
{ 
    long int    i, j, help_city;

    TRACE ( printf("pheromone evaporation nn_list\n"); );

    for ( i = 0 ; i < n ; i++ )
	{
		for ( j = 0 ; j < instance.n_near; j++ )
		{
		    help_city = instance.nn_list[i][j];
		    map[i_map].pheromone[i][help_city] = (1 - rho) * map[i_map].pheromone[i][help_city];
			map[i_map].pheromone[help_city][i] = map[i_map].pheromone[i][help_city];
		}
    }
}
void evaporation_nn_list( void )
/*    
      FUNCTION:      simulation of the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure 
                     only considers links between a city and those cities
		     of its candidate list
*/
{ 
    long int    i, j, help_city;

    TRACE ( printf("pheromone evaporation nn_list\n"); );

    for ( i = 0 ; i < n ; i++ ) {
		for ( j = 0 ; j < instance.n_near; j++ ) {
	    help_city = instance.nn_list[i][j];
	    pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
//		pheromone[help_city][i] = pheromone[i][help_city];
		}
    }
}


void as_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
}

void ras_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    long int i, k, target;
	double b;
    double *help_b;

    help_b = malloc( n_ants  * sizeof(double) );
    for ( k = 0 ; k < n_ants ; k++ )
		help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ )
	{
		b = help_b[0]; target = 0;
		for ( k = 0 ; k < n_ants ; k++ )
		{
			if ( help_b[k] < b ) 
			{
			b = help_b[k]; target = k;
			}
		}
		help_b[target] = LONG_MAX;
		global_update_pheromone_weighted( &ant[target], (ras_ranks-i-1) );
    }
    global_update_pheromone_weighted( best_so_far_ant, (ras_ranks) );	
    free ( help_b );
}


void EACO_ras_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    long int i, k, target;
	double b;
    double *help_b;

    help_b = malloc( n_ants  * sizeof(double) );
    for ( k = 0 ; k < n_ants ; k++ )
		help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ )
	{
		b = help_b[0]; target = 0;
		for ( k = 0 ; k < n_ants ; k++ )
		{
			if ( help_b[k] < b ) 
			{
			b = help_b[k]; target = k;
			}
		}
		help_b[target] = LONG_MAX;
		global_update_pheromone_weighted( &ant[target], epsilon*(ras_ranks-i-1) );/*　epsilon関係　epsilon*（rank-i-1）　*/
    }
    global_update_pheromone_weighted( best_so_far_ant, epsilon*(ras_ranks) );	  /*　epsilon関係　epsilon*（rank-i-1）　*/
    free ( help_b );
}


void global_update_pheromone( ant_struct *a )
/*    
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

    d_tau = Q / (double) (*a).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];
	pheromone[j][h] += d_tau;
	pheromone[h][j] = pheromone[j][h];
    }
}

void global_update_pheromone_weighted( ant_struct *a, double weight)
/*    
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

    d_tau = ( Q * weight ) / (double) (*a).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];
	pheromone[j][h] += d_tau;
	pheromone[h][j] = pheromone[j][h];
    }
}
void map_global_acs_pheromone_update( long int i_map, ant_struct *a )
/*    
      FUNCTION:      reinforces the edges used in ant's solution as in ACS
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    d_tau = 1.0 / (double) (*a).tour_length;

    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];

	map[i_map].pheromone[j][h] = (1. - rho) * map[i_map].pheromone[j][h] + rho * d_tau;
	map[i_map].pheromone[h][j] = map[i_map].pheromone[j][h];

	map[i_map].total[h][j] = pow(map[i_map].pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
	map[i_map].total[j][h] = map[i_map].total[h][j];
    }
}

void global_acs_pheromone_update( ant_struct *a )
/*    
      FUNCTION:      reinforces the edges used in ant's solution as in ACS
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("acs specific: global pheromone update\n"); );

    d_tau = 1.0 / (double) (*a).tour_length;

    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];

	pheromone[j][h] = (1. - rho) * pheromone[j][h] + rho * d_tau;
	pheromone[h][j] = pheromone[j][h];

	total[h][j] = pow(pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
	total[j][h] = total[h][j];
    }
}
void map_report( long int i_map, long int ntry ) 
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
*/
{
	long int i,j;
//	 checkTour( map[i_map].best_so_far_ant->tour );
	 map[i_map].best_in_try[ntry] = (*best_so_far_ant).tour_length;
	 map[i_map].time_best_found[ntry] = elapsed;
	 map[i_map].best_found_at[ntry] = iter_to_best;
     map[i_map].time_total_run[ntry] = elapsed_time();

	 total_elapsed= map_elapsed_time(0); 

	 if(map[i_map].best_so_far_ant->tour_length < (*best_so_far_ant).tour_length)
	 {
		 copy_from_tour_to_tour(map[i_map].best_so_far_ant->tour,(*best_so_far_ant).tour);
		 (*best_so_far_ant).tour_length=map[i_map].best_so_far_ant->tour_length;
	 }

	printf("IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t Time %lf\n\n",
		ntry,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
	fprintf(report,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n\n",
		 ntry,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);
	for(i=0;i<=n;i++)
	{
		fprintf(report,"%d\t%f\t%f\n",map[i_map].best_so_far_ant->tour[i],cinstance.nodeptr[map[i_map].best_so_far_ant->tour[i]].x,cinstance.nodeptr[map[i_map].best_so_far_ant->tour[i]].y);
	}
	fprintf(best_report,"IN TRY %ld:Map %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\n",
		 ntry,i_map, map[i_map].best_so_far_ant->tour_length, map[i_map].iter_to_best,total_elapsed);

	fflush(report); 
	fflush(best_report); 

	for(i=0; i<n; i++)
			for(j=0; j<n; j++)
			{
//				fprintf(best_report,"city %ld %ld ph: %0.15f\n",i,j,map[i_map].pheromone[i][j]);
			}


 }
void exit_try( long int ntry ) 
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
*/
{
  checkTour( (*best_so_far_ant).tour );
/*    printTourFile( (*best_so_far_ant).tour ); */

  best_in_try[ntry] = (*best_so_far_ant).tour_length;
  time_best_found[ntry] = elapsed;
  best_found_at[ntry] = iter_to_best;

  time_total_run[ntry] = elapsed_time();
  printf("In try %ld\t Best: %lf\t Iterations: %ld\t Time %lf\t Tot.time %lf\n",
	  ntry, best_in_try[ntry], best_found_at[ntry], time_best_found[ntry],time_total_run[ntry]);
  fprintf(report,"In try %ld\t Best: %lf\t Iterations: %ld\t\t Time %lf\t Tot.time %lf\n",
	  ntry, best_in_try[ntry], best_found_at[ntry], time_best_found[ntry],time_total_run[ntry]);
 
//  fprintf(report,"end try %ld\n\n",ntry);
 
  TRACE (output_solution();)
  fflush(report); 
 }


void exit_program( void ) 
/*    
      FUNCTION:       save some final statistical information on a trial once it finishes
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       
*/
{
  double best_tour_length, worst_tour_length;
  double   t_avgbest, t_stdbest, t_avgtotal, t_stdtotal;
  double   avg_sol_quality = 0.0, avg_iter_to_bst = 0.0, stddev_best, stddev_iterations;
  double opt_rate=0.0;
  long int i=0;
  long int opt_num=0;

  for(; i<max_tries; i++)
  {
	if(best_in_try[i]==optimal)
		opt_num++;
  }
  opt_rate=(double)opt_num/max_tries;

  best_tour_length = best_of_vector( best_in_try ,max_tries );
  worst_tour_length = worst_of_vector( best_in_try , max_tries );/**/

  avg_iter_to_bst = mean( best_found_at , max_tries );
  stddev_iterations = std_deviation( best_found_at, max_tries, avg_iter_to_bst );

  avg_sol_quality = meanr( best_in_try , max_tries );
  stddev_best = std_deviationr( best_in_try, max_tries, avg_sol_quality);

  t_avgbest = meanr( time_best_found, max_tries );
  t_stdbest = std_deviationr( time_best_found, max_tries, t_avgbest);

  t_avgtotal = meanr( time_total_run, max_tries );
  t_stdtotal = std_deviationr( time_total_run, max_tries, t_avgtotal);
 
  printf(" optimal rate = %f\n", opt_rate );
  printf(" t_avgbest = %f\n", t_avgbest );
  printf(" t_avgtotal = %f\n", t_avgtotal );

  fprintf(report,"optimal rate = %f\n", opt_rate );
  fprintf(report,"\nAverage-Best: %.2f\t Stddev-Best: %.2f", avg_sol_quality, stddev_best);
  fprintf(report,"\nAverage-Iterations: %.2f \t Stddev Iterations: %.2f", avg_iter_to_bst, stddev_iterations);
  fprintf(report,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest);  
  fprintf(report,"\nAvg.time-total: %.2f stddev.time-total: %.2f\n", t_avgtotal, t_stdtotal); 

  fflush(report);
}


void copy_ant_to_tour(ant_struct *ant, long int *tour)
{
    int   i;
    for ( i = 0 ; i < n ; i++ ) {
	tour[i] = (*ant).tour[i];
	}
	tour[n] = (*ant).tour[0];
}

void copy_from_tour_to_tour(long int *src, long int *des)
{
    int   i;
    for ( i = 0 ; i < n ; i++ ) {
	des[i] = src[i];
	}
	des[n] = src[0];
}
void copy_tour_to_ant(long int *tour, ant_struct *ant)
{
    int   i;
    for ( i = 0 ; i < n ; i++ ) {
	(*ant).tour[i] = tour[i];
	}
	(*ant).tour[n] = tour[0];
}

void add_tmp_phe( ant_struct *a )

{
    long int i, h, j;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

    d_tau = (1-epsilon) / (double) (*a).tour_length;	/* 元は 1.0 / (double) (*a).tour_length */
    for ( i = 0 ; i < n ; i++ ) {
	h = (*a).tour[i];
	j = (*a).tour[i+1];
	pheromone[h][j] += d_tau;
	pheromone[j][h] = pheromone[j][h];

	total[h][j] = pow(pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
    total[j][h] = total[h][j];
    }
}

void minus_tmp_phe( ant_struct *a, double tmp_length )
{
    long int i, h, j;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

	d_tau = (1-epsilon) / tmp_length;				/* 1.0 / tmp_length */
    for ( i = 0 ; i < n ; i++ ) {
	h = (*a).tour[i];
	j = (*a).tour[i+1];
	pheromone[h][j] -= d_tau;
	pheromone[j][h] = pheromone[j][h];

	total[h][j] = pow(pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
    total[j][h] = total[h][j];
    }
}


void minus_tmp_phe_tour( long int *tour, double tmp_length )
{
    long int i, h, j;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

	d_tau = (1-epsilon) / tmp_length;				/* 1.0 / tmp_length */
    for ( i = 0 ; i < n ; i++ ) {
	h = tour[i];
	j = tour[i+1];
	pheromone[h][j] -= d_tau;
	pheromone[j][h] = pheromone[j][h];

	total[h][j] = pow(pheromone[h][j], alpha) * pow(heuristic[h][j],beta);
    total[j][h] = total[h][j];
    }
}
void map_2opt_local_search(long int i_map)
{
    long int k;
	for ( k = 0 ; k < n_ants ; k++ )
	{
		two_opt_first( map[i_map].ant[k].tour );
		map[i_map].ant[k].tour_length = compute_tour_length( map[i_map].ant[k].tour );
	}
}
void map_local_search(long int i_map)
{
    long int k;
	for ( k = 0 ; k < n_ants ; k++ )
	{
		three_opt_first( map[i_map].ant[k].tour );
		map[i_map].ant[k].tour_length = compute_tour_length( map[i_map].ant[k].tour );
	}
}



void local_search(void)
{
    long int k;
	for ( k = 0 ; k < n_ants ; k++ )
	{
		three_opt_first( ant[k].tour );
		ant[k].tour_length = compute_tour_length( ant[k].tour );
	}
}

double individual_local_search(long int* tour )
{
		three_opt_first(tour );
		return compute_tour_length( tour );
}


void two_opt_first( long int *tour ) 
/*    
      FUNCTION:       2-opt a tour 
      INPUT:          pointer to the tour that undergoes local optimization
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 2-opt
      COMMENTS:       the neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers online available from David S. Johnson.
*/
{
    long int c1, c2;             /* cities considered for an exchange */
    long int s_c1, s_c2;         /* successor cities of c1 and c2     */
    long int p_c1, p_c2;         /* predecessor cities of c1 and c2   */   
    long int pos_c1, pos_c2;     /* positions of cities c1, c2        */
    long int i, j, h, l;
    long int improvement_flag, improve_node, help, n_improves = 0, n_exchanges=0;
    long int h1=0, h2=0, h3=0, h4=0;
    double radius;             /* radius of nn-search */
    double gain = 0;
    long int *random_vector;
    long int *pos;               /* positions of cities in tour */ 
    long int *dlb;               /* vector containing don't look bits */ 
  
    pos = malloc(n * sizeof(long int));
    dlb = malloc(n * sizeof(long int));
    for ( i = 0 ; i < n ; i++ ) {
	pos[tour[i]] = i;
	dlb[i] = FALSE;
    }

    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n );

    while ( improvement_flag ) {

	improvement_flag = FALSE;

	for (l = 0 ; l < n; l++) {

	    c1 = random_vector[l]; 
	    DEBUG ( assert ( c1 < n && c1 >= 0); )
		if ( dlb_flag && dlb[c1] )
		    continue;
	    improve_node = FALSE;
	    pos_c1 = pos[c1];
	    s_c1 = tour[pos_c1+1];
	    radius = instance.distance[c1][s_c1];

	    /* First search for c1's nearest neighbours, use successor of c1 */
	    for ( h = 0 ; h < instance.n_near ; h++ ) {
		c2 = instance.nn_list[c1][h]; /* exchange partner, determine its position */
		if ( radius > instance.distance[c1][c2] ) {
		    s_c2 = tour[pos[c2]+1];
		    gain =  - radius + instance.distance[c1][c2] + 
			instance.distance[s_c1][s_c2] - instance.distance[c2][s_c2];
		    if ( gain < 0 ) {
			h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
			improve_node = TRUE;
			goto exchange2opt;
		    }
		}
		else 
		    break;
	    }      
	    /* Search one for next c1's h-nearest neighbours, use predecessor c1 */
	    if (pos_c1 > 0)
		p_c1 = tour[pos_c1-1];
	    else 
		p_c1 = tour[n-1];
	    radius = instance.distance[p_c1][c1];
	    for ( h = 0 ; h < instance.n_near ; h++ ) {
		c2 = instance.nn_list[c1][h];  /* exchange partner, determine its position */
		if ( radius > instance.distance[c1][c2] ) {
		    pos_c2 = pos[c2];
		    if (pos_c2 > 0)
			p_c2 = tour[pos_c2-1];
		    else 
			p_c2 = tour[n-1];
		    if ( p_c2 == c1 )
			continue;
		    if ( p_c1 == c2 )
			continue;
		    gain =  - radius + instance.distance[c1][c2] + 
			instance.distance[p_c1][p_c2] - instance.distance[p_c2][c2];
		    if ( gain < 0 ) {
			h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2; 
			improve_node = TRUE;
			goto exchange2opt;
		    }
		}
		else 
		    break;
	    }      
	    if (improve_node) {
	    exchange2opt:
		n_exchanges++;
		improvement_flag = TRUE;
		dlb[h1] = FALSE; dlb[h2] = FALSE;
		dlb[h3] = FALSE; dlb[h4] = FALSE;
		/* Now perform move */
		if ( pos[h3] < pos[h1] ) {
		    help = h1; h1 = h3; h3 = help;
		    help = h2; h2 = h4; h4 = help;
		}
		if ( pos[h3] - pos[h2] < n / 2 + 1) {
		    /* reverse inner part from pos[h2] to pos[h3] */
		    i = pos[h2]; j = pos[h3];
		    while (i < j) {
			c1 = tour[i];
			c2 = tour[j];
			tour[i] = c2;
			tour[j] = c1;
			pos[c1] = j;
			pos[c2] = i;
			i++; j--;
		    }
		}
		else {
		    /* reverse outer part from pos[h4] to pos[h1] */
		    i = pos[h1]; j = pos[h4];
		    if ( j > i )
			help = n - (j - i) + 1;
		    else 
			help = (i - j) + 1;
		    help = help / 2;
		    for ( h = 0 ; h < help ; h++ ) {
			c1 = tour[i];
			c2 = tour[j];
			tour[i] = c2;
			tour[j] = c1;
			pos[c1] = j;
			pos[c2] = i;
			i--; j++;
			if ( i < 0 )
			    i = n-1;
			if ( j >= n )
			    j = 0;
		    }
		    tour[n] = tour[0];
		}
	    } else {
		dlb[c1] = TRUE;
	    }
      
	}
	if ( improvement_flag ) {
	    n_improves++;
	}
    }
    free( random_vector );
    free( dlb );
    free( pos );
}

void three_opt_first( long int *tour )

/*    
      FUNCTION:       3-opt the tour
      INPUT:          pointer to the tour that is to optimize
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 3-opt
      COMMENT:        this is certainly not the best possible implementation of a 3-opt 
                      local search algorithm. In addition, it is very lengthy; the main 
		      reason herefore is that awkward way of making an exchange, where 
		      it is tried to copy only the shortest possible part of a tour. 
		      Whoever improves the code regarding speed or solution quality, please 
		      drop me the code at stuetzle no@spam informatik.tu-darmstadt.de
		      The neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers online available from David S. Johnson.
*/
{
    /* In case a 2-opt move should be performed, we only need to store opt2_move = TRUE,
       as h1, .. h4 are used in such a way that they store the indices of the correct move */

    long int   c1, c2, c3;           /* cities considered for an exchange */
    long int   s_c1, s_c2, s_c3;     /* successors of these cities        */
    long int   p_c1, p_c2, p_c3;     /* predecessors of these cities      */   
    long int   pos_c1, pos_c2, pos_c3;     /* positions of cities c1, c2, c3    */
    long int   i, j, h, g, l;
    long int   improvement_flag, help;
    long int   h1=0, h2=0, h3=0, h4=0, h5=0, h6=0; /* memorize cities involved in a move */
    double   diffs, diffp;
    long int   between = FALSE; 
    long int   opt2_flag;  /* = TRUE: perform 2-opt move, otherwise none or 3-opt move */
    long int   move_flag;  /* 
			      move_flag = 0 --> no 3-opt move 
			      move_flag = 1 --> between_move (c3 between c1 and c2)
			      move_flag = 2 --> not_between with successors of c2 and c3
			      move_flag = 3 --> not_between with predecessors of c2 and c3
			      move_flag = 4 --> cyclic move 
			   */
    double gain, move_value, radius, add1, add2;

    double decrease_breaks;    /* Stores decrease by breaking two edges (a,b) (c,d) */
    long int val[3];
    long int n1, n2, n3;
    long int *pos;               /* positions of cities in tour */ 
    long int *dlb;               /* vector containing don't look bits */ 
    long int *h_tour;            /* help vector for performing exchange move */ 
    long int *hh_tour;           /* help vector for performing exchange move */ 
    long int *random_vector;

    pos = malloc(n * sizeof(long int));
    dlb = malloc(n * sizeof(long int));
    h_tour = malloc(n * sizeof(long int));
    hh_tour = malloc(n * sizeof(long int));

    for ( i = 0 ; i < n ; i++ ) {
	pos[tour[i]] = i;
	dlb[i] = FALSE;
    }
    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n );

    while ( improvement_flag ) {
	move_value = 0;
	improvement_flag = FALSE;

	for ( l = 0 ; l < n ; l++ ) {

	    c1 = random_vector[l];
	    if ( dlb_flag && dlb[c1] )
		continue;
	    opt2_flag = FALSE;

	    move_flag = 0;
	    pos_c1 = pos[c1];
	    s_c1 = tour[pos_c1+1];
	    if (pos_c1 > 0)
		p_c1 = tour[pos_c1-1];
	    else 
		p_c1 = tour[n-1];

	    h = 0;    /* Search for one of the h-nearest neighbours */
	    while ( h < instance.n_near ) {

		c2   = instance.nn_list[c1][h];  /* second city, determine its position */
		pos_c2 = pos[c2];
		s_c2 = tour[pos_c2+1];
		if (pos_c2 > 0)
		    p_c2 = tour[pos_c2-1];
		else 
		    p_c2 = tour[n-1];
	  
		diffs = 0; diffp = 0;

		radius = instance.distance[c1][s_c1];
		add1   = instance.distance[c1][c2];

		/* Here a fixed radius neighbour search is performed */
		if ( radius > add1 ) {
		    decrease_breaks = - radius - instance.distance[c2][s_c2];
		    diffs =  decrease_breaks + add1 + instance.distance[s_c1][s_c2];
		    diffp =  - radius - instance.distance[c2][p_c2] + 
			instance.distance[c1][p_c2] + instance.distance[s_c1][c2];
		}
		else 
		    break;
		if ( p_c2 == c1 )  /* in case p_c2 == c1 no exchange is possible */
		    diffp = 0;
		if ( (diffs < move_value) || (diffp < move_value) ) {
		    improvement_flag = TRUE; 
		    if (diffs <= diffp) { 
			h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
			move_value = diffs; 
			opt2_flag = TRUE; move_flag = 0;
			/*     	    goto exchange; */
		    } else {
			h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; 
			move_value = diffp;  
			opt2_flag = TRUE; move_flag = 0;
			/*     	    goto exchange; */
		    }
		}
		/* Now perform the innermost search */
		g = 0;
		while (g < instance.n_near) {
	  
		    c3   = instance.nn_list[s_c1][g];
		    pos_c3 = pos[c3];
		    s_c3 = tour[pos_c3+1];
		    if (pos_c3 > 0)
			p_c3 = tour[pos_c3-1];
		    else 
			p_c3 = tour[n-1];
		  
		    if ( c3 == c1 ) {
			g++;
			continue;
		    }
		    else {
			add2 = instance.distance[s_c1][c3];
			/* Perform fixed radius neighbour search for innermost search */
			if ( decrease_breaks + add1 < add2 ) {
			  
			    if ( pos_c2 > pos_c1 ) {
				if ( pos_c3 <= pos_c2 && pos_c3 > pos_c1 )
				    between = TRUE;
				else 
				    between = FALSE;
			    }
			    else if ( pos_c2 < pos_c1 )
				if ( pos_c3 > pos_c1 || pos_c3 < pos_c2 )
				    between = TRUE;
				else 
				    between = FALSE;
			    else {
				printf(" Strange !!, pos_1 %ld == pos_2 %ld, \n",pos_c1,pos_c2);
			    }
			  
			    if ( between ) {
				/* We have to add edges (c1,c2), (c3,s_c1), (p_c3,s_c2) to get 
				   valid tour; it's the only possibility */
			      
				gain = decrease_breaks - instance.distance[c3][p_c3] +
				    add1 + add2 +
				    instance.distance[p_c3][s_c2];
			      
				/* check for improvement by move */
				if ( gain < move_value ) {
				    improvement_flag = TRUE; /* g = neigh_ls + 1; */
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 1;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = p_c3; h6 = c3;
				    goto exchange;
				} 
			    }
			    else {   /* not between(pos_c1,pos_c2,pos_c3) */
			      
				/* We have to add edges (c1,c2), (s_c1,c3), (s_c2,s_c3) */
			      
				gain = decrease_breaks - instance.distance[c3][s_c3] +
				    add1 + add2 + 
				    instance.distance[s_c2][s_c3];
			      
				if ( pos_c2 == pos_c3 ) {
				    gain = 20000;
				}
			      
				/* check for improvement by move */
				if ( gain < move_value ) {
				    improvement_flag = TRUE; /* g = neigh_ls + 1; */
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 2;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = c3; h6 = s_c3;
				    goto exchange;
				}
			      
				/* or add edges (c1,c2), (s_c1,c3), (p_c2,p_c3) */
				gain = - radius - instance.distance[p_c2][c2] 
				    - instance.distance[p_c3][c3] +
				    add1 + add2 + 
				    instance.distance[p_c2][p_c3];
			      
				if ( c3 == c2 || c2 == c1 || c1 == c3 || p_c2 == c1 ) {
				    gain = 2000000;
				}
			      
				if ( gain < move_value ) {
				    improvement_flag = TRUE;
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 3;
				    h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = p_c3; h6 = c3;
				    goto exchange;
				}
			      
				/* Or perform the 3-opt move where no subtour inversion is necessary 
				   i.e. delete edges (c1,s_c1), (c2,p_c2), (c3,s_c3) and 
				   add edges (c1,c2), (c3,s_c1), (p_c2,s_c3) */
			      
				gain = - radius - instance.distance[p_c2][c2] - 
				    instance.distance[c3][s_c3]
				    + add1 + add2 + instance.distance[p_c2][s_c3];
			      
				/* check for improvement */
				if ( gain < move_value ) {
				    improvement_flag = TRUE;
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 4;
				    improvement_flag = TRUE;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = c3; h6 = s_c3; 
				    goto exchange;
				}
			    }
			}
			else
			    g = instance.n_near + 1;
		    }
		    g++;
		}
		h++;
	    }
	    if ( move_flag || opt2_flag ) {
	    exchange:
		move_value = 0;

		/* Now make the exchange */
		if ( move_flag ) {
		    dlb[h1] = FALSE; dlb[h2] = FALSE; dlb[h3] = FALSE; 
		    dlb[h4] = FALSE; dlb[h5] = FALSE; dlb[h6] = FALSE;
		    pos_c1 = pos[h1]; pos_c2 = pos[h3]; pos_c3 = pos[h5];
		  
		    if ( move_flag == 4 ) {

			if ( pos_c2 > pos_c1 ) 
			    n1 = pos_c2 - pos_c1;
			else
			    n1 = n - (pos_c1 - pos_c2);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c1 > pos_c3 ) 
			    n3 = pos_c1 - pos_c3;
			else
			    n3 = n - (pos_c3 - pos_c1);
		      
			/* n1: length h2 - h3, n2: length h4 - h5, n3: length h6 - h1 */
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Now order the partial tours */
			h = 0;
			help = LONG_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
		      
			/* order partial tours according length */
			if ( h == 0 ) {
			    /* copy part from pos[h4] to pos[h5]
			       direkt kopiert: Teil von pos[h6] to pos[h1], it
			       remains the part from pos[h2] to pos[h3] */
			    j = pos[h4];
			    h = pos[h5];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    /* First copy partial tour 3 in new position */
			    j = pos[h4];
			    i = pos[h6];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c1) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {
			  
			    /* copy part from pos[h6] to pos[h1]
			       direkt kopiert: Teil von pos[h2] to pos[h3], it
			       remains the part from pos[h4] to pos[h5] */
			    j = pos[h6];
			    h = pos[h1];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    /* First copy partial tour 3 in new position */
			    j = pos[h6];
			    i = pos[h2];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c2) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 2 ) {
			    /* copy part from pos[h2] to pos[h3]
			       direkt kopiert: Teil von pos[h4] to pos[h5], it
			       remains the part from pos[h6] to pos[h1] */
			    j = pos[h2];
			    h = pos[h3];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
	      
			    /* First copy partial tour 3 in new position */
			    j = pos[h2];
			    i = pos[h4];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c3) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 1 ) {
		      
			if ( pos_c3 < pos_c2 ) 
			    n1 = pos_c2 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c2);
			if ( pos_c3 > pos_c1 ) 
			    n2 = pos_c3 - pos_c1 + 1;
			else
			    n2 = n - (pos_c1 - pos_c3 + 1);
			if ( pos_c2 > pos_c1 ) 
			    n3 = n - (pos_c2 - pos_c1 + 1);
			else
			    n3 = pos_c1 - pos_c2 + 1;
		      
			/* n1: length h6 - h3, n2: length h5 - h2, n2: length h1 - h3 */
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Now order the partial tours */
			h = 0;
			help = LONG_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h5] to pos[h2]
			       (inverted) and from pos[h4] to pos[h1] (inverted)
			       it remains the part from pos[h6] to pos[h3] */
			    j = pos[h5];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h1];
			    h = pos[h4];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h4];
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if (j >= n)
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {
			  
			    /* copy part from h3 to h6 (wird inverted) erstellen : */
			    j = pos[h3];
			    h = pos[h6];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j  < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h6];
			    i = pos[h4];
			  
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c1) {
				i++;
				j++;
				if ( j >= n)
				    j = 0;
				if ( i >= n)
				    i = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    i = 0;
			    tour[j] = h_tour[i];
			    pos[h_tour[i]] = j; 
			    while ( j != pos_c1 ) {
				j++;
				if ( j >= n )
				    j = 0;
				i++;
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h2] to pos[h5] and
			       from pos[h3] to pos[h6] (inverted), it
			       remains the part from pos[h4] to pos[h1] */
			    j = pos[h2];
			    h = pos[h5];
			    i = 0;
			    h_tour[i] =  tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			    j = pos_c2;
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h2];
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n)
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 2 ) {
		      
			if ( pos_c3 < pos_c1 ) 
			    n1 = pos_c1 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c1);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c2 > pos_c1 ) 
			    n3 = pos_c2 - pos_c1;
			else
			    n3 = n - (pos_c1 - pos_c2);
		      
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Determine which is the longest part */
			h = 0;
			help = LONG_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h3] to pos[h2]
			       (inverted) and from pos[h5] to pos[h4], it
			       remains the part from pos[h6] to pos[h1] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h5];
			    h = pos[h4];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h2];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
	      
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			    /*  	      getchar(); */
			}
			else if ( h == 1 ) {
			  
			    /* copy part from pos[h2] to pos[h3] and
			       from pos[h1] to pos[h6] (inverted), it
			       remains the part from pos[h4] to pos[h5] */
			    j = pos[h2];
			    h = pos[h3];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n  )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j =  n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			    j = pos[h6];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h1] to pos[h6]
			       (inverted) and from pos[h4] to pos[h5],
			       it remains the part from pos[h2] to
			       pos[h3] */
			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h4];
			    h = pos[h5];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n  )
				    j = 0;
				hh_tour[i] = tour[j];
				n2++;
			    }

			    j = pos[h4];
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 3 ) {
		      
			if ( pos_c3 < pos_c1 ) 
			    n1 = pos_c1 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c1);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c2 > pos_c1 ) 
			    n3 = pos_c2 - pos_c1;
			else
			    n3 = n - (pos_c1 - pos_c2);
			/* n1: length h6 - h1, n2: length h4 - h5, n2: length h2 - h3 */
		      
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Determine which is the longest part */
			h = 0;
			help = LONG_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h2] to pos[h3]
			       (inverted) and from pos[h4] to pos[h5]
			       it remains the part from pos[h6] to pos[h1] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h2];
			    h = pos[h5];
			    i = pos[h4];
			    tour[j] = h4;
			    pos[h4] = j;
			    while ( i != h ) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j;
			    }
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i < n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {

			    /* copy part from pos[h3] to pos[h2]
			       (inverted) and from  pos[h6] to pos[h1],
			       it remains the part from pos[h4] to pos[h5] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0  )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h6];
			    h = pos[h1];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n )
				    j = 0;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h6];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }

			    for ( i = 0 ; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h4] to pos[h5]
			       (inverted) and from pos[h6] to pos[h1] (inverted)
			       it remains the part from pos[h2] to pos[h3] */
			    j = pos[h5];
			    h = pos[h4];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }

			    j = pos[h4];
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else {
			printf(" Some very strange error must have occurred !!!\n\n");
			exit(0);
		    }
		}
		if (opt2_flag) {

		    /* Now perform move */
		    dlb[h1] = FALSE; dlb[h2] = FALSE;
		    dlb[h3] = FALSE; dlb[h4] = FALSE;
		    if ( pos[h3] < pos[h1] ) {
			help = h1; h1 = h3; h3 = help;
			help = h2; h2 = h4; h4 = help;
		    }
		    if ( pos[h3]-pos[h2] < n / 2 + 1) {
			/* reverse inner part from pos[h2] to pos[h3] */
			i = pos[h2]; j = pos[h3];
			while (i < j) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i++; j--;
			}
		    }
		    else {
			/* reverse outer part from pos[h4] to pos[h1] */
			i = pos[h1]; j = pos[h4];
			if ( j > i )
			    help = n - (j - i) + 1;
			else 
			    help = (i - j) + 1;
			help = help / 2;
			for ( h = 0 ; h < help ; h++ ) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i--; j++;
			    if ( i < 0 )
				i = n - 1;
			    if ( j >= n )
				j = 0;
			}
			tour[n] = tour[0];
		    }
		}
	    }
	    else {
		dlb[c1] = TRUE;
	    }
	}
    }
    free( random_vector );
    free( h_tour );
    free( hh_tour );
    free( pos );
    free( dlb );
}


long int * generate_random_permutation( long int n )
/*    
      FUNCTION:       generate a random permutation of the integers 0 .. n-1
      INPUT:          length of the array
      OUTPUT:         pointer to the random permutation
      (SIDE)EFFECTS:  the array holding the random permutation is allocated in this 
                      function. Don't forget to free again the memory!
      COMMENTS:       only needed by the local search procedures
*/
{
   long int  i, help, node, tot_assigned = 0;
   double    rnd;
   long int  *r;

   r = malloc(n * sizeof(long int));  

   for ( i = 0 ; i < n; i++) 
     r[i] = i;

   for ( i = 0 ; i < n ; i++ ) {
     /* find (randomly) an index for a free unit */ 
     rnd  = ran01 ( &seed );
     node = (long int) (rnd  * (n - tot_assigned)); 
     assert( i + node < n );
     help = r[i];
     r[i] = r[i+node];
     r[i+node] = help;
     tot_assigned++;
   }
   return r;
}

long int check_convergency(long int i_map)
{
	double sum=0.0;
	long int i=0;
	for(i=0; i<n_ants; i++)
	{
		sum+=map[i_map].ant[i].tour_length;
	}

	if(sum==n_ants*(map[i_map].ant[0].tour_length) )
	{
//		printf("converge to %lf\n", map[i_map].ant[0].tour_length);
		return 1;
	}
	else return 0;

}

void swap_pheromone(double**a, double**b)
{
	double**tem=generate_double_matrix(n,n);
	long int i=0,j=0;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			tem[i][j]=a[i][j];
			a[i][j]=b[i][j];
			b[i][j]=tem[i][j];
		}
	}
	free(tem);
}
void average_pheromone(double**a, double**b)
{
	long int i=0,j=0;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			a[i][j]=(a[i][j]+b[i][j])/2;
			b[i][j]=a[i][j];
		}
	}
}

void island_crossover(long int island_num)
{
	double ***phe;
	long int pa1,pa2;
	long int is;
	long int i,j;
	phe=malloc(sizeof(double**)*island_num);
	for(is=0; is<island_num; is++)
	{
		pa1=rand()%island_num;
		pa2=rand()%island_num;
		phe[is]=generate_double_matrix( n, n );
		for(i=0; i<n; i++)
			for(j=0; j<n; j++)
			{
				phe[is][i][j]=(map[pa1].pheromone[i][j] + map[pa2].pheromone[i][j])/2;
			}
	}
	for(is=0; is<island_num; is++)
	{
		for(i=0; i<n; i++)
			for(j=0; j<n; j++)
			{
				map[is].pheromone[i][j]=phe[is][i][j];
			}
	}	
	for(i=0; i<island_num; i++)
	{
		free(phe[i]);
	}
	free(phe);

}

void island_crossover_c(long int island_num)
{
	double **phe;
	long int is;
	long int i,j;
	long int inew,jnew;

	printf("cross over the pheromone from each map\n");

	phe=generate_double_matrix( tn, tn );

	for(i=0; i<tn; i++)
		{
			for(j=0; j<tn; j++)
			{
				phe[i][j]=map[island_num].pheromone[i][j];
			}
		}


	for(is=0; is<island_num; is++)
	{
		n=map[is].mn;
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
				inew=relationship[is][i];
				jnew=relationship[is][j];
				phe[inew][jnew]=map[is].pheromone[i][j];
//				printf("inew:%d jnew:%d is:%d i:%d j:%d\n",inew,jnew,is,i,j);
//				printf("phe: %0.15f\n",phe[inew][jnew]);
			}
		}
	}

		for(i=0; i<tn; i++)
		{
			for(j=0; j<tn; j++)
			{
				map[island_num].pheromone[i][j]=phe[i][j];
			}
		}

	n=tn;



//	instance=cinstance;

	
	
	map_compute_total_information(island_num);
	
	free(phe);

}