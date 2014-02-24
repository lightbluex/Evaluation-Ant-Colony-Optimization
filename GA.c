#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"AS.h"
#include"GA.h"
#include"utilities.h"
#include"Map.h"

void map_GA(long int i_map)
{
	/*best N of offsprings substide worst N of parents*/
	long int gen_n=500;
	long int die_n= n_ants/15;/*set default*/

	long int** offsprings = generate_int_matrix(n_ants,n+1);/*put offsprings into*/
	double  *  offspring_length = malloc( n_ants*sizeof(double) );
	double  *  offspring_fitness = malloc( n_ants*sizeof(double) );
	long int*  offspring_rank   = malloc( n_ants*sizeof(long int) );
	double  *  parent_length    = malloc( n_ants*sizeof(double) );
	double  *  parent_fitness    = malloc( n_ants*sizeof(double) );

	long int*  parent_rank      = malloc( n_ants*sizeof(long int) );

	long int*  parent_a;
	long int*  parent_b;

	long int a,b;
	long int generation=0,j=0;
	double sum_len=0;
	/*copy parents,and check if converge*/
	for(j=0; j<n_ants; j++)
	{
		parent_length[j]=map[i_map].ant[j].tour_length;
		parent_fitness[j]=1/parent_fitness[j];
	}
	for(generation=0; generation<gen_n; generation++)
	{
		sum_len=0;/*if converged, break*/
		for(j=0; j<n_ants; j++)
		{
			a = rand()%n_ants;/*randomly get parents*/
			b = rand()%n_ants;
		//	a = roulette(parent_fitness,n_ants); 
		//	b = roulette(parent_fitness,n_ants); 
			parent_a=map[i_map].ant[a].tour;
			parent_b=map[i_map].ant[b].tour;
			crossover(parent_a,parent_b,offsprings[j]);
			checkTour(offsprings[j]);
			offspring_length[j]=compute_tour_length(offsprings[j]);/*=fitness*/
			offspring_fitness[j]=1/offspring_length[j];
			offspring_rank[j] =j;
//			parent_rank   [j] =j;
			sum_len += map[i_map].ant[j].tour_length;/*break if converged */
		}
		if(sum_len == n_ants*parent_length[0] )
		{
			printf("generation=%d\n",generation);
			break;
		}
		/*sort*/
		sort2( offspring_length, offspring_rank, 0, n_ants-1 );
//		sort2( parent_length   , parent_rank   , 0, n_ants-1 );

		/*natural select*/
		for(j=0; j< die_n;j++)
		{			
			/*copy tour*/
		  //a=parent_rank[n_ants-1-j];
			a=roulette(parent_length,n_ants);/*if length is large, then selected probability is large*/
	      //a = rand()%n_ants;/*randomly get parents*/
		  //b=offspring_rank[j];
			b=roulette(offspring_fitness,n_ants);
			copy_from_tour_to_tour(offsprings[b], map[i_map].ant[a].tour);
			map[i_map].ant[a].tour_length = offspring_length[j];/*renew parents=new generation*/
		//	map[i_map].ant[a].tour_length = offspring_length[b];
		//    parent_length[a]= offspring_length[b];
		    parent_length[a]= offspring_length[j];/*renew copy*/
			parent_fitness[a] = 1/offspring_length[j];
		}
		/*mutate by LocalSearch*/
		a = rand()%n_ants;
		//a=roulette(parent_length,n_ants);
		map[i_map].ant[a].tour_length = individual_local_search( map[i_map].ant[a].tour );
		parent_length[a]=map[i_map].ant[a].tour_length;
		parent_fitness[a]=1/parent_length[a];
	}	//end for
	free(offsprings);
	free(offspring_length);
	free(offspring_fitness);
	free(parent_fitness);
	free(offspring_rank);
	free(parent_length);
	free(parent_rank);
}


long int locate(long int *tour, long int e)
{
	long int i=0;
	for(; i<n; i++)
	{
		if(tour[i]==e)
			return i;
	}
	return -1;//if return -1, error!
}



void crossover(long int *parent_a, long int *parent_b, long int* offspring)
{
	long int *tmp_offspring   = calloc( 2*n-1, sizeof(long int) );/*set to 1, if city is selected*/
	long int *in_flag   = calloc( n, sizeof(long int) );/*set to 1, if city is selected*/
	long int * last;
	long int x=0;
	long int px_a=0;
	long int px_b=0;
	long int a_flag=1;
	long int b_flag=1;
	long int center=n-1;
	long int step = 0;
	long int side_a=0;
	long int side_b=0;
	long int pos=0;/*position of offspring*/
	long int i,j;
	x = rand()%n;//get a random city
	px_a = locate(parent_a,x);/*get postion of x from parent a and b*/
	px_b = locate(parent_b,x);

	tmp_offspring[ center ]=x;/*put x into middle bit*/
	in_flag[x]  =1;
	do
	{
		if(a_flag==1)
		{
			px_a = px_a-1;
			if(px_a<0)
			px_a = n;

			if(in_flag[ parent_a[px_a] ]==0)
			/*put parent_a[px_a] into g*/
			{
				side_a++;
				pos=center - side_a;
				tmp_offspring[pos]=parent_a[ px_a];
				in_flag[ tmp_offspring[pos] ] = 1;
			}
			else
				a_flag=0;
		}
		if(b_flag==1)
		{
			px_b = px_b+1;
			if(px_b>n-1)
    			px_b = 0;

			if(in_flag[ parent_b[px_b] ]==0)
			{
				side_b++;
				pos=center+side_b;
				tmp_offspring[pos]=parent_b[ px_b];
				in_flag[ tmp_offspring[pos] ]=1;
			}
			else
				b_flag=0;
		}
	}while(  (a_flag || b_flag)  ); //&& (side_a+side_b+1)<=n  
	step = side_a+side_b+1;
	last=generate_random_permutation(n-step);/*0~(n-step-1) -> step~n-1*/
	for(i=0;i<step;i++)
	{
		offspring[i]=tmp_offspring[ center-side_a+i ];
	}
	for(i=0;i<n-step;i++)
	{		
		for(j=0;j<n;j++)
		{
			if(in_flag[j]==0)
			{
				offspring[ last[i]+ step ]=j;/*0~(n-step-1) -> step~n-1*/
				in_flag[j]=1;
				break;
			}
		}
	}
	offspring[n]=offspring[0];
	free(tmp_offspring);
	free(in_flag);
	free(last);
}