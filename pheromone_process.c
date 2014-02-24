#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#include"AS.h"
#include"Map.h"
#include"pheromone_process.h"
#include"utilities.h"




void set_relation()
{
	long int i;

	for ( i = 0 ; i < n ; i++ ) {
	    relationship[instance.nodeptr[i].group][(instance.nodeptr[i].subcitynumber-1)]=(instance.nodeptr[i].citynumber-1) ;
//		printf("%d,%d  %d\n",instance.nodeptr[i].group,instance.nodeptr[i].subcitynumber-1,(instance.nodeptr[i].citynumber-1));
	}
}




void generate_cinstance( long int n_map )
{
	tn=n;

	/*
	if((cinstance = malloc(sizeof( problem ))) == NULL)
	{
		printf("Out of memory, exit.");
		exit(1);
    }
	*/

	strcpy(cinstance.name, instance.name);
	cinstance.n = instance.n;
	cinstance.n_near = instance.n_near;
	cinstance.optimum = instance.optimum;

}

void set_map(long int i_map)
{
	long int i,j;
	j=0;
	printf("set map: %d\n",i_map);

	for ( i = 0 ; i < tn ; i++ )
		{
//			printf("\n%ld %ld %lf %lf %ld %ld",i, cinstance.nodeptr[i].citynumber, cinstance.nodeptr[i].x, cinstance.nodeptr[i].y, cinstance.nodeptr[i].group, cinstance.nodeptr[i].subcitynumber );
		}

	for ( i = 0 ; i < tn ; i++ ) {
	    if(cinstance.nodeptr[i].group==i_map)
		{
			instance.nodeptr[j]=cinstance.nodeptr[i];
			printf("%d\n",instance.nodeptr[j].citynumber);
			j++;
		}
		else;
		
	}
	n=j;
	printf("n=%d\n",n);

	///////////change by 2014.1.24////////////////////
	i=i_map;

	map[i].pheromone = generate_double_matrix( n, n );
	map[i].total = generate_double_matrix( n, n );

	map[i].best_in_try = calloc(max_tries, sizeof(double));
	map[i].best_found_at = calloc(max_tries, sizeof(long int));
	map[i].time_best_found = calloc(max_tries, sizeof(double));
	map[i].time_total_run = calloc(max_tries, sizeof(double));

	map_allocate_ants(i);
	//////////////here ends///////////////////////////


	map[i_map].mn=j;
	instance.distance = compute_distances();
	instance.nn_list  = compute_nn_lists();

	compute_heuristic();

/*	for ( i = 0 ; i < n ; i++ ) {
	printf("%d %d\n",instance.nodeptr[i].subcitynumber,instance.nodeptr[i].citynumber);
	}
*/

	
}


void indi_map(long int x_num,long int y_num) 
/*    
      FUNCTION: divide the cities to several groups
      INPUT:    the dividing number of x-axis and y axis
      OUTPUT:   get change the [group] and [subcitynumber] in [point]
      COMMENTS: none
*/
{
	
	long int middle_size,small_size;
	
	long int i,j,k,l;
	long int min_i;
	double min_x;
	double min_y;
	
	tn=n;
	
	middle_size=(tn/y_num);
	small_size=(middle_size/x_num);

	printf("\nMap dividing begins\n");
	
	
	for(i=0;i<tn;i++)
	{
		instance.nodeptr[i].group=-1;
		instance.nodeptr[i].citynumber=i+1;
	}
	
	
	for( k=0 ; k<y_num ; k++ )								//get middle group for y_num times
	{
		if( k == y_num-1 )									//when get last group, all the last cities belong to this group
		{
			for( i=0 ; i<tn ; i++ )
			{
				if( instance.nodeptr[i].group != -1 )
				{;}
				else
				{
					instance.nodeptr[i].group=k*x_num;
				}
			}
		}
		
		for( j=0 ; j<middle_size ; j++ )					//every middle group have middle_size cities
		{
			min_y=INFTY;
			for( i=0 ; i<tn ; i++ )
			{
				if(instance.nodeptr[i].group!=-1)
				{;}
				else if((instance.nodeptr[i].y)<=min_y)
				{
					min_y=instance.nodeptr[i].y;
					min_i=i;
				}
				else
				{;}
			}
			instance.nodeptr[min_i].group=x_num*k;
		}
		
	}


	
	for( k=0 ; k<y_num ; k++ )								//divid every middle groups
	{
		for(i=0;i<tn;i++)
		{
			if( instance.nodeptr[i].group == x_num*k )
			{
				instance.nodeptr[i].group=-1;
			}
			else
			{;}
		
		}

		for( l=0 ; l<x_num ; l++)
		{
			if( l == x_num-1 )
			{
				for( i=0 ; i<tn ; i++ )
				{
					if( instance.nodeptr[i].group != -1 )
					{;}
					else
					{
						instance.nodeptr[i].group=k*x_num+l;
					}
				}
			}
			
			for( j=0 ; j<small_size ; j++)					//individ each middle group into small group
			{
				min_x=INFTY;
				for( i=0 ; i<tn ; i++ )
				{
					if(instance.nodeptr[i].group!=-1)
					{;}
					else if((instance.nodeptr[i].x)<=min_x)
					{
						min_x=instance.nodeptr[i].x;
						min_i=i;
					}
					else
					{;}
				}
				instance.nodeptr[min_i].group=x_num*k+l;
			}
			
		}
		
	}

	for(k=0 ; k<x_num*y_num ; k++)
	{	
		j=1;
		for( i=0 ; i<tn ; i++ )
		{
			if(instance.nodeptr[i].group == k)
			{
				instance.nodeptr[i].subcitynumber=j;
				j++;
			}

			else;
		}
	}


	for ( i = 0 ; i < n ; i++ ) 
	{
//		printf("\n%ld %lf %lf %ld %ld", instance.nodeptr[i].citynumber, instance.nodeptr[i].x, instance.nodeptr[i].y, instance.nodeptr[i].group, instance.nodeptr[i].subcitynumber );
	}
	printf("\n\nMap dividing ends\n");
	

	for(i=0;i<tn;i++)
	{
		cinstance.nodeptr[i]=instance.nodeptr[i];
	}


}

void pheromone_divide( long int island_num)
{
//	long int is;
	long int i,j;
	long int inew,jnew;

	printf("cross over the pheromone from each map\n");


		n=map[island_num].mn;
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
				inew=relationship[island_num][i];
				jnew=relationship[island_num][j];
				map[island_num].pheromone[i][j]=phex[inew][jnew];
//				printf("inew:%d jnew:%d is:%d i:%d j:%d\n",inew,jnew,is,i,j);
//				printf("phe: %0.15f\n",phe[inew][jnew]);
			}
		}

	map_compute_total_information(island_num);
	
}


void redivide(long int rx_num, long int ry_num)
/*    
      FUNCTION: redivide the cities to several groups
      INPUT:    the redividing number of rx_axis and ry_axis
      OUTPUT:   get change the [group] and [subcitynumber] in [point]
				put the pheromone from total map to divided maps
      COMMENTS: none
*/
{


	/****************************************************************/
	/*the redivide start*/
	/****************************************************************/
	
	indi_map(rx_num,ry_num);
	
	/****************************************************************/
	/* the redivide ends*/
	/****************************************************************/
		
		
	set_relation();
		
		



		
		
}

