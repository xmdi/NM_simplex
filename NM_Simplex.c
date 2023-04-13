#include <stdlib.h>
#include <stdio.h>
#include "NM_Simplex.h"

void populate_simplex_iteration_history(struct NM_Simplex_Problem* problem, struct NM_Simplex_History_Iteration* iteration, double stopping_criteria){
	double* points=malloc(problem->num_vars*(problem->num_vars+1)*sizeof(double));
	double* values=malloc((problem->num_vars+1)*sizeof(double));
	for (int i=0; i<(problem->num_vars+1); i++){
		for (int j=0; j<problem->num_vars; j++){
			points[i*(problem->num_vars)+j]=problem->points[i*(problem->num_vars)+j];
		}
		values[i]=problem->values[i];
	}
	iteration->points=points;
	iteration->values=values;
	iteration->stopping_criteria=stopping_criteria;	
	iteration->next_iteration=NULL;
}

struct NM_Simplex_Problem* setup_simplex(int num_vars,double points[], double (*obj_function)(double[]), double tolerance){
	struct NM_Simplex_Problem* problem=malloc(sizeof(struct NM_Simplex_Problem));
	problem->num_vars=num_vars;
	problem->tolerance=tolerance;
	double* pts=malloc(num_vars*(num_vars+1)*sizeof(double));
	for (int i=0; i<(num_vars+1); i++){
		for (int j=0; j<num_vars; j++){
			pts[i*(num_vars)+j]=points[i*(num_vars)+j];
		}
	}
	double* values=calloc((num_vars+1),sizeof(double));
	problem->values=values;
	problem->points=pts;
	problem->obj_function=obj_function;
	problem->reflection_coefficient=1.0;
	problem->expansion_coefficient=2.0;
	problem->contraction_coefficient=0.5;
	problem->shrink_coefficient=0.5;
	
	struct NM_Simplex_History_Iteration* initial_iteration=malloc(sizeof(struct NM_Simplex_History_Iteration));
	problem->initial_iteration=initial_iteration;

	return problem;
}

void kill_simplex_history_iteration(struct NM_Simplex_History_Iteration* iteration){
	if (iteration->next_iteration!=NULL){
		kill_simplex_history_iteration(iteration->next_iteration);
	}
	free(iteration->points);
	free(iteration->values);
	free(iteration);
}

void kill_simplex(struct NM_Simplex_Problem* problem){
	free(problem->points);
	free(problem->values);
	kill_simplex_history_iteration(problem->initial_iteration);
	free(problem);
}

void optimize_simplex(struct NM_Simplex_Problem* problem, unsigned char output_flags){

	// compute initial y_i's
	for (int i=0; i<(problem->num_vars+1); i++){
		problem->values[i]=problem->obj_function(&problem->points[i*(problem->num_vars)]);
	}
	int iteration=0;
	int fevals=problem->num_vars+1;
	populate_simplex_iteration_history(problem,problem->initial_iteration,1.0/0.0);
	struct NM_Simplex_History_Iteration* current_iteration=problem->initial_iteration;

	double stopping_criteria=1+problem->tolerance;
	double tolerance=problem->tolerance*problem->tolerance*problem->num_vars; // avoid sqrt and division every iteration

	if (output_flags&1){
		printf("\x1B[1;32m> Nelder-Mead Simplex Start >\n\x1B[0m");
		if ((~output_flags)&2){
			printf("\x1B[1;34m\titer\tfevals\tfmin\t\tstopping_criteria (<=%.20f)\n\x1B[0m",tolerance);
		}
	}
	
	while (stopping_criteria>tolerance){
		iteration++;
		// identify the best & worst P_i
		double y_w=problem->values[0];
		double y_b=problem->values[0];
		int i_w=0;
		int i_b=0;
		for (int i=1; i<(problem->num_vars+1); i++){
			if (problem->values[i]>y_w){
				y_w=problem->values[i];
				i_w=i;
			}
			else if (problem->values[i]<y_b){
				y_b=problem->values[i];
				i_b=i;
			}
		}
		// compute centroid P_bar of non-worst point cloud
		double* P_bar=calloc(problem->num_vars,sizeof(double));
		for (int i=0; i<(problem->num_vars+1); i++){
			for (int j=0; j<(problem->num_vars); j++){
				if (i!=i_w){
					P_bar[j]+=problem->points[i*(problem->num_vars)+j];
				}
			}
		}
		for (int i=0; i<(problem->num_vars); i++){
			P_bar[i]=P_bar[i]/(problem->num_vars);
		}
		// REFLECT STEP
		double* P_star=malloc(problem->num_vars*sizeof(double));
		for (int i=0; i<problem->num_vars; i++){
			P_star[i]=(1.0+problem->reflection_coefficient)*P_bar[i]-problem->reflection_coefficient*problem->points[i_w*problem->num_vars+i];
		}
		// evaluate y_star at P_star
		double y_star=problem->obj_function(P_star);
		fevals++;
		if (y_star<problem->values[i_b]){
			// EXPAND STEP
			double* P_starstar=malloc(problem->num_vars*sizeof(double));
			for (int i=0; i<problem->num_vars; i++){
				P_starstar[i]=(1.0-problem->expansion_coefficient)*P_bar[i]-problem->expansion_coefficient*P_star[i];
			}
			double y_starstar=problem->obj_function(P_starstar);
			fevals++;
			if (y_starstar<problem->values[i_b]){	
				// replace P_w with P_starstar
				problem->values[i_w]=y_starstar;
				for (int i=0; i<problem->num_vars; i++){
					problem->points[i_w*problem->num_vars+i]=P_starstar[i];
				}
			}
			else {
				// replace P_w with P_star
				problem->values[i_w]=y_star;
				for (int i=0; i<problem->num_vars; i++){
					problem->points[i_w*problem->num_vars+i]=P_star[i];
				}
			}
			free(P_starstar);
		}
		else {
			// CONTRACT STEP
			int rank_of_y_star=0;
			for (int i=0; i<(problem->num_vars+1); i++){
				if ((i!=i_w)&&(y_star>problem->values[i])){
					rank_of_y_star++;
				}
			}
			if (rank_of_y_star==problem->num_vars){
				if (y_star<=problem->values[i_w]){
					// replace P_w with P_star
					problem->values[i_w]=y_star;
					for (int i=0; i<problem->num_vars; i++){
						problem->points[i_w*problem->num_vars+i]=P_star[i];
					}
				}
				double* P_starstar=malloc(problem->num_vars*sizeof(double));
				for (int i=0; i<problem->num_vars; i++){
					P_starstar[i]=(1.0-problem->contraction_coefficient)*P_bar[i]-problem->contraction_coefficient*problem->points[i_w*problem->num_vars+i];
				}
				double y_starstar=problem->obj_function(P_starstar);
				fevals++;
				if (y_starstar>problem->values[i_w]){
					// SHRINK STEP
					for (int i=0; i<(problem->num_vars+1); i++){
						if (i!=i_b){
							for (int j=0; j<problem->num_vars; j++){
								problem->points[i*problem->num_vars+j]=0.5*(problem->points[i*problem->num_vars+j]+problem->points[i_b*problem->num_vars+j]);
							}
							problem->values[i]=problem->obj_function(&problem->points[i*problem->num_vars]);
							fevals++;
						}
					}
				}
				else {
					// replace P_w with P_starstar
					problem->values[i_w]=y_starstar;
					for (int i=0; i<problem->num_vars; i++){
						problem->points[i_w*problem->num_vars+i]=P_starstar[i];
					}
				}
				free(P_starstar);
			}
			else {
				// replace P_w with P_star
				problem->values[i_w]=y_star;
				for (int i=0; i<problem->num_vars; i++){
					problem->points[i_w*problem->num_vars+i]=P_star[i];
				}
			}
		}

		free(P_bar);
		free(P_star);

		// compute stopping criteria
		stopping_criteria=0;
		double y_bar=0;
		for (int i=0; i<(problem->num_vars+1); i++){
			y_bar+=problem->values[i];
		}
		y_bar=y_bar/(problem->num_vars+1);
		for (int i=0; i<(problem->num_vars+1); i++){
			stopping_criteria+=(problem->values[i]-y_bar)*(problem->values[i]-y_bar);
		}

		// recompute y_b and i_b
		y_b=problem->values[0];
		i_b=0;
		for (int i=1; i<(problem->num_vars+1); i++){
			if (problem->values[i]<y_b){
				y_b=problem->values[i];
				i_b=i;
			}
		}

		if (output_flags&1){
			if (output_flags&2){
				printf("\x1B[1;34m\titer\tfevals\tfmin\t\tstopping_criteria (<=%.20f)\n\x1B[0m",tolerance);
			}
			printf("\t%d\t%d\t%f\t%.20f\n",iteration,fevals,y_b,stopping_criteria);
			if (output_flags&2){
				printf("\n\t\tP_best = (");
				for (int i=0; i<problem->num_vars; i++){
					if (i){
						printf(",");
					}
					printf("%f",problem->points[i_b*problem->num_vars+i]);
				}
				printf(")\n\n");
			}
		}

		if (output_flags&4){
			struct NM_Simplex_History_Iteration* next_iteration=malloc(sizeof(struct NM_Simplex_History_Iteration));
			current_iteration->next_iteration=next_iteration;
			current_iteration=next_iteration;
			populate_simplex_iteration_history(problem,current_iteration,stopping_criteria);
		}

	}
	
	if (output_flags&1){
		printf("\x1B[1;31m< Nelder-Mead Simplex End <\n\x1B[0m");
	}

}
