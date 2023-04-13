#include <stdlib.h>
#include <stdio.h>
#include "NM_Simplex.h"

// sample analytical objective function
// you can replace this with anything you want to optimize
double treasure_seeking(double point[]){ 
	return 0.006*(point[0]-1)*(point[0]-1)+0.009*(point[1]-2)*(point[1]-2)+.005*point[0]*point[1]-10;
}

int main(void){
	int num_vars=2;
	double points[3*2]={-30,30,  // P_0
			-20,20,  // P_1
			-10,30}; // P_2
	double (*obj_function)(double[])=&treasure_seeking;
	double tolerance=1e-10;
	unsigned char output_flags=0b001; // bit 1 (rightmost): print convergence history during run
									  // bit 2: if bit 1 high, also track simplex minimum vertex
									  // bit 3: save simplex history to data structure for later use
	struct NM_Simplex_Problem* problem=setup_simplex(num_vars,points,obj_function,tolerance);
	optimize_simplex(problem,output_flags);

	// below code prints out the points of the "converged" simplex, you don't need this.
	printf("\x1B[1;33m\nConverged Simplex:\n");
	for (int j=0;j<(problem->num_vars+1);j++){
			printf("(");
			for (int k=0;k<problem->num_vars;k++){
				if (k){
					printf(",");
				}
				printf("%f",problem->points[j*problem->num_vars+k]);
			}
			printf(",%f)\n",problem->values[j]);
		}
	printf("\n\x1B[0m");

	// below code prints out the full simplex history linked list, you don't need this.
	// this is only possible if the 3rd bit in the output_flags is high
	if (output_flags&4){
		printf("\nSimplex History:\n");
		int i=0;
		struct NM_Simplex_History_Iteration* iter=problem->initial_iteration;
		while (iter!=NULL){
			printf("iteration #%d\t",i++);
			for (int j=0;j<(problem->num_vars+1);j++){
				printf("(");
				for (int k=0;k<problem->num_vars;k++){
					if (k){
						printf(",");
					}
					printf("%f",iter->points[j*problem->num_vars+k]);
				}
				printf(",%f)\n",iter->values[j]);
			}
			printf("\n");
			iter=iter->next_iteration;
		}
	}

	kill_simplex(problem);
	return 0;
}
