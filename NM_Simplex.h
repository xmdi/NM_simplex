#ifndef NM_SIMPLEX_H
#define NM_SIMPLEX_H

struct NM_Simplex_History_Iteration{
	double* points;
	double* values;
	double stopping_criteria;
	struct NM_Simplex_History_Iteration* next_iteration;
};

struct NM_Simplex_Problem{
	int num_vars; // number of variables (n)
	double* points; // 1D array encoding (n+1) points in (n) variables
	double* values; // 1D array encoding objective function values for each of (n+1) points
	double (*obj_function)(double[]); // function pointer to objective function
	double tolerance; // flatness tolerance for stopping criteria
	double reflection_coefficient;
	double expansion_coefficient;
	double contraction_coefficient;
	double shrink_coefficient;
	struct NM_Simplex_History_Iteration* initial_iteration;
};

struct NM_Simplex_Problem* setup_simplex(int num_vars, double points[], double (*obj_function)(double[]), double tolerance);
void kill_simplex(struct NM_Simplex_Problem* problem);
void optimize_simplex(struct NM_Simplex_Problem* problem, unsigned char output_flags);

#endif
