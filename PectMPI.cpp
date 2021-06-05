# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "mpi.h"
int main(int argc, char* argv[]);
double f(double x);


int main(int argc, char* argv[])
{
	double end_time;
	double h;
	int i;
	int ierr;
	int m;
	int master = 0;
	int n;
	int process;
	int process_id;
	int process_num;
	double q_global;
	double q_local;
	int received;
	int source;
	double start_time;
	MPI_Status status;
	int tag;
	int target;
	double x;
	double xb[2];
	double x_max = 100;
	double x_min = 0;
	ierr = MPI_Init(&argc, &argv);
	
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
	
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	
	if (process_id == master)
	{
		start_time = MPI_Wtime();
		if (process_num <= 1)
		{
			ierr = MPI_Finalize();
			
		}
	}
	
	
	if (process_id == master)
	{
		for (process = 1; process <= process_num - 1; process++)
		{
			xb[0] = ((double)(process_num - process) * x_min
				+ (double)(process - 1) * x_max)
				/ (double)(process_num - 1);
			xb[1] = ((double)(process_num - process - 1) * x_min
				+ (double)(process)*x_max)
				/ (double)(process_num - 1);
			target = process;
			tag = 1;
			printf("Interval points %f \n", xb[0]);
			ierr = MPI_Send(xb, 2, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
		}
	}
	else
	{
		source = master;
		tag = 1;
		ierr = MPI_Recv(xb, 2, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	}
	
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	
	
	m = 100;
	source = master;
	ierr = MPI_Bcast(&m, 1, MPI_INT, source, MPI_COMM_WORLD);
	
	if (process_id!= master)
	{
		q_local = 0.0;
		for (i = 1; i <= m; i++)
		{
			x = ((double)(2 * m - 2 * i + 1) * xb[0]
				+ (double)(2 * i - 1) * xb[1])
				/ (double)(2 * m);
			q_local = q_local + f(x);
		}
		q_local = q_local * (xb[1] - xb[0]) / (double)(m);
		target = master;
		tag = 2;
		ierr = MPI_Send(&q_local, 1, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
	}
	
	else
	{
		received = 0;
		q_global = 0.0;
		while (received < process_num - 1)
		{
			source = MPI_ANY_SOURCE;
			tag = 2;
			ierr = MPI_Recv(&q_local, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
				&status);
			q_global = q_global + q_local;
			received = received + 1;
		}
	}
	
	if (process_id == master)
	{
		printf("\n");
		printf("Function integral  F (x) = %f\n", q_global);
		printf("Calculation error  %f\n", q_global - 0.624);
		end_time = MPI_Wtime();
		printf("\n");
		printf("Time spent on calculations = %f\n",
			end_time - start_time);
	}
	
	ierr = MPI_Finalize();
	
	
}

double f(double x)
{
	double value;
	value = 2.0 * cos(x) / (2.0 + x * x);
	return value;
}

