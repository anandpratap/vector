#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "assert.h"

#include "functions.hpp"
#include "utils.hpp"
#include "queue.hpp"
#include "omp.h"

#define A 1.0
#define B 100.0
#define TOL 1E-6
#define S 12.0
#define VERYSMALL -1E10

int NUM_THREADS;

// maximum length of queue = 116144
// maxima = 4.1924720855452948

int checkmax(double a, double b, double *maxval){
	double s = S;
	double tol = TOL;
	double fa, fb, check; 
	double localmaxval;

#pragma omp atomic read
	localmaxval = *maxval;


	fa = func(a);
	fb = func(b);
	localmaxval = max3(fa, fb, localmaxval);

#pragma omp atomic write
	*maxval = localmaxval;

	check = 0.5*(fa + fb + (b-a)*s);
	
	if(check > *maxval + tol){
		return 0;
	}

	return 1;

}

int get_total_length(struct queue* workqueue){
	int num_threads = NUM_THREADS;
	int size = 0;
	for(int i = 0; i< num_threads; i++){
		size += workqueue[i].length;
	}
	
	return size;
}

int get_min_thread(struct queue* workqueue){
	int minval= 1000000;
	int id;
	for(int i = 0; i < NUM_THREADS; i++){
		if(workqueue[i].length < minval){
			minval = workqueue[i].length;
			id = i;
		}
	}
	return id;
}

void find_maxima(struct queue *workqueue, double *maxval, int *run){
	int id = omp_get_thread_num();
	struct queue *q;
	q = &workqueue[id];
	int minid;
	printf("Thread id %d \n", id);
	
	int status = 1;
	double a, b;
	unsigned int max_len = 0;
	struct node *denode;
	while(1){
		if(q->length != 0){
#pragma omp atomic write
			*run = 1;
			denode = dequeue(q);
			a = denode->interval->a;
			b = denode->interval->b;
			
			status = checkmax(a, b, maxval);
			delete(denode);
			//printf("%d %d\n", q->length, id);
			if(!status){
#pragma omp critical
				{
					minid = get_min_thread(workqueue);
					enqueue(&workqueue[minid], a, 0.5*(a+b));
					enqueue(&workqueue[minid], 0.5*(a+b), b);
				}
			}
			max_len = fmax(max_len, q->length);
#pragma omp atomic write
			*run = 0;

		}
#pragma omp barrier
		if(get_total_length(workqueue) == 0 && *run == 0){
			break;
		}
	}
	printf("maximum length of queue = %d\n", max_len);
}

int main(int argc, char **argv){
	NUM_THREADS = atoi(argv[1]);
	// setup problem of maxima
	double a, b;
	a = A;
	b = B;
	struct queue *workqueue;
	workqueue = new queue[NUM_THREADS];
	enqueue(&workqueue[0], a, b);
	int run=0;
	double maxval = VERYSMALL;
	double start_time = omp_get_wtime();
#pragma omp parallel num_threads(NUM_THREADS)
	{
		find_maxima(workqueue, &maxval, &run);
	}
	double time = omp_get_wtime() - start_time;
	printf("maxima = %1.6f\n", maxval);
	printf("Time Taken= %10.10f\n", time);
	
	delete(workqueue);
	return 0;
}
