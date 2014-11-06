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

 void find_maxima(struct queue *q, double *maxval){
	 int id = omp_get_thread_num();
	 printf("Thread id %d \n", id);
	 int status;
	 double a, b;
	 unsigned int max_len = 0;
	 struct node *denode;
	 while(q->length > 0){
		 denode = dequeue(q);
		 a = denode->interval->a;
		 b = denode->interval->b;
		 
		 status = checkmax(a, b, maxval);
		 delete(denode);
		 if(!status){
			 enqueue(q, a, 0.5*(a+b));
			 enqueue(q, 0.5*(a+b), b);
		 }
		 max_len = fmax(max_len, q->length);
	 }
	 printf("maximum length of queue = %d\n", max_len);
}

int main(void){

	// setup problem of maxima
	double a, b;
	a = A;
	b = B;
	struct queue *workqueue;
	workqueue = new queue;
	enqueue(workqueue, a, b);
	
	double maxval = VERYSMALL;
#pragma omp parallel num_threads(1)
	{
		find_maxima(workqueue, &maxval);
	}
	printf("maxima = %1.16f\n", maxval);
	
	delete(workqueue);
	return 0;
}
