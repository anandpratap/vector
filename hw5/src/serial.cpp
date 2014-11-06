#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "assert.h"

#include "functions.hpp"
#include "utils.hpp"
#include "queue.hpp"

#define A 1.0
#define B 100.0
#define TOL 1E-6
#define S 12.0
#define VERYSMALL -1E10

void find_maxima(struct queue *q, double *maxval){
	double s = S;
	double tol = TOL;
	double a, b, check;
	double fa, fb;
	
	unsigned int max_len = 0;
	struct node *denode;
	while(q->length > 0){
		denode = dequeue(q);
		a = denode->interval->a;
		b = denode->interval->b;
		delete(denode);
		fa = func(a);
		fb = func(b);
		
		*maxval = max3(fa, fb, *maxval);
		check = 0.5*(fa + fb + (b-a)*s);
		if(check > *maxval + tol){
			enqueue(q, a, 0.5*(a+b));
			enqueue(q, 0.5*(a+b), b);
		}
		max_len = fmax(max_len, q->length);
	}
	printf("maximum length of queue = %d\n", max_len);
}

void test(void){
	double a, b;
	a = 1.0;
	b = 1000.0;
	
	struct queue *workqueue;
	workqueue = (struct queue*) malloc(sizeof(struct queue));

	assert(workqueue->first == NULL);
	assert(workqueue->last == NULL);
	assert(workqueue->length == 0);


	// Add two nodes in queue
	enqueue(workqueue, a, b);
	enqueue(workqueue, a, b-10);

	assert(workqueue->length == 2);
	display(workqueue);

	struct node *dnode;
	// remove one and this should be the first added element
	dnode = dequeue(workqueue);
	assert(workqueue->length == 1);
	assert(fabs(dnode->interval->a - a) <1E-15);
	assert(fabs(dnode->interval->b - b) <1E-15);
    
	// remove another 
	dnode = dequeue(workqueue);

	assert(workqueue->first == NULL);
	assert(workqueue->last == NULL);
	assert(workqueue->length == 0);

	display(workqueue);
}

int main(void){
	// run tests
	test();

	// setup problem of maxima
	double a, b;
	a = A;
	b = B;
	struct queue *workqueue;
	workqueue = new queue;
	enqueue(workqueue, a, b);
	
	double maxval = VERYSMALL;
	find_maxima(workqueue, &maxval);
	printf("maxima = %1.16f\n", maxval);
	
	delete(workqueue);
	return 0;
}
