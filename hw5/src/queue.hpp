#ifndef __QUEUE_H
#define __QUEUE_H
#include "assert.h"
#include "utils.hpp"


// enqueue(*q, a, b)
// *node = dequeue(*q)
// n = get_length(*q)
// remove from first add to last
// while adding update current last->next to new and last to new

struct queue{
	struct node *first;
	struct node *last;
	unsigned int length;
	queue(){
		first = NULL;
		last = NULL;
		length = 0;
	};
};

unsigned int get_length(struct queue *q){
	return q->length;
}


void enqueue(struct queue *q, double a, double b){
	assert(b > a);
	struct node *temp = new node;
	struct domain *interval = new domain;
	interval->a = a;
	interval->b = b;
	temp->interval = interval;
	temp->next = NULL;
    
	// if the queue is empty
	if(q->first == NULL){
		q->first = temp;
		q->last = temp;
	}
	// if the queue is not empty
	else{
		q->last->next = temp;
		q->last = temp;
	}
	// increment length counter
	q->length++;
}

struct node* dequeue(struct queue *q){
	// check for empty queue
	assert(q->first != NULL);

	// create a node to return
	struct node *temp;
	temp = q->first;
	q->first = q->first->next;
    
	// decrement length counter
	q->length--;
    
	// if the resulting queue is empty set last pointer to null as well
	if(q -> length == 0){
		q->last = NULL;
	}

	return temp;
}

void display(struct queue *q){
	// display starting from first to last
	struct node *root;
	if(q->first == NULL){
		printf("The queue is empty!\n");
	}
	else{
		root = q->first;
		while(true){ 
			printf("%f %f\n", root->interval->a, root->interval->b);
			if(root->next == NULL){
				break;
			}
			root = root->next;
		}
	}
}


#endif
