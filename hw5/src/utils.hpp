#ifndef __UTILS_H
#define __UTILS_H

#include "math.h"

template <class T>
inline T max3(T a, T b, T c){
	return fmax(fmax(a, b), c);
}


struct domain{
	double a, b;
};


struct node{
	struct domain *interval;
	struct node *next;
};

#endif
