#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

template <class T>
T func(T x){
	T fval = 0.0;
	T local_sum;
    
	int n = 10;
    
	for(int i=n; i>=1; i--){
		local_sum = 0.0;
	
		for(int j=i; j>=1; j--){
			local_sum += pow(x + j, -3.1);
		}
	
		fval += sin(x + local_sum)/pow(1.2, i);
	}
    
	return fval;
}

#endif
