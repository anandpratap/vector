#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#define ERROR -1

// 1000        499.697897
// 5000      2499.400123


// For N = 1000:
// Sum = 490.632874
// Time (for p = 1) = 4.96 secs

// For N = 5000:
// Sum = 2485.381836
// Time(for p = 1) = 116.69 secs

void parseInput(int argc, char **argv, int* N){
  // Currently there is only one command line input which is
  // N, the size of the matrix
  
  // exit with ERROR flag if any exceptions are found
  
  if(argc != 2){
    printf("Please provide the size of the matrix!\n");
    exit(ERROR);
  }
  else{
    *N = atoi(argv[1]);
  }
  if(*N <= 0){
    printf("Please provide a valid integer size!\n");
    exit(ERROR);
  }
  else{
    printf("Using N = %d\n", *N);
  }
}

template <class dtype>
dtype sumDiag(int N, dtype A[N][N]){
  dtype sum = 0.0;
  for(int i=0; i<N;i++){
    sum += A[i][i];
  }
  return sum;
}

template <class dtype>
void init(int N, dtype A[N][N]){
  A[0][0:N] = 0.0;
  A[1:N-2][0:N] = 0.5;
  double fac = M_PI/pow(N,2);
  for(int j=0; j<N; j++){
    A[N-1][j] = 5.0*sin(fac*pow(j,2));
  }
}

int main(int argc, char **argv){
  int N;
  parseInput(argc, argv, &N);
 
  double A[N][N], B[N][N+2];
  init(N, A);
  
  int t = 0, tf = 500;

  // local variables used inside the loop
  int jp, jm, jn;

  // stores sum of the diagonal entries
  double sum;

  // print the initial sum
  sum = sumDiag(N, A);
  printf("sum after %d iterations = %f\n",sum, t);

  // entries which remains constant are copied to B


  while(true){
    t++;

    B[0:N][1:N] = A[0:N][0:N];
    // apply boundary condition
    B[0:N][0] = A[0:N][N-1];
    B[0:N][N+1] = A[0:N][0];


    for(int i=1; i<N-1; i++){
      for(int j=0; j<N; j++){
	// j indices of B are inflated by 1
	jp = j+1 + 1;
	jm = j-1 + 1;
	jn = j + 1;
	A[i][j] = 0.125*(B[i][jp] + B[i][jm] + B[i-1][jp] + B[i-1][jm] + B[i+1][jp] + B[i+1][jm] + B[i+1][jn] + B[i-1][jn]);
      }
    }
    // check for number of iterations
    if( t == tf){
      break;
    }
    
  }
  

  // print the final results
  sum = sumDiag(N, A);
  printf("sum after %d iterations = %f\n",sum, t);
  
  return 0;
}
