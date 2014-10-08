#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include "mpi.h"

#define ERROR -1

int get_rank(int x, int y, int z){
  // return processor rank based on x and y coordinates
  return x*z + y;
}

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
dtype sumDiag(int m, dtype A[m][m], int x, int y, int size){
  // return sum of the diagonal entries of the original matrices contains in processor (x,y)
  int io = x*m;
  int jo = y*m;
  dtype sum = 0.0;
  for(int i=0; i<m;i++){
    for(int j=0; j<m;j++){
      if( io + i == jo + j){
	sum += A[i][i];
      }
    }
  }
  return sum;
}

template <class dtype>
void init(int m, dtype A[m][m], int x, int y, int size){
  int N = sqrt(size)*m;
  int jstart = y*m;
  dtype fac = M_PI/pow(N,2);
  // bottom boundary
  if(x == 0){
    A[0][0:m] = 0.0;
  }
  else{
    A[0][0:m] = 0.5;
  }
  // interior points
  A[1:m-2][0:m] = 0.5;
  if(x == (int)sqrt(size) - 1){
    for(int j=0; j<m; j++){
      A[m-1][j] = 5.0*sin(fac*pow(jstart+j,2));
    }
  }
  else{
    for(int j=0; j<m; j++){
      A[m-1][j] = 0.5;
    }
  }
}


int main(int argc, char **argv){
  int N;
  // get the size of the matrix
  parseInput(argc, argv, &N);

  // MPI init stuff
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int x, y;
  x = rank / (int) sqrt(size);
  y = rank % (int) sqrt(size);
  
  int xp, yp, xm, ym;
  int z = (int) sqrt(size);
  xp = (x + 1 >= z) ? 0 : x + 1;
  xm = (x - 1 < 0) ? z - 1 : x - 1;
  yp = (y + 1 >= z) ? 0 : y + 1;
  ym = (y - 1 < 0) ? z - 1 : y - 1;
  
  int south, southeast, east, northeast, north, northwest, west, southwest;
  south = get_rank(xm, y, z);
  southeast = get_rank(xm, yp, z);
  east = get_rank(x, yp, z);
  northeast = get_rank(xp, yp, z);
  north = get_rank(xp, y, z);
  northwest = get_rank(xp, ym, z);
  west = get_rank(x, ym, z);
  southwest = get_rank(xm, ym, z);
    
  // size of the submatrix
  int m;
  m = N/(int)sqrt(size);
  
  double A[m][m], B[m+2][m+2];
  init(m, A, x, y, size);
  int t = 0, tf = 500;

  // local variables used inside the loop
  int jp, jm, jn, ip, im, in;
  int istart, iend;

  // stores sum of the diagonal entries
  double sum, total_sum;
  
  // print the initial sum
  sum = sumDiag(m, A, x, y, size);
  MPI_Reduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank == 0){
    printf("sum after %d iterations = %10.10f\n",total_sum, t);
  }
  
  // entries which remains constant are copied to B
  double send_buffer[m], recv_buffer[m];
  MPI_Status status;
  // startime

  MPI_Barrier(MPI_COMM_WORLD);
  double starttime = MPI_Wtime();
  while(true){
    t++;
    B[1:m][1:m] = A[:][:];
    if(size != 1){
      // send to south and receive from north
      send_buffer[0:m] = A[0][0:m];
      MPI_Send(&send_buffer, m, MPI_DOUBLE, south, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, m, MPI_DOUBLE, north, 1, MPI_COMM_WORLD, &status);
      B[m+1][1:m] = recv_buffer[0:m];
      
      // send to north receive from south
      send_buffer[0:m] = A[m-1][0:m];
      MPI_Send(&send_buffer, m, MPI_DOUBLE, north, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, m, MPI_DOUBLE, south, 1, MPI_COMM_WORLD, &status);
      B[0][1:m] = recv_buffer[0:m];
      
      // send to east and receive from west
      send_buffer[0:m] = A[0:m][m-1];
      MPI_Send(&send_buffer, m, MPI_DOUBLE, east, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, m, MPI_DOUBLE, west, 1, MPI_COMM_WORLD, &status);
      B[1:m][0] = recv_buffer[0:m];
      
      // send to west and receive from east
      send_buffer[0:m] = A[0:m][0];
      MPI_Send(&send_buffer, m, MPI_DOUBLE, west, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, m, MPI_DOUBLE, east, 1, MPI_COMM_WORLD, &status);
      B[1:m][m+1] = recv_buffer[0:m];
      
      // send to south east and receive from north west
      send_buffer[0] = A[0][m-1];
      MPI_Send(&send_buffer, 1, MPI_DOUBLE, southeast, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, northwest, 1, MPI_COMM_WORLD, &status);
      B[m+1][0] = recv_buffer[0];
      
    // send to north east and receive from south west
      send_buffer[0] = A[m-1][m-1];
      MPI_Send(&send_buffer, 1, MPI_DOUBLE, northeast, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, southwest, 1, MPI_COMM_WORLD, &status);
      B[0][0] = recv_buffer[0];
      
      // send to north west and receive from south east
      send_buffer[0] = A[m-1][0];
      MPI_Send(&send_buffer, 1, MPI_DOUBLE, northwest, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, southeast, 1, MPI_COMM_WORLD, &status);
      B[0][m+1] = recv_buffer[0];
      
      // send to south west and and receive from north east
      send_buffer[0] = A[0][0];
      MPI_Send(&send_buffer, 1, MPI_DOUBLE, southwest, 1, MPI_COMM_WORLD);
      MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, northeast, 1, MPI_COMM_WORLD, &status);
      B[m+1][m+1] = recv_buffer[0];
    }
    
    else{
      B[m+1][1:m] = A[0][0:m];
      B[0][1:m] = A[m-1][0:m];
      B[1:m][0] = A[0:m][m-1];
      B[1:m][m+1] = A[0:m][0];
      B[m+1][0] = A[0][m-1];
      B[0][0] = A[m-1][m-1];
      B[0][m+1] = A[m-1][0];
      B[m+1][m+1] = A[0][0];
    }
    
    
    if(x == 0){
      A[1:m-1][0:m] = 0.125*(B[2:m-1][2:m] + B[2:m-1][0:m] + B[1:m-1][2:m] + B[1:m-1][0:m] + B[3:m-1][2:m] + B[3:m-1][0:m] + B[3:m-1][1:m] + B[1:m-1][2:m]);
    }
    else if(x == z - 1){
      A[0:m-1][0:m] = 0.125*(B[1:m-1][2:m] + B[1:m-1][0:m] + B[0:m-1][2:m] + B[0:m-1][0:m] + B[2:m-1][2:m] + B[2:m-1][0:m] + B[2:m-1][1:m] + B[0:m-1][1:m]);
    }
    else{
      A[0:m][0:m] = 0.125*(B[1:m][2:m] + B[1:m][0:m] + B[0:m][2:m] + B[0:m][0:m] + B[2:m][2:m] + B[2:m][0:m] + B[2:m][1:m] + B[0:m][1:m]);
    }
    
    if( t == tf){
      break;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double endtime = MPI_Wtime();
  if(rank == 0)
    printf("time taken = %10.10f\n", endtime-starttime);
  
  // print the final results
  sum = sumDiag(m, A, x, y, size);
  MPI_Reduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank == 0){
    printf("sum after %d iterations = %10.10f m=%d\n",total_sum, t, m);
  }
  MPI_Finalize();
  return 0;
}
