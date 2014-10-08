#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include "mpi.h"

#define ERROR -1

inline int get_rank(int x, int y, int z){
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
    printf("sum after %d iterations = %10.10f\n", t, total_sum);
  }
  
  // entries which remains constant are copied to B
  double send_buffer_0[m], recv_buffer_0[m];
  double send_buffer_1[m], recv_buffer_1[m];
  double send_buffer_2[m], recv_buffer_2[m];
  double send_buffer_3[m], recv_buffer_3[m];
  double send_buffer_4[1], recv_buffer_4[1];
  double send_buffer_5[1], recv_buffer_5[1];
  double send_buffer_6[1], recv_buffer_6[1];
  double send_buffer_7[1], recv_buffer_7[1];
  MPI_Request s[8], r[8];

  MPI_Status status[8];
  // startime
  int flag;

  MPI_Barrier(MPI_COMM_WORLD);
  double starttime = MPI_Wtime();
  while(true){
    t++;
    
    if(size != 1){
      // send to south and receive from north
      send_buffer_0[0:m] = A[0][0:m];
      MPI_Isend(&send_buffer_0, m, MPI_DOUBLE, south, 1, MPI_COMM_WORLD, &s[0]);
            
      send_buffer_1[0:m] = A[m-1][0:m];
      MPI_Isend(&send_buffer_1, m, MPI_DOUBLE, north, 1, MPI_COMM_WORLD, &s[1]);
      
      send_buffer_2[0:m] = A[0:m][m-1];
      MPI_Isend(&send_buffer_2, m, MPI_DOUBLE, east, 1, MPI_COMM_WORLD, &s[2]);
      
      send_buffer_3[0:m] = A[0:m][0];
      MPI_Isend(&send_buffer_3, m, MPI_DOUBLE, west, 1, MPI_COMM_WORLD, &s[3]);
      
      send_buffer_4[0] = A[0][m-1];
      MPI_Isend(&send_buffer_4, 1, MPI_DOUBLE, southeast, 1, MPI_COMM_WORLD, &s[4]);
      
      send_buffer_5[0] = A[m-1][m-1];
      MPI_Isend(&send_buffer_5, 1, MPI_DOUBLE, northeast, 1, MPI_COMM_WORLD, &s[5]);

      send_buffer_6[0] = A[m-1][0];
      MPI_Isend(&send_buffer_6, 1, MPI_DOUBLE, northwest, 1, MPI_COMM_WORLD, &s[6]);
      
      send_buffer_7[0] = A[0][0];
      MPI_Isend(&send_buffer_7, 1, MPI_DOUBLE, southwest, 1, MPI_COMM_WORLD, &s[7]);

      MPI_Irecv(&recv_buffer_0, m, MPI_DOUBLE, north, 1, MPI_COMM_WORLD, &r[0]);
      // send to north receive from south
      MPI_Irecv(&recv_buffer_1, m, MPI_DOUBLE, south, 1, MPI_COMM_WORLD, &r[1]);
      // send to east and receive from west
      MPI_Irecv(&recv_buffer_2, m, MPI_DOUBLE, west, 1, MPI_COMM_WORLD, &r[2]);
      // send to west and receive from east
      MPI_Irecv(&recv_buffer_3, m, MPI_DOUBLE, east, 1, MPI_COMM_WORLD, &r[3]);
      // send to south east and receive from north west
      MPI_Irecv(&recv_buffer_4, 1, MPI_DOUBLE, northwest, 1, MPI_COMM_WORLD, &r[4]);
      // send to north east and receive from south west
      MPI_Irecv(&recv_buffer_5, 1, MPI_DOUBLE, southwest, 1, MPI_COMM_WORLD, &r[5]);
      // send to north west and receive from south east
      MPI_Irecv(&recv_buffer_6, 1, MPI_DOUBLE, southeast, 1, MPI_COMM_WORLD, &r[6]);
      // send to south west and and receive from north east
      MPI_Irecv(&recv_buffer_7, 1, MPI_DOUBLE, northeast, 1, MPI_COMM_WORLD, &r[7]);


    }

    B[1:m][1:m] = A[:][:];
    // update interior points
    A[1:m-2][1:m-2] = 0.125f*(B[2:m-2][3:m-2] + B[2:m-2][1:m-2] + B[1:m-2][3:m-2] + B[1:m-2][1:m-2] + B[3:m-2][3:m-2] + B[3:m-2][1:m-2] + B[3:m-2][2:m-2] + B[1:m-2][2:m-2]);
    
    if(size != 1){

      MPI_Waitall(8, r, status);
      B[m+1][1:m] = recv_buffer_0[0:m];
      B[0][1:m] = recv_buffer_1[0:m];
      B[1:m][0] = recv_buffer_2[0:m];
      B[1:m][m+1] = recv_buffer_3[0:m];
      B[m+1][0] = recv_buffer_4[0];
      B[0][0] = recv_buffer_5[0];
      B[0][m+1] = recv_buffer_6[0];
      B[m+1][m+1] = recv_buffer_7[0];
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
    
    // update boundary points

    // i = 0
    if(x != 0)
      A[0][0:m] = 0.125f*(B[1][2:m] + B[1][0:m] + B[0][2:m] + B[0][0:m] + B[2][2:m] + B[2][0:m] + B[2][1:m] + B[0][1:m]);
    
    // i = m-1
    if(x != z - 1)
      A[m-1][0:m] = 0.125f*(B[m][2:m] + B[m][0:m] + B[m-1][2:m] + B[m-1][0:m] + B[m+1][2:m] + B[m+1][0:m] + B[m+1][1:m] + B[m-1][1:m]);

    // j = 0
    A[1:m-2][0] = 0.125f*(B[2:m][2] + B[2:m][0] + B[1:m][2] + B[1:m][0] + B[3:m][2] + B[3:m][0] + B[3:m][1] + B[1:m][1]);
    // j = m-1
    A[1:m-2][m-1] = 0.125f*(B[2:m][m+1] + B[2:m][m-1] + B[1:m][m+1] + B[1:m][m-1] + B[3:m][m+1] + B[3:m][m-1] + B[3:m][m] + B[1:m][m]);

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
    printf("sum after %d iterations = %10.10f m=%d\n", t, total_sum, m);
    
    FILE *fp;
    fp=fopen("timings.dat", "a+");
    fprintf(fp, "%d %d %f %10.16f\n", size, N, endtime-starttime, total_sum);
    fclose(fp);


  }

  MPI_Finalize();
  return 0;
}
