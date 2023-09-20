/* PURPOSE:
 * Compute the best column addition to a matrix to lower the singular value the most
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define BIG_NONSENSE 1e12

void write_field( char *filename, double *data, int N );
void read_field(  char *filename, double *data, int N );


double secant_method( double l, double u, int n, double *s, double *wi, double al, double tau, double threshold );







int main( int argc, char **argv ){

  if( argc != 9 ){
    printf("Error: wrong number of arguments.\n");
    return -1;
  }


  int n = atoi( argv[1] ); //Current size of library
  int p = atoi( argv[2] ); //number of available modifications
  
  double threshold = atof( argv[3] ); //number of available modifications

  double *s   = (double *) malloc( sizeof(double) * n );
  double *w   = (double *) malloc( sizeof(double) * n * p );
  double *al  = (double *) malloc( sizeof(double) * p );
  double *tau = (double *) malloc( sizeof(double) * p );

  double *s_min = (double*) malloc( sizeof(double) * p );

  //Read in all the data
  read_field( argv[4], s,   n   );
  read_field( argv[5], w,   n*p );
  read_field( argv[6], al,  p   );
  read_field( argv[7], tau, p   );

  //Identify the minimal singular value for each possible modification
  //This is done via secant method

  double l = 0;             //lower bound
  double u = s[n-1]*s[n-1]; //upper bound
  
  double best = BIG_NONSENSE; //smallest singular value so far!

  for( int i=0; i<p; i++){
    double *wi = &w[ i*n ];

    s_min[i] = secant_method( l, u, n, s, wi, al[i], tau[i], threshold );

    s_min[i] = sqrt( s_min[i] );
  }

  //Write out the new singular values
  write_field( argv[8], s_min, p );
}


void write_field( char *filename, double *data, int N ){
  FILE *stream = fopen(filename, "w"); //open the file
  fwrite( &N,   sizeof(int),    1, stream );
  fwrite( data, sizeof(double), N, stream );
  fclose(stream);
}

void read_field( char *filename, double *data, int N ){
  //Inverse of write field
  FILE *stream = fopen(filename, "r"); //open the file
  int N2;  
  fread( &N2,   sizeof(int),          1,    stream );
  //printf("%s %d %d\n", filename, N, N2 );
  assert( N == N2 );

  //Now read it
  fread( data, sizeof(double), N, stream );
  fclose(stream);
}


double objective( double x, int n, double *s, double *wi, double al, double tau ){
  /* INPUT:
   * x  - the value of the function you want to investigate
   * n  - see main
   * s  - singular values of the existing library matrix
   * wi - pointer to vector of interest: a column of w
   * al - value of 1/alpha^2
   * tau- see paper 
   *
   * OUTPUT:
   * value of secular equation with poles from 0 and smallest sigma removed
   */
 
  //From working MATLAB code:
  // f = @(k)  k*(s(n).^2 - k) + (k/alpha^2)*wbar(n).^2 + (k*(s(n).^2 - k))/alpha^2*sum( wbar(1:n-1).^2 ./ (s(1:n-1).^2 - k) )  - (s(n).^2 - k).*tau_sq/alpha^2;

 
  double val = 0; //return value
  
  double reg = (s[n-1]*s[n-1] - x); //regularization factor

  //First term
  val =  x*reg;

  //second term
  val += x*al*wi[n-1]*wi[n-1];

  //third term
  for(int i=0; i<n-1; i++){
    val += x*al*reg/(s[i]*s[i] - x)*wi[i]*wi[i];
  }

  //fourth term
  val -= reg*tau*al;

  return val;
}



double secant_method( double l, double u, int n, double *s, double *wi, double al, double tau, double threshold ){
  //Compute function value at bounds
  double fl = objective( l, n, s, wi, al, tau);
  double fu = objective( u, n, s, wi, al, tau);

  //Check if we can use secant method. If not, exit early!
  if( fl*fu > 0 ){
    return BIG_NONSENSE;
  }

  //Allocate memory for the function value of our guess-point
  double g; //g for guess
  double fg = BIG_NONSENSE; //Make up nonsense value (large so the loop starts)

  int maxit = 1024;
  int it = 0;
  while( fabs(fg) > threshold && it < maxit  ){
    it++;

    g = (u+l)/2; //bijection
    
    //g = u - l*(fu-fl)/(u-l); //Use the secant method to guess a root
    fg = objective( g, n, s, wi, al, tau); //Compute the function

    //printf("%d: g = %e\t fg = %e\t fu = %e\t fl = %e\n", it, g, fg, fu, fl);

    if( fg*fu > 0 ){
      u = g;
      fu=fg;
    }
    else{
      l = g;
      fl=fg;
    }
  }

  return g;
}
