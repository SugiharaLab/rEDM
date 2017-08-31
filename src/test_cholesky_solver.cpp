#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

MatrixXd stable_cholesky_solver( LDLT<MatrixXd> ldltDecomp,
				 MatrixXd A,
				 bool transpose = false )
{

  // Preparations:

  // For some reason if I sub it below I get error
  MatrixXd L = ldltDecomp.matrixL();

  // Number of rows is all that matters, regardless if rhs is a
  // matrix or a vector
  int k = A.rows(); 

  // Manually inverting D. This procedure has the advantage that
  // D^{-1/2} can also be applied to matrices.
  VectorXd diag;
  diag.resize(k);
  for( int i = 0 ; i < k ; ++i ) 
    diag(i) = 1. / sqrt( ldltDecomp.vectorD()(i) ) ; // Manual inversion
  DiagonalMatrix<double, Dynamic > sqrtInvD = diag.asDiagonal();

  // The permutation "matrix" P
  Transpositions<Dynamic> P = ldltDecomp.transpositionsP(); 

  // Holds all the computations
  MatrixXd x;
  
  // Now the actual computation

  if( !transpose ) {
      
      // x = PA
      x = P * A;
      
      // x = L^{-1}PA
      x = L.triangularView<Lower>().solve(x);
      
      // x = D^{-1/2}L^{-1}PA
      x = sqrtInvD * x;
      
    } else {

    // x = D^{-1/2}A
    x = sqrtInvD * A;

    // x = L^{-t}D^{-1/2}A
    x = L.triangularView<Lower>().transpose().solve(x);

    // x = P^tL^{-t}D^{-1/2}A
    x = P.transpose() * x; 
  }
  return x;
  
}



int main()
{

  int k = 3; // Dimensionality

  // Define, declare and enter the problem's data
  MatrixXd A;
  A.resize(k, k);
  MatrixXd b;
  b.resize(k, 2 );

  A <<
    13, 5, 7 ,
    5 , 9, 3 ,
    7 , 3, 11;
  b <<
    3, 3, 4,
    1,-2, 9;
 
  cout << "Here is the " << A.rows() << " by " << A.cols() << " matrix A:\n" << A << endl;
  cout << "Here is the " << b.rows() << " by " << b.cols() << " matrix b:\n" << b << endl;
  cout << "Let's solve Ax = b using different methods.\n" <<endl;

  // Two placeholders that will be used throughout
  MatrixXd L;
  MatrixXd x;

  // ldlt()
  cout << "\n\nUsing the stable Cholesky decompostion ldlt()" << endl;

  // The object that actually holds the entire decomposition
  LDLT<MatrixXd> ldltDecomp = A.ldlt();

  // Direct solution, using Eigen's routines
  x = ldltDecomp.solve(b);
  cout << "Direct x =\n" << x << endl;
  cout << "Direct b =\n" << A*x << endl;

  // Manual solution - implementing the process Eigen is taking, step
  // by step (in the function defined above). 
  x = stable_cholesky_solver( ldltDecomp, b );
  x = stable_cholesky_solver( ldltDecomp, x, true );
 
  cout << "Manual x =\n" << x << endl;
  cout << "Manual b =\n" << A*x << endl;


  // llt()
  cout << "\n\nUsing the less stable, but faster Cholesky decomposition " <<
    "without pivoting llt()" << endl;

  // Here A.llt() is the object that actually holds the decomposition
  // (like ldltDecomp before) but we only need the matrix L.
  L = A.llt().matrixL();
  
  x = L.triangularView<Lower>().solve(b);
  x = L.triangularView<Lower>().transpose().solve(x);
  cout << "Manual x =\n" << x << endl;
  cout << "Manual b =\n" << A*x << endl;
   
}

// https://stackoverflow.com/questions/24442850/solving-a-sparse-upper-triangular-system-in-eigen

// https://stackoverflow.com/questions/21645604/solve-for-inverse-square-root

// https://stackoverflow.com/questions/20181940/most-efficient-way-to-solve-a-system-of-linear-equations
