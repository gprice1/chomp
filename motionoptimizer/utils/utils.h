/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>


#ifdef DEBUG
    #define debug_assert assert
    #define debug std::cout
    #define debug_status( tag, function, point ) std::cout << tag << " - " << function << " - " << point << std::endl
#else 
    #define debug_assert if(0) assert
    #define debug if(0) std::cout
    #define debug_status( tag, function, point ) if(0) std::cout << tag 
#endif

namespace chomp {

/////////////////////typedefs///////////////////////////////

typedef Eigen::MatrixXd MatX;
typedef Eigen::Map<MatX>  MatMap;
typedef const Eigen::Map< const MatX>  ConstMatMap;

//this matrix map is used for subsampling a matrix. 
typedef Eigen::Stride<Eigen::Dynamic,2> SubMatMapStride;
typedef Eigen::Map<MatX, 0, SubMatMapStride > SubMatMap; 
typedef const Eigen::Map<const MatX, 0, SubMatMapStride > ConstSubMatMap;

typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> DynamicStride;
typedef Eigen::Map<MatX, 0, DynamicStride > DynamicMatMap; 
typedef const Eigen::Map<const MatX, 0, DynamicStride > ConstDynamicMatMap;

typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, 
                       Eigen::RowMajor> MatXR;
typedef Eigen::Map<MatXR> MatMapR;
typedef const Eigen::Map<const MatXR> ConstMatMapR;
typedef Eigen::Map<MatXR, 0, Eigen::OuterStride<> > SubMatMapR; 

typedef Eigen::Block<DynamicMatMap, 1, Eigen::Dynamic> Row;
typedef const Eigen::Block<const DynamicMatMap, 1, Eigen::Dynamic> ConstRow;
typedef Eigen::Block<DynamicMatMap, Eigen::Dynamic, 1, Eigen::Dynamic> Col;
typedef const Eigen::Block< const DynamicMatMap, Eigen::Dynamic,
                            1, Eigen::Dynamic> ConstCol;


/////////////////forward class declarations///////////////
class ChompGradient;
class ChompOptimizerBase;
class ConstraintFactory;
class Constraint;
class HMC;
class Chomp;
class Trajectory;

class OptimizerBase;
class ChompOptimizerBase;
class ChompOptimizer;
class ChompLocalOptimizer;

enum ChompEventType { 
    CHOMP_INIT,
    CHOMP_GLOBAL_ITER,
    CHOMP_LOCAL_ITER,
    NLOPT_ITER,
    CHOMP_FINISH,
    CHOMP_TIMEOUT,
    CHOMP_GOALSET_ITER,
    CHOMP_COVARIANT_ITER

};



/////////////////Utility types /////////////////////////
// These classes interface with chomp, to do something useful
class ChompObserver {
  public:
    virtual ~ChompObserver(){}
    virtual int notify(const OptimizerBase& c, 
                       ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
        return 0;
    }
};

class DebugChompObserver: public ChompObserver {
  public:
    virtual ~DebugChompObserver(){}
    virtual int notify(const OptimizerBase& c, 
                       ChompEventType e,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
            
         std::string event_string;

         switch (e) {
            case CHOMP_INIT:
                 event_string = "CHOMP_INIT"; break;
            case CHOMP_GLOBAL_ITER:
                 event_string = "CHOMP_GLOBAL_ITER"; break;
            case CHOMP_LOCAL_ITER:
                 event_string = "CHOMP_LOCAL_ITER"; break;
            case NLOPT_ITER:
                 event_string = "NLOPT_ITER"; break;
            case CHOMP_FINISH:
                 event_string = "CHOMP_FINISH"; break;
            case CHOMP_TIMEOUT:
                 event_string = "CHOMP_TIMEOUT"; break;
            case CHOMP_GOALSET_ITER: 
                 event_string = "CHOMP_GOALSET_ITER"; break;
            case CHOMP_COVARIANT_ITER:
                 event_string = "CHOMP_COVARIANT_ITER"; break;
            default:
                 event_string = "[INVALID]"; break;
         }

         std::cout << "chomp debug: "
              << "event=" << event_string << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";
        if ( e != CHOMP_INIT && e != CHOMP_FINISH && 
            (std::isnan(curObjective) || std::isinf(curObjective) ||
             std::isnan(lastObjective) || std::isinf(lastObjective)) )
        {
            return 1;
        }
        return 0;
    }
};

enum ChompObjectiveType {
    MINIMIZE_VELOCITY = 0,
    MINIMIZE_ACCELERATION = 1,
};

///////////////////////////////////////////////////////////
/////////////////////Inline Functions//////////////////////

//Copy a MatX into a std::vector of doubles
inline void matToVec( const MatX & mat, std::vector<double> & vec ){
    const double * data = mat.data();
    vec = std::vector<double>(data, data + mat.size());
}

//Copy a vector of doubles into a column-major Matrix.
inline void vecToMat( const std::vector<double> & vec, MatX & mat )
{
    if ( mat.size() == int(vec.size()) ){
        mat = ConstMatMap( vec.data(), mat.rows(), mat.cols() );
    }
    else {
        mat = ConstMatMap( vec.data(), 1, vec.size());
    }
}

//Copy a vector of doubles into a MatX.
inline void vecToMat( const std::vector<double> & vec,
                      int N, int M, MatX & mat)
{
    assert( int( vec.size() ) == N * M );
    mat = ConstMatMap( vec.data(), N, M );
}

//this is a function to compute the dot product of matrices
template <class Derived1, class Derived2>
static inline double mydot(const Eigen::MatrixBase<Derived1>& a,
                           const Eigen::MatrixBase<Derived2>& b)
{
  return a.cwiseProduct(b).sum();
}


template <class Derived>
inline MatX getPos(const Eigen::MatrixBase<Derived>& x, double h){
    double fac = 1;
    double hn = 1;

    MatX rval = x.row(0);
    for (int i=1; i<x.rows(); ++i) {
      fac *= i;
      hn *= h;
      rval += (hn / fac) * x.row(i);
    }

    return rval;
}



///////////////////////////////////////////////////////////////
//////////////////////////Linear Algebra routines//////////////

//coeffs is a vector of coefficients that correspond to the
//    the values in the chomp A matrix,
//      e.g. [-1,2] or [1,-4,6]
template <class Derived1, class Derived2, class Derived3>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs,                                  const Eigen::MatrixBase<Derived2>& x,
             const Eigen::MatrixBase<Derived3>& Ax_const );

//Creates a Lower triangular matrix that corresponds to the
//  cholesky decomposition of the A matrix.
//  L is stored as a skyline matrix for size and speed efficiency.
//coeffs is a vector of coefficients that correspond to the
//    the values in the chomp A matrix,
//      e.g. [-1,2] or [1,-4,6]
template <class Derived1, class Derived2>
void skylineChol(int n, const Eigen::MatrixBase<Derived1>& coeffs,
                 Eigen::PlainObjectBase<Derived2>& L);


//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived1, class Derived2>
void skylineCholMultiply(const Eigen::MatrixBase<Derived1>& L,
                         const Eigen::MatrixBase<Derived2>& x_const);



//Perform an inverse transpose multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L^-T * x_const;
//This is used by the HMC class to generate random smooth momenta.
template <class Derived1, class Derived2>
void skylineCholMultiplyInverseTranspose(
                               const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const);

//Perform an inverse multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L^-1 * x_const;
template <class Derived1, class Derived2>
void skylineCholMultiplyInverse(
                               const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const);

//Solve for the inverse multiplication of the A matrix using
//  the cholesky decomposition.  
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L^-T * L^-1 * x_const, which is equivalent to:
//      x_const = A^-1 * x_const
template <class Derived1, class Derived2>
void skylineCholSolve(const Eigen::MatrixBase<Derived1>& L,
                      const Eigen::MatrixBase<Derived2>& x_const);

//Perform the above cholesky decomposition on a matrix that has
//  the height several times the size of L.
//  This is a utility function to perform cholesky decomposition on
//      a large amount of data.
template <class Derived1, class Derived2>
void skylineCholSolveMulti(
                const Eigen::MatrixBase<Derived1>& L, 
                const Eigen::MatrixBase<Derived2>& xx_const);

template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived2>& x1,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt);


/////////////////////////////////////////////////////////////////////
//////////////////////Matrix Methods for Goalset-CHOMP///////////////
  
//this is a diag mul for goal set chomp
template <class Derived1, class Derived2, class Derived3, class Derived4>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs,
             const Eigen::MatrixBase<Derived2>& gs_coeffs,
             const Eigen::MatrixBase<Derived3>& x,
             const Eigen::MatrixBase<Derived4>& Ax_const );


//this is a skyline chol for goal set chomp
template <class Derived1, class Derived2, class Derived3>
void skylineChol(int n,
                 const Eigen::MatrixBase<Derived1>& coeffs,
                 const Eigen::MatrixBase<Derived2>& gs_coeffs,
                 Eigen::PlainObjectBase<Derived3>& L);



//This version of createBMatrix is used for goal set chomp.
template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt);

////////////////////////////////////////////////////////////////////
/////////////Linear Algebra Implementations ////////////////////////

template <class Derived1, class Derived2, class Derived3>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs, // e.g. [1, -4, 6]
             const Eigen::MatrixBase<Derived2>& x,
             const Eigen::MatrixBase<Derived3>& Ax_const ) {

    assert( Ax_const.rows() == x.rows() && Ax_const.cols() == x.cols() );

    Eigen::MatrixBase<Derived3>& Ax = 
      const_cast<Eigen::MatrixBase<Derived3>&>(Ax_const);
    
    //in the case that we are doing goal chomp, n != x.rows(),
    //  so if goal chomp is being done, then, the correct value of n,
    //  should be passed in.
    int n = x.rows();

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    for (int i=0; i<x.rows(); ++i) {

      int j0 = std::max(i-o, int(0));
      int j1 = std::min(i+nc, n );

      Ax.row(i) = x.row(j0) * coeffs(j0-i+o);

      for (int j=j0+1; j<=i; ++j) {
        Ax.row(i) += x.row(j) * coeffs(j-i+o);
      }
    
      for (int j=i+1; j<j1; ++j) {
        Ax.row(i) += x.row(j) * coeffs(i-j+o);
      }

    }
}

template <class Derived1, class Derived2>
void skylineChol(int n,
                 // e.g. [-1,2] or [1,-4,6]
                 const Eigen::MatrixBase<Derived1>& coeffs,
                 Eigen::PlainObjectBase<Derived2>& L) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    L.resize(n, nc);
  
    for (int j=0; j<n; ++j) {
      
      //i1 is the current row, forwarded by the amount of coeffs.
      int i1 = std::min(j+nc, n);
      
      for (int i=j; i<i1; ++i) {

        double sum = 0;

        int k0 = std::max(0,i-o);

        for (int k=k0; k<j; ++k) {
          sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
        }

        if (i == j) {
          L(j,o) = sqrt(coeffs(o) - sum);
        } else {
          L(i,j-i+o) = (coeffs(j-i+o) - sum) / L(j,o);
        }
      }
    }

}


//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived1, class Derived2>
void skylineCholMultiply(const Eigen::MatrixBase<Derived1>& L,
                         const Eigen::MatrixBase<Derived2>& x_const)
{
    const int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
        const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    const int nc = L.cols();
    const int o = nc-1;
    
    // i is the column index 
    // j is the 
    for (int i=n-1; i>=0; --i) {
        
        x.row(i) *= L(i, o);
        
        const int j1 = std::max( 0, i - nc);
        for (int j = i-1; j > j1; j--) {
            x.row(i) += L(i, o+j-i) * x.row( j ); //
        }
    }
}

//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived1, class Derived2>
void skylineCholMultiplyTranspose(const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const)
{
    const int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
        const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    const int nc = L.cols();
    const int final_column_index = nc-1;
     
    // i is the column index 
    // j indexes into the correct row of xi
    for (int i=0 ; i < n ; i++ )
    {
        x.row(i) *= L(i, final_column_index);

        const int j1 = std::min( n, i+nc);
        for (int j = i+1; j < j1; j++) {
            const int L_col_index = final_column_index - j + i ;
            x.row(i) += L(j, L_col_index ) * x.row( j ); //
        }
    }
}

//////////////////////////////////////////////////////////////////////
//This is used by the HMC class to generate random smooth momenta.
template <class Derived1, class Derived2>
void skylineCholMultiplyInverseTranspose(
                               const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const)
{

    const int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    const int nc = L.cols();
    const int o = nc-1;

    for (int i=n-1; i>=0; --i) {
      const int j1 = std::min(i+nc, n);
      for (int j=i+1; j<j1; ++j) {
        x.row(i) -= L(j, i-j+o) * x.row(j); // here j > i so col < row
      }
      x.row(i) /= L(i,o);
    }
}

template <class Derived1, class Derived2>
void skylineCholMultiplyInverse(const Eigen::MatrixBase<Derived1>& L,
                                const Eigen::MatrixBase<Derived2>& x_const)
{
    int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    int nc = L.cols();
    int o = nc-1;

    for (int i=0; i<n; ++i) {
      int j0 = std::max(0, i-o);
      for (int j=j0; j<i; ++j) {
        x.row(i) -= L(i,j-i+o)*x.row(j); // here j < i so col < row
      }
      x.row(i) /= L(i,o);
    }
}


template <class Derived1, class Derived2>
void skylineCholSolve(const Eigen::MatrixBase<Derived1>& L,
                      const Eigen::MatrixBase<Derived2>& x_const)
{

    int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    int nc = L.cols();
    int o = nc-1;

    for (int i=0; i<n; ++i) {
      int j0 = std::max(0, i-o);
      for (int j=j0; j<i; ++j) {
        x.row(i) -= L(i,j-i+o)*x.row(j); // here j < i so col < row
      }
      x.row(i) /= L(i,o);
    }

    for (int i=n-1; i>=0; --i) {
      int j1 = std::min(i+nc, n);
      for (int j=i+1; j<j1; ++j) {
        x.row(i) -= L(j,i-j+o) * x.row(j); // here j > i so col < row
      }
      x.row(i) /= L(i,o);
    }

}


template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived2>& x1,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt)
{

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    assert(b_const.rows() == n);
    assert(b_const.cols() == x0.cols());
    assert(b_const.cols() == x1.cols());

    Eigen::MatrixBase<Derived3>& b = 
      const_cast<Eigen::MatrixBase<Derived3>&>(b_const);

    b.setZero();

    double c = 0;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      int i1 = n-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        int t1 = -t0;
        
        b.row(i0) += coeffs(j)*getPos(x0, t0*dt);
        b.row(i1) += coeffs(j)*getPos(x1, t1*dt);
      }
      c += mydot(b.row(i0), b.row(i0));

      if (i0 != i1) { c += mydot(b.row(i1), b.row(i1)); }
    }

    return 0.5*c;
}
  
  
//this is a diag mul for goal set chomp
template <class Derived1, class Derived2, class Derived3, class Derived4>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs,
             const Eigen::MatrixBase<Derived2>& gs_coeffs,
             const Eigen::MatrixBase<Derived3>& x,
             const Eigen::MatrixBase<Derived4>& Ax_const ){

    assert( Ax_const.rows() == x.rows() && Ax_const.cols() == x.cols() );

    Eigen::MatrixBase<Derived4>& Ax = 
      const_cast<Eigen::MatrixBase<Derived4>&>(Ax_const);
    
    //in the case that we are doing goal chomp, n != x.rows(),
    //  so if goal chomp is being done, then, the correct value of n,
    //  should be passed in.
    int n = x.rows();

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    const int start_gs = x.rows() - gs_coeffs.rows();

    for (int i=0; i<x.rows(); ++i) {

      int j0 = std::max(i-o, int(0));
      int j1 = std::min(i+nc, n );
      
      Ax.row(i) = x.row(j0) * coeffs(j0-i+o);

      for (int j=j0+1; j<=i; ++j) {
        double coeff;
        if ( i >= start_gs && j >= start_gs ){ 
            coeff = gs_coeffs( i - start_gs, j - start_gs );
        }else { coeff = coeffs(j-i+o); }

        Ax.row(i) += x.row(j) * coeff;
      }
    
      for (int j=i+1; j<j1; ++j) {
        double coeff;
        if ( i >= start_gs && j >= start_gs ){ 
            coeff = gs_coeffs( i - start_gs, j - start_gs );
        }else { coeff = coeffs(i-j+o); }

        Ax.row(i) += x.row(j) * coeff;
      }
    }
}

//this is a skyline chol for goal set chomp
template <class Derived1, class Derived2, class Derived3>
void skylineChol(int n,
                 const Eigen::MatrixBase<Derived1>& coeffs,
                 const Eigen::MatrixBase<Derived2>& gs_coeffs,
                 Eigen::PlainObjectBase<Derived3>& L) {

    int nc = coeffs.size();
    int o = nc-1;

    const int start_gs = n - gs_coeffs.rows();

    L.resize(n, nc);
  
    for (int j=0; j<n; ++j) {
      
      //i1 is the current row, forwarded by the amount of coeffs.
      int i1 = std::min(j+nc, n);
      
      for (int i=j; i<i1; ++i) {

        double sum = 0;

        int k0 = std::max(0,i-o);

        for (int k=k0; k<j; ++k) {
          sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
        }
        
        if (i == j) {
          double coeff;
          if ( i >= start_gs ){
            coeff = gs_coeffs( i - start_gs , i - start_gs);
          } else { coeff = coeffs(o); }

          L(j,o) = sqrt(coeff - sum);

        } else {
          double coeff;
          if ( i >= start_gs && j >= start_gs){
            coeff = gs_coeffs( i - start_gs, j - start_gs);
          } else { coeff = coeffs(j-i+o); }

          L(i,j-i+o) = (coeff - sum) / L(j,o);

        }
      }
    }
}



//This version of createBMatrix is used for goal set chomp.
template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    assert(b_const.rows() == n);
    assert(b_const.cols() == x0.cols());

    Eigen::MatrixBase<Derived3>& b = 
      const_cast<Eigen::MatrixBase<Derived3>&>(b_const);

    b.setZero();

    double c = 0;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        b.row(i0) += coeffs(j)*getPos(x0, t0*dt);
      }
      c += mydot(b.row(i0), b.row(i0));
    }

    return 0.5*c;
}

}//namespace

#endif
