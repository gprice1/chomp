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

#ifndef _CHOMPUTIL_H_
#define _CHOMPUTIL_H_

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

namespace chomp {

/////////////////////typedefs///////////////////////////////

typedef Eigen::MatrixXd MatX;
typedef Eigen::Map<MatX>  MatMap;
//this matrix map is used for subsampling a matrix. 
typedef Eigen::Map<MatX, 0, Eigen::InnerStride<2> > SubMatMap; 

typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, 
                       Eigen::RowMajor> MatXR;
typedef Eigen::Map<MatXR> MatMapR;
typedef Eigen::Map<MatXR, 0, Eigen::OuterStride<> > SubMatMapR; 



/////////////////forward class declarations///////////////

class ConstraintFactory;
class Constraint;
class HMC;
class Chomp;

enum ChompEventType { 
    CHOMP_INIT,
    CHOMP_GLOBAL_ITER,
    CHOMP_LOCAL_ITER,
    CHOMP_FINISH,
    CHOMP_TIMEOUT,
    CHOMP_GOALSET_ITER,
};

//prints out a string that corresponds to the variable name
//  of the event type.
const char* eventTypeString(int eventtype);


/////////////////Utility types /////////////////////////
// These classes interface with chomp, to do something useful
class ChompObserver {
  public:
    virtual ~ChompObserver();
    virtual int notify(const Chomp& c, 
                       ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation);
};

class DebugChompObserver: public ChompObserver {
  public:
    virtual ~DebugChompObserver();
    virtual int notify(const Chomp& c, 
                       ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation);
};

class ChompGradientHelper {
  public:
    virtual ~ChompGradientHelper();
    virtual double addToGradient(const Chomp& c, MatX& g) =0;
  };

class ChompCollisionHelper {
  public:

    // nbodies = number of bodies considered in gradient term
    // nwkspace = workspace dimension (2 or 3 probably)
    // ncspace = config. space dimension

    ChompCollisionHelper(size_t ncspace, size_t nwkspace, size_t nbodies);

    virtual ~ChompCollisionHelper();

    // return the cost for a given configuration/body, along with jacobians
    virtual double getCost(const MatX& q,         // configuration
                           size_t body_index,     // which body
                           MatX& dx_dq, // Jacobian of workspace
                                        // pos (nwkspace-by-ncspace)
                           MatX& cgrad) =0; // gradient (Jacobian transpose)
                                            // of cost wrt workspace pos
                                            //  (ncspace-by-1)
    
    
    size_t ncspace;
    size_t nwkspace;
    size_t nbodies;

};

class ChompCollGradHelper: public ChompGradientHelper {
  public:

    ChompCollGradHelper(ChompCollisionHelper* h, double gamma=0.1);
    virtual ~ChompCollGradHelper();

    virtual double addToGradient(const Chomp& c, MatX& g);

    ChompCollisionHelper* chelper;
    double gamma;

    MatX dx_dq;
    MatX cgrad;

    MatX q0, q1, q2;
    MatX cspace_vel, cspace_accel;
    MatX wkspace_vel, wkspace_accel;
    MatX P;
    MatX K; 
    
};

enum ChompObjectiveType {
    MINIMIZE_VELOCITY = 0,
    MINIMIZE_ACCELERATION = 1,
};

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

///////////////////////////////////////////////////////////
/////////////////////Inline Functions//////////////////////

//this is a function to compute the dot product of matrices
template <class Derived1, class Derived2>
static inline double mydot(const Eigen::MatrixBase<Derived1>& a,
                           const Eigen::MatrixBase<Derived2>& b) {

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

#endif
