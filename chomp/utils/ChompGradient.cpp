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

#include <float.h>
#include <cmath>
#include <iomanip>

#include "ChompGradient.h"

#define DEBUG_PRINTING 1
#if DEBUG_PRINTING 
    #define debug std::cout
    #define debug_assert assert
#else
    #define debug if (0) std::cout
    #define debug_assert if (0) assert
#endif

namespace chomp {

ChompGradientHelper::~ChompGradientHelper() {}

ChompCollisionHelper::ChompCollisionHelper(size_t nc,
                                           size_t nw,
                                           size_t nb):
    ncspace(nc), nwkspace(nw), nbodies(nb) {}

ChompCollisionHelper::~ChompCollisionHelper() {}

ChompCollGradHelper::ChompCollGradHelper(ChompCollisionHelper* h,
                                           double g):
    chelper(h), gamma(g)
{
    dx_dq = MatX(h->nwkspace, h->ncspace);
    cgrad = MatX(h->nwkspace, 1);
}

ChompCollGradHelper::~ChompCollGradHelper() {}


template< class Derived1, class Derived2, class Derived3>
inline double ChompCollGradHelper::computeGradient(
                    const Eigen::MatrixBase<Derived1> & xi,
                    const Eigen::MatrixBase<Derived2> & pinit,
                    const Eigen::MatrixBase<Derived2> & pgoal,
                    double dt,
                    const Eigen::MatrixBase<Derived3> & g_const)
{
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived3>& g = 
        const_cast<Eigen::MatrixBase<Derived3>&>(g_const);

    q1 = getTickBorderRepeat(-1, xi, pinit, pgoal, dt).transpose();
    q2 = getTickBorderRepeat(0,  xi, pinit, pgoal, dt).transpose();
    
    const double inv_dt = 1/dt;
    const double inv_dt_squared = inv_dt * inv_dt;

    double total = 0.0;

    for (int t=0; t < xi.rows() ; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = getTickBorderRepeat(t+1, xi, pinit, pgoal, dt).transpose();

      cspace_vel = 0.5 * (q2 - q0) * inv_dt;        
      cspace_accel = (q0 - 2.0*q1 + q2) * inv_dt_squared;

      for (size_t u=0; u < chelper->nbodies; ++u) {

        float cost = chelper->getCost(q1, u, dx_dq, cgrad);
        if (cost > 0.0) {

          wkspace_vel = dx_dq * cspace_vel;

          //this prevents nans from propagating. Several lines below, 
          //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
          if (wkspace_vel.isZero()){ continue; }

          wkspace_accel = dx_dq * cspace_accel;
          
          float wv_norm = wkspace_vel.norm();
          wkspace_vel /= wv_norm;

          // add to total
          double scl = wv_norm * gamma * dt;

          total += cost * scl;
          
          P = MatX::Identity(chelper->nwkspace, chelper->nwkspace)
              - (wkspace_vel * wkspace_vel.transpose());

          K = (P * wkspace_accel) / (wv_norm * wv_norm);
         

          //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
          g.row(t) += (scl * (dx_dq.transpose() *
                      (P * cgrad - cost * K)).transpose());
         

        }
      }
    }
    return total;

}
double ChompCollGradHelper::addToGradient(const MatX& xi,
                                          const MatX& pinit,
                                          const MatX& pgoal,
                                          double dt,
                                          MatX& g)
{
    return computeGradient( xi, pinit, pgoal, dt, g );
}

double ChompCollGradHelper::addToGradient(ConstMatMap& xi,
                                          const MatX& pinit,
                                          const MatX& pgoal,
                                          double dt,
                                          MatMap& g)
{
    return computeGradient( xi, pinit, pgoal, dt, g );
}

ChompGradient::ChompGradient(const MatX& pinit, 
                             const MatX& pgoal, 
                             ChompObjectiveType objective_type,
                             double total_time) :
    ghelper(NULL),
    objective_type( objective_type ),
    q0( pinit ), q1( pgoal ),
    t_total( total_time )
{   
    M = q0.size();

    if (objective_type == MINIMIZE_VELOCITY) {
        coeffs.resize(1,2);
        coeffs_sub.resize(1,1);
        coeffs_goalset.resize(1,1);

        coeffs << -1, 2;
        coeffs_sub << 2;
        coeffs_goalset << 1;
    } else {
        coeffs.resize(1,3);
        coeffs_sub.resize(1,2);
        coeffs_goalset.resize(2,2);
        
        coeffs << 1, -4, 6;
        coeffs_sub << 1, 6;
        coeffs_goalset << 6, -3,
                         -3,  2 ;
    }

}

void ChompGradient::prepareRun(int N,
                               bool use_goalset,
                               bool subsample)
{
    iteration = 0;
    
    this->N = N;
    this->use_goalset = use_goalset;

    //resize the g, b, and ax matrices.
    g.resize(N,M);
    Ax.resize(N,M);
    b.resize(N,M);
    
    //set b to zero to prepare for creating the b matrix
    b.setZero();

    //get the b matrix, and get its contribution to the
    //  objective function
    if (use_goalset){
        dt = t_total/N; 
        skylineChol(N, coeffs, coeffs_goalset, L);
        c = createBMatrix(N, coeffs, q0, b, dt);

    } else{
        dt = t_total/(N+1);
        skylineChol(N, coeffs, L);
        c = createBMatrix(N, coeffs, q0, q1, b, dt);
    }
    
    inv_dt = 1/dt;

    if ( objective_type == MINIMIZE_VELOCITY ){ fscl = inv_dt * inv_dt;}
    else { fscl = inv_dt * inv_dt * inv_dt;  }


    if (subsample) {
        int N_sub = (N+1)/2;
        skylineChol(N_sub, coeffs_sub, L_sub); 
    }
}


MatX& ChompGradient::getInvAMatrix( bool subsample){
    return (subsample ? L_sub : L );
}

MatX& ChompGradient::getCollisionGradient( const MatX & xi )
{
    g.setZero();
    computeCollisionGradient( xi, g );
    return g;

}

MatX& ChompGradient::getGradient( const MatX & xi )
{
    
    computeSmoothnessGradient( xi, g );
    computeCollisionGradient( xi, g );

    return g;

}

MatX& ChompGradient::getSmoothnessGradient( const MatX & xi )
{
    
    computeSmoothnessGradient( xi, g );
    return g;
}

MatX& ChompGradient::getSubsampledGradient(int N_sub)
{   
    g_sub.resize( N_sub, M );

    for ( int i = 0; i < g_sub.rows(); i ++ ){
        g_sub.row( i ) = g.row( i * 2 );
    }

    return g_sub;
}

double ChompGradient::getGradient( unsigned n_by_m,
                                   const double * xi,
                                   double * grad)
{
    iteration ++;

    ConstMatMap xi_mat( xi, N, M );
    
    if ( grad != NULL ){

        assert( unsigned(N*M) == n_by_m );
        MatMap g_mat( grad, N, M);
        
        g_mat.setZero();

        computeSmoothnessGradient( xi_mat, g_mat );
        computeCollisionGradient(  xi_mat, g_mat );
        
        //skylineCholSolve( L, g_mat );

    }else{
        computeSmoothnessGradient( xi_mat, g );
        computeCollisionGradient(  xi_mat, g );
    }


    const double cost = evaluateObjective( xi_mat );
    
    debug << "Iteration[" << iteration << "] -- Cost: " << cost << "\n";
    return cost;
}


template <class Derived1, class Derived2>
inline void ChompGradient::computeSmoothnessGradient(
                    const Eigen::MatrixBase<Derived1> & xi,
                    const Eigen::MatrixBase<Derived2> & g_const)
{
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived2>& grad = 
        const_cast<Eigen::MatrixBase<Derived2>&>(g_const);

    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    if( use_goalset ){ diagMul(coeffs, coeffs_goalset, xi, Ax); }
    else { diagMul(coeffs, xi, Ax); }
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    grad = Ax + b;
}

void ChompGradient::computeCollisionGradient(ConstMatMap & xi,
                                             MatMap & grad)
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(xi, q0, q1, dt, grad);
    } else {
        fextra = 0;
    }
}

void ChompGradient::computeCollisionGradient(const MatX & xi,
                                             MatX & grad)
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(xi, q0, q1, dt, grad);
    } else {
        fextra = 0;
    }
}

}// namespace

