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


template< class Derived >
inline double ChompCollGradHelper::computeGradient(
                    const Trajectory & traj,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
    q1 = traj.getTick( -1 ).transpose();
    q2 = traj.getTick( 0  ).transpose();
    
    const double inv_dt = 1/traj.getDt();
    const double inv_dt_squared = inv_dt * inv_dt;

    double total = 0.0;

    for (int t=0; t < traj.rows() ; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = traj.getTick( t+1 ).transpose();

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
          double scl = wv_norm * gamma * traj.getDt();

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
double ChompCollGradHelper::addToGradient( const Trajectory & traj,
                                          MatX& g)
{
    return computeGradient( traj, g );
}

double ChompCollGradHelper::addToGradient( const Trajectory & traj,
                                          MatMap& g)
{
    return computeGradient( traj, g );
}

ChompGradient::ChompGradient(
            Trajectory & traj,
            ChompObjectiveType objective_type) :
    trajectory( traj ),
    ghelper(NULL),
    objective_type( objective_type )

{   

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

const char* ChompGradient::TAG = "ChompGradient";

void ChompGradient::prepareRun(const Trajectory & traj,
                               bool use_goalset)
{
    debug_status( TAG, "prepareRun", "start" );

    iteration = 0;
    
    this->use_goalset = use_goalset;

    //resize the g, b, and ax matrices.
    const int N = traj.isSubsampled() ? traj.sampledN() :
                                        traj.N();
    const int M = traj.M();

    g.resize(N,M);
    Ax.resize(N,M);
    b.resize(N,M);
    
    //set b to zero to prepare for creating the b matrix
    
    const double dt = traj.getDt();
    //get the b matrix, and get its contribution to the
    //  objective function.
    //  Also, compute the L matrix (lower triangluar cholesky 
    //  decomposition). 
    if (use_goalset){
        skylineChol( traj.N(), coeffs, coeffs_goalset, L);
        c = createBMatrix( N, coeffs, traj.getQ0(), b, dt);
    } else{
        const MatX & current_coeffs = (traj.isSubsampled() ?
                                       coeffs_sub :
                                       coeffs ); 
        skylineChol( traj.N(), current_coeffs, L); 
        c = createBMatrix(N, coeffs, traj.getQ0(), traj.getQ1(), b, dt);
    }
    
    double inv_dt = 1/dt;

    if ( objective_type == MINIMIZE_VELOCITY ){ fscl = inv_dt * inv_dt;}
    else { fscl = inv_dt * inv_dt * inv_dt;  }

    debug_status( TAG, "prepareRun", "start" );
}

MatX& ChompGradient::getInvAMatrix(){ return L; }

MatX& ChompGradient::getCollisionGradient( const Trajectory & traj )
{
    g.setZero();
    computeCollisionGradient( traj , g );
    return g;
}

MatX& ChompGradient::getGradient( const Trajectory & traj )
{
    debug_status( TAG, "getGradient", "start" );
    
    computeSmoothnessGradient( traj, g );
    computeCollisionGradient( traj, g );
    
    if ( traj.isSubsampled() ){ return getSubsampledGradient(traj.N());}
    debug_status( TAG, "getGradient", "end" );
    
    return g;

}

MatX& ChompGradient::getSmoothnessGradient( const Trajectory & traj )
{
    
    computeSmoothnessGradient( traj, g );
    return g;
}

MatX& ChompGradient::getSubsampledGradient(int N_sub)
{   
    debug_status( TAG, "getSubsampledGradient", "start" );
    
    g_sub.resize( N_sub, trajectory.M() );

    for ( int i = 0; i < g_sub.rows(); i ++ ){
        g_sub.row( i ) = g.row( i * 2 );
    }

    debug_status( TAG, "getSubsampledGradient", "end" );
    
    return g_sub;
}

double ChompGradient::getGradient( unsigned n_by_m,
                                   const double * xi,
                                   double * grad)
{
    iteration ++;
    trajectory.setData( xi );
    
    const int N = trajectory.N();
    const int M = trajectory.M();

    if ( grad != NULL ){

        assert( unsigned(N*M) == n_by_m );
        MatMap g_mat( grad, N, M);
        
        g_mat.setZero();

        computeSmoothnessGradient( trajectory, g_mat );
        computeCollisionGradient(  trajectory, g_mat );
        
        //skylineCholSolve( L, g_mat );

    }else{
        computeSmoothnessGradient( trajectory, g );
        computeCollisionGradient(  trajectory, g );
    }


    const double cost = evaluateObjective( trajectory );
    
    debug << "Iteration[" << iteration << "] -- Cost: " << cost << "\n";
    return cost;
}

double ChompGradient::evaluateObjective( Trajectory & traj ) const
{
    if ( traj.isSubsampled() ){
        return 0.5 * mydot(traj.getSampledXi(), Ax) + 
               mydot(traj.getSampledXi(), b) +
               c + fextra;
    }
    
    return 0.5 * mydot(traj.getXi(), Ax) + 
           mydot(traj.getXi(), b) +
           c + fextra;
}


template <class Derived>
inline void ChompGradient::computeSmoothnessGradient(
                    const Trajectory & traj,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    debug_status( TAG, "computeSmoothnessGradient", "start" );
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& grad = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);

    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    //
    const MatMap & xi = (traj.isSubsampled() ? 
                         traj.getSampledXi() : 
                         traj.getXi() );

    if( use_goalset ){
        diagMul(coeffs, coeffs_goalset, xi, Ax);
    } else { 
        diagMul(coeffs, xi, Ax);
    }
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    grad = Ax + b; 

    debug_status( TAG, "computeSmoothnessGradient", "end" );
}

void ChompGradient::computeCollisionGradient(const Trajectory & traj,
                                             MatMap & grad)
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(traj, grad);
    } else {
        fextra = 0;
    }
}

void ChompGradient::computeCollisionGradient(const Trajectory & traj,
                                             MatX & grad)
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(traj, grad);
    } else {
        fextra = 0;
    }
}

}// namespace

