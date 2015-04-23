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


#include "Gradient.h"


namespace mopt {

GradientHelper::~GradientHelper() {}

CollisionHelper::CollisionHelper(size_t nc,
                                           size_t nw,
                                           size_t nb):
    ncspace(nc), nwkspace(nw), nbodies(nb) {}

CollisionHelper::~CollisionHelper() {}

CollGradHelper::CollGradHelper(CollisionHelper* h,
                                           double g):
    chelper(h), gamma(g)
{
    dx_dq = MatX(h->nwkspace, h->ncspace);
    cgrad = MatX(h->nwkspace, 1);
}

CollGradHelper::~CollGradHelper() {}


template< class Derived >
inline double CollGradHelper::computeGradient(
                    const Trajectory & trajectory,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
    q1 = trajectory.getTick( -1 ).transpose();
    q2 = trajectory.getTick( 0  ).transpose();
    
    const double inv_dt = 1/trajectory.getDt();
    const double inv_dt_squared = inv_dt * inv_dt;

    double total = 0.0;

    for (int t=0; t < trajectory.rows() ; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = trajectory.getTick( t+1 ).transpose();

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
          double scl = wv_norm * gamma * trajectory.getDt();

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
double CollGradHelper::addToGradient( const Trajectory & trajectory,
                                          MatX& g)
{
    return computeGradient( trajectory, g );
}

double CollGradHelper::addToGradient( const Trajectory & trajectory,
                                          MatMap& g)
{
    return computeGradient( trajectory, g );
}

Gradient::Gradient() :
    ghelper(NULL)
{}


const char* Gradient::TAG = "Gradient";

void Gradient::prepareRun(const Trajectory & trajectory,
                               bool use_goalset,
                               bool is_covariant)
{
    debug_status( TAG, "prepareRun", "start" );

    //resize the g, b, and ax matrices.
    const int N = trajectory.fullN();
    const int M = trajectory.M();

    Ax.resize(N,M);
    b.resize(N,M);
    
    //get the b matrix, and get its contribution to the
    //  objective function.
    //  Also, compute the L matrix (lower triangluar cholesky 
    //  decomposition). 
    metric.initialize( N, trajectory.getObjectiveType(),
                       false,
                       use_goalset );
    c = metric.createBMatrix(trajectory.getStartPoint(),
                             trajectory.getEndPoint(),
                             b,
                             trajectory.getDt());
    
    if( trajectory.isSubsampled() ){
        subsampled_metric.initialize( N,
                                      trajectory.getObjectiveType(),
                                      true);
    }

    if( is_covariant ){ metric.multiplyLowerInverse( b ); }
    
    debug_status( TAG, "prepareRun", "end" );
}

    


double Gradient::evaluateObjective(const Trajectory & trajectory,
                                   bool is_covariant) const
{
    if ( is_covariant ){
        return 0.5 * mydot( trajectory.getFullXi(), trajectory.getFullXi() )
               + mydot( trajectory.getFullXi(), b ) + c + fextra;
    }
    return 0.5 * mydot(trajectory.getFullXi(), Ax) + 
           mydot(trajectory.getFullXi(), b) +
           c + fextra;
}

}// namespace

