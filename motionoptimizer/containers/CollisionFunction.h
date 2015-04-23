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

#ifndef _COLLISION_FUNCTION_H_
#define _COLLISION_FUNCTION_H_

#include "../utils/utils.h"
#include "Trajectory.h"
#include "Metric.h"

namespace mopt {

class CollisionFunction {

  private:
    static const char* TAG;

    size_t configuration_space_DOF;
    size_t workspace_DOF;
    size_t number_of_bodies;

    double gamma;

    //the jacobian that maps between work and configuration space
    MatX dx_dq;

    //The vector that is the collision gradient in workspace.
    MatX cgrad;
    
    //Working variables for the collision gradient computation
    MatX q0, q1, q2;
    MatX cspace_vel,  cspace_accel,
         wkspace_vel, wkspace_accel;
    MatX P, K; 

  public:
    // return the cost for a given configuration/body, along with jacobians
    // q is the current configuration,
    // body_index is the index of the current body element being
    //      collision checked.
    // dx_dq is the jacobian of workspace position
    //       (workspace_DOFs X configuration_space_DOF)
    // cgrad is the gradient (Jacobian transpose)
    //       of cost with respect to workspace position
    //       it should be a vector of shape: (workspace_DOF X 1)
    typedef double (*CostFunction)(const MatX&, size_t,
                                    MatX&,  MatX&, void * );
  private:
    CostFunction cost_function;
    void * cost_function_data; 

  public:

    CollisionFunction( size_t cspace_dofs,
                       size_t workspace_dofs, 
                       size_t n_bodies,
                       double gamma,
                       CostFunction func,
                       void * data = NULL );
                        
    //evaluate the gradient of the objective function
    //  at the current trajectory
    template <class Derived>
    double evaluate( const Trajectory & trajectory,
                     const Eigen::MatrixBase<Derived> & g_const);

    double evaluate( const Trajectory & trajectory );

};

//Include this single templated function because otherwise we would need
//  an entire extra file
template< class Derived >
double CollisionFunction::evaluate(
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

        for (size_t u=0; u < number_of_bodies; ++u) {

            float cost = (*cost_function)(q1, u, dx_dq, cgrad,
                                      cost_function_data);

            debug_assert( size_t(dx_dq.rows()) == workspace_DOF );
            debug_assert( size_t(dx_dq.cols()) == configuration_space_DOF );

            if (cost > 0.0) {

                wkspace_vel = dx_dq * cspace_vel;

                //this prevents nans from propagating.
                //   Several lines below,  wkspace_vel /= wv_norm
                //   if wv_norm is zero, nans propogate.
                if (wkspace_vel.isZero()){ continue; }

                wkspace_accel = dx_dq * cspace_accel;
              
                float wv_norm = wkspace_vel.norm();
                wkspace_vel /= wv_norm;

                // add to total
                double scl = wv_norm * gamma * trajectory.getDt();

                total += cost * scl;
              
                P = MatX::Identity(workspace_DOF, workspace_DOF)
                    - (wkspace_vel * wkspace_vel.transpose());

                K = (P * wkspace_accel) / (wv_norm * wv_norm);

                // scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
                g.row(t) += (scl * (dx_dq.transpose() *
                          (P * cgrad - cost * K)).transpose());
            }
        }
    }

    return total;
}

}//namespace 


#endif
