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

    double gamma, dt;

    //the jacobian that maps between work and configuration space
    MatX dx_dq;

    //The vector that is the collision gradient in workspace.
    MatX collision_gradient;
    MatX gradient_t;
    
    //Working variables for the collision gradient computation
    MatX q0, q1, q2;
    MatX cspace_vel,  cspace_accel,
         wkspace_vel, wkspace_accel;
    MatX P, K; 

  public:

    CollisionFunction( size_t cspace_dofs,
                       size_t workspace_dofs, 
                       size_t n_bodies,
                       double gamma);
                        
    //evaluate the gradient of the objective function
    //  at the current trajectory
    template <class Derived>
    double evaluate( const Trajectory & trajectory,
                     const Eigen::MatrixBase<Derived> & g_const);
    double evaluate( const Trajectory & trajectory );

  private:
    virtual double evaluateTimestep( int t,
                                     const Trajectory & trajectory,
                                     bool set_gradient = true );
    
    // return the cost for a given configuration/body, along with jacobians
    // q is the current configuration,
    // body_index is the index of the current body element being
    //      collision checked.
    // dx_dq is the jacobian of workspace position
    //       (workspace_DOFs X configuration_space_DOF)
    // collision_gradient is the gradient (Jacobian transpose)
    //       of cost with respect to workspace position
    //       it should be a vector of shape: (workspace_DOF X 1)
    virtual double getCost( const MatX& state,
                            size_t current_index,
                            MatX& dx_dq, 
                            MatX& collision_gradient ) = 0;

    double projectCost( double cost, bool setGradient = true);

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
    
    dt = trajectory.getDt();
    
    double total = 0.0;

    for (int t=0; t < trajectory.rows() ; ++t) {
        total += evaluateTimestep( t, trajectory );
        g.row(t) += gradient_t;
    }

    return total;
}

}//namespace 


#endif
