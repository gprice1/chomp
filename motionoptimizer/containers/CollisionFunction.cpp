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


#include "CollisionFunction.h"


namespace mopt {

const char* CollisionFunction::TAG = "CollisionFunction";

CollisionFunction::CollisionFunction( size_t cspace_dofs,
                                      size_t workspace_dofs, 
                                      size_t n_bodies,
                                      double gamma) :
    configuration_space_DOF( cspace_dofs ),
    workspace_DOF( workspace_dofs ),
    number_of_bodies( n_bodies ),
    gamma( gamma ),
    dx_dq( workspace_dofs, cspace_dofs ),
    collision_gradient( workspace_dofs, 1 ),
    gradient_t( 1, cspace_dofs )
{
}

double CollisionFunction::evaluate( const Trajectory & trajectory )
{
    
    q1 = trajectory.getTick( -1 ).transpose();
    q2 = trajectory.getTick( 0  ).transpose();
    
    dt = trajectory.getDt();

    double total = 0.0;

    for (int t=0; t < trajectory.rows() ; ++t) {

        total += evaluateTimestep( t, trajectory, false );

    }
    return total;
}


double CollisionFunction::evaluateTimestep( int t, 
                                            const Trajectory & trajectory,
                                            bool set_gradient )
{

    q0 = q1;
    q1 = q2;
    q2 = trajectory.getTick( t+1 ).transpose();

    cspace_vel = 0.5 * (q2 - q0) / dt;        
    cspace_accel = (q0 - 2.0*q1 + q2) / (dt * dt);

    gradient_t.setZero();
    
    double total = 0;

    for (size_t u=0; u < number_of_bodies; ++u) {

        float cost = getCost(q1, u, dx_dq, collision_gradient);

        debug_assert( size_t(dx_dq.rows()) == workspace_DOF );
        debug_assert( size_t(dx_dq.cols()) == configuration_space_DOF );

        if (cost > 0.0) {
            total += projectCost( cost, set_gradient );
        }
    }
    
    return total;
}


double CollisionFunction::projectCost( double cost, bool set_gradient ){
    
    wkspace_vel = dx_dq * cspace_vel;

    //this prevents nans from propagating.
    //   Several lines below,  wkspace_vel /= wv_norm
    //   if wv_norm is zero, nans propogate.
    if (wkspace_vel.isZero()){ return 0; }
    
    float wv_norm = wkspace_vel.norm();
    
    double scl = wv_norm * gamma * dt;

    if ( set_gradient ){
        wkspace_vel /= wv_norm;
        
        wkspace_accel = dx_dq * cspace_accel;

        P = MatX::Identity(workspace_DOF, workspace_DOF)
            - (wkspace_vel * wkspace_vel.transpose());

        K = (P * wkspace_accel) / (wv_norm * wv_norm);

        // scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
        gradient_t += (scl * (dx_dq.transpose() *
                      (P * collision_gradient - cost * K)).transpose());
    }

    return cost * scl;
}

}// namespace

