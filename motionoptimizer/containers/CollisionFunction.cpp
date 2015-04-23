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
                                      double gamma,
                                      CostFunction func,
                                      void * data ) :
    configuration_space_DOF( cspace_dofs ),
    workspace_DOF( workspace_dofs ),
    number_of_bodies( n_bodies ),
    gamma( gamma ),
    dx_dq( workspace_dofs, cspace_dofs ),
    cgrad( workspace_dofs, 1 ),
    cost_function( func ),
    cost_function_data( data )
{
}

double CollisionFunction::evaluate( const Trajectory & trajectory )
{
    
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
        if (cost > 0.0) {

          wkspace_vel = dx_dq * cspace_vel;

          //this prevents nans from propagating. Several lines below, 
          //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
          if (wkspace_vel.isZero()){ continue; }

          float wv_norm = wkspace_vel.norm();

          // add to total
          double scl = wv_norm * gamma * trajectory.getDt();

          total += cost * scl;
         
        }
      }
    }
    return total;
}






}// namespace

