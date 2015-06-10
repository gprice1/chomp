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
    
  protected:
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

    virtual ~CollisionFunction(){}

    //evaluate the gradient of the objective function
    //  at the current trajectory
    template <class Derived>
    double evaluate( const Trajectory & trajectory,
                     const Eigen::MatrixBase<Derived> & g_const);
    double evaluate( const Trajectory & trajectory );

    size_t getNumberOfBodies() const { return number_of_bodies; }
    size_t getWorkspaceDOF() const { return workspace_DOF; }
    size_t getConfigurationSpaceDOF() const { return configuration_space_DOF; }

    void setNumberOfBodies( size_t size ){ number_of_bodies = size; }

  protected:
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
                            MatX& collision_gradient ){ return 0;};

    template< class Derived1, class Derived2 >
    double projectCost( double cost,
                        const Eigen::MatrixBase<Derived1> & jacobian,
                        const Eigen::MatrixBase<Derived2> & coll_grad,
                        bool set_gradient );

  private:
    void prepareEvaluation( const Trajectory & trajectory );

};

//Include this single templated function because otherwise we would need
//  an entire extra file
template< class Derived >
double CollisionFunction::evaluate(
                    const Trajectory & trajectory,
                    const Eigen::MatrixBase<Derived> & g_const)
{

    debug_status( TAG, "evaluate", "start");

    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
    prepareEvaluation( trajectory );

    double total = 0.0;

    for (int t=0; t < trajectory.rows() ; ++t) {
        q0 = q1;
        q1 = q2;
        q2 = trajectory.getTick( t+1 ).transpose();

        cspace_vel = 0.5 * (q2 - q0) / dt;        
        cspace_accel = (q0 - 2.0*q1 + q2) / (dt * dt);

        gradient_t.setZero();

        total += evaluateTimestep( t, trajectory );
        g.row(t) += gradient_t;
    }

    debug_status( TAG, "evaluate", "start");
    
    return total;
}

}//namespace 


#endif
