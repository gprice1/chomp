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

#ifndef _CHOMP_GRADIENT_H_
#define _CHOMP_GRADIENT_H_

#include "chomputil.h"
#include <vector>

namespace chomp {

class ChompGradientHelper {
  public:
    virtual ~ChompGradientHelper();

    virtual double addToGradient(const MatX& xi, const MatX& pinit,
                                 const MatX& pgoal, double dt,
                                 MatX& g) =0;
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

    virtual double addToGradient(const MatX& xi, const MatX& pinit,
                                 const MatX& pgoal, double dt,
                                 MatX& g);

    ChompCollisionHelper* chelper;
    double gamma;

    MatX dx_dq, cgrad;

    MatX q0, q1, q2;
    MatX cspace_vel,  cspace_accel,
         wkspace_vel, wkspace_accel;
    MatX P, K; 
    
};

class ChompGradient {
public:
    
    const Chomp & chomper; // this is for sending to the gradient helper

    ChompGradientHelper* ghelper;

    ChompObjectiveType objective_type;
    
    int M; // degrees of freedom
    int N; // number of timesteps

    MatX q0; // initial point of size M
    MatX q1; // end point of size M

    MatX Ax; // A*x of size N-by-M
  
    MatX coeffs; // coeffs for A e.g. [1, -4, 6] of size D-by-1 for accel
    MatX coeffs_sub; // coeffs for downsampled A e.g. [1, 6] for accel
    MatX coeffs_goalset;

    double fscl; // dynamic scaling factor for f, e.g. (N+1)*(N+2) for accel
    double fextra; // extra objective function from gradient helper

    double t_total; // total time for (N+1) timesteps
    double dt; // computed automatically from t_total and N
    double inv_dt; // computed automatically from t_total and N

    MatX L, L_sub; // skyline Cholesky coeffs of A of size N-by-D

    MatX g; // gradient terms (Ax + b) of size N-by-M
    MatX g_sub;

    // working variables
    MatX H_trans, P, P_trans, HP, Y, W, g_trans, delta, delta_trans; 

    MatX b; // endpoint vectors for this problem of size N-by-M


    double c; // c value for objective function

    bool use_goalset;
    
    ChompGradient( const Chomp & chomper,
                   const MatX& pinit, const MatX& pgoal, 
                   ChompObjectiveType objective_type ,
                   double total_time );

    ~ChompGradient(){}
    
    //prepares chomp to be run at a resolution level
    void prepareRun( int N,
                     bool use_goalset=false,
                     bool subsample=false );

    MatX& getInvAMatrix( bool subsample=false);
   
    MatX& getGradient( const MatX & xi);

    MatX& getGradientThroughLMatrix(const MatX & xi );

    MatX& getSmoothnessGradient(const MatX & xi);

    MatX& getCollisionGradient( const MatX & xi);

    MatX& getSubsampledGradient(int N_sub);

    // evaluates the objective function for cur. thing.
    // only works if prepareChompIter has been called since last
    // modification of xi.
    double evaluateObjective( const MatX & xi) const ;

    inline void addCollisionGradient(const MatX & xi){

        //If there is a gradient helper, add in the contribution from
        //  that source, and set the fextra variable to the cost
        //  associated with the additional gradient.
        if (ghelper) {
            fextra = ghelper->addToGradient(xi, q0, q1, dt, g);
        } else {
            fextra = 0;
        }

    }

    inline MatX getTick( int tick, const MatX & xi ){
        return getTickBorderRepeat( tick, xi, q0, q1, dt );
    }
    

};

}//namespace chomp


#endif
