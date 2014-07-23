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

#include <iostream>
#include <vector>
#include <pthread.h>
#include "../mzcommon/TimeUtil.h"
#include "Chomp.h"

namespace chomp {


class ChompGradient {
public:

    ChompGradientHelper* ghelper;

    ChompObjectiveType objective_type;
    
    int M; // degrees of freedom

    int N; // number of timesteps

    MatX Ax; // A*x of size N-by-M
  
    MatX coeffs; // coeffs for A e.g. [1, -4, 6] of size D-by-1 for accel
    MatX coeffs_sub; // coeffs for downsampled A e.g. [1, 6] for accel
    MatX coeffs_goalset;

    double fscl; // dynamic scaling factor for f, e.g. (N+1)*(N+2) for accel

    double fextra; // extra objective function from gradient helper

    MatX L; // skyline Cholesky coeffs of A of size N-by-D

    MatX g; // gradient terms (Ax + b) of size N-by-M

    // working variables
    MatX H_trans, P, P_trans, HP, Y, W, g_trans, delta, delta_trans; 

    MatX b; // endpoint vectors for this problem of size N-by-M

    MatX q0; // initial point of size M
    MatX q1; // end point of size M

    double c; // c value for objective function

    bool use_goalset;
    
    ChompGradient( const MatX& pinit,const MatX& pgoal, 
                    ChompObjectiveType objective_type);

    ~ChompGradient(){}
    
    //prepares chomp to be run at a resolution level
    void prepareRun( size_t N,
                     bool subsample=false,
                     bool usingGoalSet=false);

    MatX& getGradient( const MatX& trajectory );
    MatX& getGradientThroughLMatrix( const MatX& trajectory);
    MatX& getSmoothnessGradient( const MatX& trajectory );
    MatX& getCollisionGradient( const MatX& trajectory );

    // get the tick, respecting endpoint repetition
    MatX getTickBorderRepeat(int tick) const;

    // evaluates the objective function for cur. thing.
    // only works if prepareChompIter has been called since last
    // modification of xi.
    double evaluateObjective( const MatX & xi) const ;

};

}


#endif