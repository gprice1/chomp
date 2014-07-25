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

#define debug if (0) std::cout
#define debug_assert if (0) assert

namespace chomp {

ChompGradient::ChompGradient( const MatX& pinit, 
                              const MatX& pgoal, 
                              ChompObjectiveType objective_type ){
    q0( pinit ), q1( pgoal ),
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

void ChompGradient::prepareRun(int N,
                               bool using_goalset,
                               bool subsample)
{
    
    //if we are using a goal set, increase the size of 
    //  N by one.
    if (using_goalset){ N++; }
    
    //resize, and set the b matrix to zero
    b.resize(N,M);
    b.setZero();

    //get the b matrix, and get their contribution to the
    //  objective function
    if (using_goalset){
        skylineChol(N, coeffs, L);
        c = createBMatrix(N, coeffs, q0, q1, b, dt);
    }else {
        prepareNonGoalSetRun();
        c = createBMatrix(N, coeffs, q0, b, dt);
    }
    
    //resize the g and ax matrices.
    g.resize(N,M);
    Ax.resize(N,M);

    if (subsample) {
        int N_sub = (N+1)/2;
        skylineChol(N_sub, coeffs_sub, L_sub); 
    }
}


MatX& ChompGradient::getInvAMatrix( bool subsample){
    return (subsample ? L_sub : L );
}

template <class Derived>
MatX& ChompGradient::getCollisionGradient( 
                    const Eigen::MatrixBase<Derived>& xi )
{
        //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    g.setZero();

    if (ghelper) {
        fextra = ghelper->addToGradient(*this, g);
    } else {
        fextra = 0;
    }
    return g;

}


template <class Derived>
MatX& ChompGradient::getGradient( 
                           const Eigen::MatrixBase<Derived>& xi )
{
    
    getSmoothnessGradient( xi );

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(*this, g);
    } else {
        fextra = 0;
    }

    return g;

}

template <class Derived>
MatX& ChompGradient::getSmoothnessGradient( 
                     const Eigen::MatrixBase<Derived>& xi )
{
    
    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    if( use_goalset ){ diagMul(coeffs, coeffs_goalset, xi, Ax); }
    else { diagMul(coeffs, xi, Ax); }
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    g = Ax + b;

    return g;
}

template <class Derived>
SubMatMap& ChompGradient::getSubsampledGradient( 
                          const Eigen::MatrixBase<Derived>& xi,
                          int N_sub)
{
    
    //get the contribution from the smoothness.
    getSmoothnessGradient( xi );

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        fextra = ghelper->addToGradient(*this, g );
        
        //TODO make this line below possible
        //fextra = ghelper->addToGradient(*this, g, sub_factor );
    } else {
        fextra = 0;
    }
    
    g_sub = SubMatMap( g.data(), N_sub, M );

    return g_sub;

}


// evaluates the objective function for cur. thing.
//
// only works if prepareChompIter has been called since last
// modification of xi.
double Chomp::evaluateObjective( const MatX & xi ) const {

  /*

    K is (n+1)-by-n

    K = [  1  0  0  0 ... 0  0  0
    -2  1  0  0 ... 0  0  0 
    1 -2  1  0 ... 0  0  0
    0  1 -2  1 ... 0  0  0
    ...
    0  0  0  0 ... 1 -2  1 
    0  0  0  0 ... 0  1 -2 
    0  0  0  0 ... 0  0  1 ]

    e = [ -x0  x0   0 ... 0  x1 -x1 ]^T

    ||Kx + e||^2 = 
     
  */
  
    if ( false ){
        const double xi_Ax = 0.5 * mydot( xi, Ax );
        const double xi_b = mydot( xi, b );
        const double smoothness = ( xi_Ax + xi_b + c ) * fscl;
        std::cout <<  "xi * Ax    = " << xi_Ax << std::endl;
        std::cout <<  "xi * b     = " << xi_b << std::endl;
        std::cout <<  "c          = " << c << std::endl;
        std::cout <<  "fextra     = " << fextra << std::endl;
        std::cout <<  "fscale     = " << fscl << std::endl;
        std::cout <<  "smoothness = " << smoothness << std::endl;

        return smoothness + fextra;
    }

    return (0.5 * mydot( xi, Ax ) + mydot( xi, b ) + c ) + fextra;

}

}// namespace

