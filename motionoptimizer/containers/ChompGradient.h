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

#include "../utils/utils.h"
#include "Trajectory.h"

namespace chomp {

class ChompGradientHelper {
  public:
    virtual ~ChompGradientHelper();

    virtual double addToGradient( const Trajectory & traj, MatX& g) =0;
    virtual double addToGradient( const Trajectory & traj, MatMap& g) =0;
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

    virtual double addToGradient(const Trajectory & traj, MatX& g);
    virtual double addToGradient(const Trajectory & traj, MatMap& g);

  private:
    //this function does all of the actual work.
    template< class Derived>
    double computeGradient(const Trajectory & traj,
                           const Eigen::MatrixBase<Derived> & g_const);

  public:
    ChompCollisionHelper* chelper;
    double gamma;

    MatX dx_dq, cgrad;

    MatX q0, q1, q2;
    MatX cspace_vel,  cspace_accel,
         wkspace_vel, wkspace_accel;
    MatX P, K; 
    
};

class ChompGradient {

  private:
    Trajectory & trajectory;
    static const char* TAG;
    
  private:
    
    ChompGradientHelper* ghelper;
    
    //TODO decide to keep this here and take it out of trajectory,
    //  or remove this. 
    
    MatX Ax; // A*x of size N-by-M
  
    MatX coeffs; // coeffs for A e.g. [1, -4, 6] of size D-by-1 for accel
    MatX coeffs_sub; // coeffs for downsampled A e.g. [1, 6] for accel
    MatX coeffs_goalset;

    double fscl; // dynamic scaling factor for f, e.g. (N+1)*(N+2) for accel
    double fextra; // extra objective function from gradient helper

    MatX L, L_sub; // skyline Cholesky coeffs of A of size N-by-D

    MatX g; // gradient terms (Ax + b) of size N-by-M
    MatX g_sub;

    MatX b; // endpoint vectors for this problem of size N-by-M

    double c; // c value for objective function

    bool use_goalset;
    
  public:
    ChompGradient( Trajectory & traj );

    ~ChompGradient(){}
    
    inline void setGradientHelper( ChompGradientHelper* help ){
        ghelper = help;
    }

    inline ChompGradientHelper * getGradientHelper(){ return ghelper; }

    //prepares chomp to be run at a resolution level
    void prepareRun(bool use_goalset=false);

    const MatX& getInvAMatrix();
   
    MatX& getGradient();

    MatX& getSmoothnessGradient();

    MatX& getCollisionGradient();

    MatX& getSubsampledGradient(int N_sub);

    double getGradient( unsigned n_by_m, 
                        const double * xi,
                        double * grad);
    
    // evaluates the objective function for cur. thing.
    // only works if prepareChompIter has been called since last
    // modification of xi.
    double evaluateObjective() const;
    
  private:
    template <class Derived>
    void computeSmoothnessGradient(const Eigen::MatrixBase<Derived> & g );
    
    void computeCollisionGradient( MatMap & grad);
    void computeCollisionGradient( MatX & grad  );

};


}//namespace chomp


#endif
