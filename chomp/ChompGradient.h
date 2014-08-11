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
    virtual double addToGradient(ConstMatMap& xi, const MatX& pinit,
                                 const MatX& pgoal, double dt,
                                 MatMap& g) =0;
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
    virtual double addToGradient(ConstMatMap& xi, const MatX& pinit,
                                 const MatX& pgoal, double dt,
                                 MatMap& g);

  private:
    //this function does all of the actual work.
    template< class Derived1, class Derived2, class Derived3>
    double computeGradient(const Eigen::MatrixBase<Derived1> & xi,
                           const Eigen::MatrixBase<Derived2> & pinit,
                           const Eigen::MatrixBase<Derived2> & pgoal,
                           double dt,
                           const Eigen::MatrixBase<Derived3> & g_const);

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
public:
    
    ChompGradientHelper* ghelper;

    ChompObjectiveType objective_type;
    
    int M; // degrees of freedom
    int N; // number of timesteps
    
    int iteration;

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
    
    ChompGradient( const MatX& pinit, const MatX& pgoal, 
                   ChompObjectiveType objective_type=MINIMIZE_ACCELERATION,
                   double total_time=1.0);

    ~ChompGradient(){}
    
    //prepares chomp to be run at a resolution level
    void prepareRun( int N,
                     bool use_goalset=false,
                     bool subsample=false );

    MatX& getInvAMatrix( bool subsample=false);
   
    MatX& getGradient( const MatX & xi);

    MatX& getSmoothnessGradient(const MatX & xi);

    MatX& getCollisionGradient( const MatX & xi);

    MatX& getSubsampledGradient(int N_sub);

    double getGradient( unsigned n_by_m, const double * xi, double * grad);


    static double NLoptFunction(unsigned n_by_m, const double* x,
                                double* grad, void *data) 
    {
        return reinterpret_cast<ChompGradient*>(data)
                            ->getGradient(n_by_m, x, grad);
    }
    
    template<class Derived>
    inline MatX getTick( int tick,
                         const Eigen::MatrixBase<Derived> & xi) const
    {
        return getTickBorderRepeat( tick, xi, q0, q1, dt );
    }

    // evaluates the objective function for cur. thing.
    // only works if prepareChompIter has been called since last
    // modification of xi.
    template<class Derived>
    inline double evaluateObjective(
                    const Eigen::MatrixBase<Derived> & xi) const
    {
        return (0.5 * mydot(xi, Ax) + mydot(xi, b) + c) + fextra;
    }

  private:
    template <class Derived1, class Derived2>
    void computeSmoothnessGradient( const Eigen::MatrixBase<Derived1> & xi,
                                const Eigen::MatrixBase<Derived2> & g );
    
    void computeCollisionGradient(const MatX & xi, MatX & grad);
    void computeCollisionGradient(ConstMatMap & xi, MatMap & grad);

};


}//namespace chomp


#endif
