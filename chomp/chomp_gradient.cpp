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

#include "Chomp.h"
#include "ConstraintFactory.h"
#include "Constraint.h"
#include "HMC.h"
#include <float.h>
#include <cmath>
#include <iomanip>

#define debug if (0) std::cout
#define debug_assert if (0) assert

namespace chomp {

  bool isnan( MatX & mat ){
      for ( int i = 0; i < mat.rows() ; i ++ ){
          for ( int j = 0; j < mat.cols(); j ++ ){
              if ( mat(i,j) != mat(i,j) ){
                  return true;
              }
          }
      }
      return false;
  }

  const char* eventTypeString(int eventtype) {
    switch (eventtype) {
    case CHOMP_INIT: return "CHOMP_INIT";
    case CHOMP_GLOBAL_ITER: return "CHOMP_GLOBAL_ITER";
    case CHOMP_LOCAL_ITER: return "CHOMP_LOCAL_ITER";
    case CHOMP_FINISH: return "CHOMP_FINISH";
    case CHOMP_TIMEOUT: return "CHOMP_TIMEOUT";
    case CHOMP_GOALSET_ITER: return "CHOMP_GOALSET_ITER";
    default: return "[INVALID]";
    }
  }

  ChompObserver::~ChompObserver() {}

  int ChompObserver::notify(const Chomp&, 
                            ChompEventType,
                            size_t,
                            double, double, double) { 
    return 0; 
  }

  DebugChompObserver::~DebugChompObserver() {}

  int DebugChompObserver::notify(const Chomp& c, 
                                 ChompEventType e,
                                 size_t iter,
                                 double curObjective, 
                                 double lastObjective,
                                 double constraintViolation) { 
    std::cout << "chomp debug: "
              << "event=" << eventTypeString(e) << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10) << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10) << constraintViolation << "\n";

    if (std::isnan(curObjective) || std::isinf(curObjective) ||
        std::isnan(lastObjective) || std::isinf(lastObjective)) {
      return 1;
    }
    return 0;
  }

  ChompGradientHelper::~ChompGradientHelper() {}

  ChompCollisionHelper::ChompCollisionHelper(size_t nc, size_t nw, size_t nb):
    ncspace(nc), nwkspace(nw), nbodies(nb) {}

  ChompCollisionHelper::~ChompCollisionHelper() {}

  ChompCollGradHelper::ChompCollGradHelper(ChompCollisionHelper* h,
                                           double g):
    chelper(h), gamma(g)
  {
    dx_dq = MatX(h->nwkspace, h->ncspace);
    cgrad = MatX(h->nwkspace, 1);
  }

  ChompCollGradHelper::~ChompCollGradHelper() {}

  double ChompCollGradHelper::addToGradient(const Chomp& c,
                                            MatX& g) {

    q1 = c.getTickBorderRepeat(-1).transpose();
    q2 = c.getTickBorderRepeat(0).transpose();

    double total = 0.0;

    for (int t=0; t<c.N; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = c.getTickBorderRepeat(t+1).transpose();

      cspace_vel = 0.5 * (q2 - q0) * c.inv_dt;        
      cspace_accel = (q0 - 2.0*q1 + q2) * (c.inv_dt * c.inv_dt);

      for (size_t u=0; u < chelper->nbodies; ++u) {

        float cost = chelper->getCost(q1, u, dx_dq, cgrad);
        if (cost > 0.0) {

          wkspace_vel = dx_dq * cspace_vel;

          //this prevents nans from propagating. Several lines below, 
          //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
          if (wkspace_vel.isZero()){ continue; }

          wkspace_accel = dx_dq * cspace_accel;
          
          float wv_norm = wkspace_vel.norm();
          wkspace_vel /= wv_norm;

          // add to total
          double scl = wv_norm * gamma / c.inv_dt;

          total += cost * scl;
          
          P = MatX::Identity(chelper->nwkspace, chelper->nwkspace)
              - (wkspace_vel * wkspace_vel.transpose());

          K = (P * wkspace_accel) / (wv_norm * wv_norm);
         

          //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
          g.row(t) += (scl * (dx_dq.transpose() *
                      (P * cgrad - cost * K)).transpose());
         

        }
      }
    }
    return total;

  }


 


  //this sets up an iteration of 'normal' chomp, as opposed to
  //    goal set chomp.
  void Chomp::prepareStandardChomp(){
    skylineChol(N, coeffs, L);

    b.resize(N,M);
    b.setZero();

    c = createBMatrix(N, coeffs, q0, q1, b, dt);

    g.resize(N,M);
  }

  void Chomp::prepareChomp() {
    
    if (objective_type == MINIMIZE_VELOCITY) {

      coeffs.resize(1,2);
      coeffs_sub.resize(1,1);

      coeffs << -1, 2;
      coeffs_sub << 2;

    } else {
    
      coeffs.resize(1,3);
      coeffs_sub.resize(1,2);
      
      coeffs << 1, -4, 6;
      coeffs_sub << 1, 6;

    }
    
    dt = t_total / (N+1);
    inv_dt = (N+1) / t_total;

    if (objective_type == MINIMIZE_VELOCITY) {
      fscl = inv_dt*inv_dt;
    } else {
      fscl = inv_dt*inv_dt*inv_dt;
    }

    cur_global_iter = 0;
    cur_local_iter = 0;

    clearConstraints();

    if (factory) {
      factory->getAll(N, constraints);
    }
    

    if( use_goalset ){ prepareGoalSet(); }
    else { prepareStandardChomp(); }


    // decide whether base case or not
    bool subsample = ( N > minN && !use_goalset );
    if (full_global_at_final && N >= maxN) {
      subsample = false;
    }

    if (!subsample) {
      N_sub = 0;
      if ( use_momentum ){
        momentum.resize( N, M );
        momentum.setZero();
      }
      if( hmc ){ hmc->setupRun(); }

    } else {
      N_sub = (N+1)/2;
      g_sub.resize(N_sub, M);
      xi_sub.resize(N_sub, M);
      skylineChol(N_sub , coeffs_sub, L_sub); 
    }


  }

  // evaluates the objective function for cur. thing.
  //
  // only works if prepareChompIter has been called since last
  // modification of xi.
  double Chomp::evaluateObjective() {

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

  void Chomp::prepareGoalSet(){
    
    //resize xi, and add q1 into it.
    xi.conservativeResize( xi.rows() + 1, xi.cols() );
    xi.row( xi.rows() - 1 ) = q1;
    int n = xi.rows();

    if (objective_type == MINIMIZE_VELOCITY) {
      goalset_coeffs.resize(1,1);
      goalset_coeffs << 1;

    } else {
      goalset_coeffs.resize(2,2);
      goalset_coeffs << 6, -3,
                       -3,  2 ;
    }
    
    constraints.push_back( goalset );

    skylineChol( n, coeffs, goalset_coeffs, L);

    b.resize(n,M);
    b.setZero();

    c = createBMatrix(n, coeffs, q0, b, dt);

    g.resize(n,M);

    //if we are doing goal set chomp, we should not subsample
    N_sub = 0;

    N = n;
  }

  void Chomp::finishGoalSet(){
    
    std::cout << "Finishing goal set" << std::endl;

    use_goalset = false;

    q1 = xi.row( xi.rows() - 1 );
    xi.conservativeResize( xi.rows() -1, xi.cols() );
    
    N = xi.rows();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    constraints.pop_back();
    
    std::cout << "Preparing chomp" <<std::endl;

    prepareChomp();
    std::cout << "Done Preparing chomp" <<std::endl;
  }



}// namespace

