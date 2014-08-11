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

#include "chomputil.h"
#include <iomanip>

namespace chomp {

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

int ChompObserver::notify(const ChompOptimizerBase&, 
                            ChompEventType,
                            size_t,
                            double, double, double) { 
    return 0; 
}

DebugChompObserver::~DebugChompObserver() {}

int DebugChompObserver::notify(const ChompOptimizerBase& c, 
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
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";

    if (std::isnan(curObjective) || std::isinf(curObjective) ||
        std::isnan(lastObjective) || std::isinf(lastObjective)) {
        return 1;
    }
    return 0;
}

void upsampleTrajectory(const MatX & xi,
                        const MatX & q0,
                        const MatX & q1,
                        double dt,
                        ChompObjectiveType objective_type,
                        MatX & xi_up)
{

  const int N = xi.rows();
  const int M = xi.cols();
  const int N_up = 2*N+1; // e.g. size 3 goes to size 7
  
  xi_up.resize(N_up, M);

  // q0    d0    d1    d2    q1   with n = 3
  // q0 u0 u1 u2 u3 u4 u5 u6 q1   with n = 7
  //
  // u0 = 0.5*(q0 + d0)
  // u1 = d0
  // u2 = 0.5*(d0 + d1)
  // u3 = d1 
  // u4 = 0.5*(d1 + d2)
  // u5 = d2
  // u6 = 0.5*(d2 + q1)

  for (int t=0; t<N_up; ++t) { // t is timestep in new (upsampled) regime

    if (t % 2 == 0) {

      assert(t == N_up-1 || (t/2) < xi.rows());
      assert(t < xi_up.rows());

      if (objective_type == MINIMIZE_VELOCITY) {

        MatX qneg1 = getTickBorderRepeat(t/2-1, xi, q0, q1, dt);
        MatX qpos1 = getTickBorderRepeat(t/2,   xi, q0, q1, dt);
        xi_up.row(t) = 0.5 * (qneg1 + qpos1);

      } else { 

        MatX qneg3 = getTickBorderRepeat(t/2-2, xi, q0, q1, dt);
        MatX qneg1 = getTickBorderRepeat(t/2-1, xi, q0, q1, dt);
        MatX qpos1 = getTickBorderRepeat(t/2,   xi, q0, q1, dt);
        MatX qpos3 = getTickBorderRepeat(t/2+1, xi, q0, q1, dt);

        const double c3 = -1.0/160;
        const double c1 = 81.0/160;

        xi_up.row(t) = c3*qneg3 + c1*qneg1 + c1*qpos1 + c3*qpos3;

      }


    } else {
      xi_up.row(t) = xi.row(t/2);
    }

  }

}

void createInitialTraj( const MatX & q0, const MatX & q1, 
                        int N, ChompObjectiveType objective_type,
                        MatX & xi)
{
    assert( q0.size() == q1.size() );
    xi.resize( N, q0.size() );
    
    createInitialTraj( q0, q1, objective_type, xi );
}


void createInitialTraj( const MatX & q0, const MatX & q1, 
                        ChompObjectiveType objective_type,
                        MatX & xi)
{
    assert( q0.size() == q1.size() );
    assert( xi.cols() == q0.size() );
    
    const int N = xi.rows();

    if ( objective_type == MINIMIZE_VELOCITY ){
        //if the goal is Minimize Velocity, simply linearly interpolate
        //  from the start to goal state.
        for ( int i=0; i < N; i ++ ){

            double t = double(i+1) / double(N+1);
            xi.row( i ) = q0 + (q1-q0)*t;
        }
    }
    //we are minimizing acceleration, so do some math to get
    //   a low acceleration initial traj.
    else{

        MatX q_mid = (q0 + q1)/2;

        //apply a constant positive acceleration from 
        //  the start to the midpoint.
        for ( int i=0; i < (N+1)/2; i ++ ){
            double t = double(i+1)/double(N/2 + 1);
            xi.row( i ) = q0 + (q_mid-q0)*t*t;
        }

        //apply a constant negative acceleration from
        //  the midpoint to the endpoint.
        for ( int i = N-1; i > N/2; i -- ){
            double t = double(N-i)/double(N/2 + 1);
            xi.row( i ) = q1 + (q_mid-q1)*t*t;
        }
    }
}



}//namespace
