/*
* Copyright (c) 2008-2015, Matt Zucker and Temple Price
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

#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "mzcommon/TimeUtil.h"
#include "OptimizerBase.h"
#include "HMC.h"

namespace mopt {

class ChompOptimizerBase : public OptimizerBase{
    
  public:
    
    //either global or local iteration, depending
    //on the type of the optimization, this should be assigned in 
    // inheriting class constructor    
    EventType event; 

    static const char* TAG;

    double alpha;       // the gradient step size

    MatX bounds_violations;
    
    //timeout_seconds : the amount of time from the start of chomp
    //                  to a forced timeout.
    //canTimeout : is timing out a possible termination condition?
    // didTimeout : has chomp timed out ?
    // stop_time : the time at which chomp will timeout.
    bool canTimeout, didTimeout;
    TimeStamp stop_time;

    //TODO set to zero.
    size_t min_iter;
   
    bool use_momentum;
    MatX g, momentum;
    
    //an HMC object for performing the Hamiltonian Monte Carlo method
    HMC * hmc;
    
    ChompOptimizerBase(ProblemDescription & problem,
                       Observer * observer = NULL,
                       double obstol = 1e-8,
                       double timeout_seconds = 0,
                       size_t max_iter = size_t(-1)); 

    virtual ~ChompOptimizerBase(){};
    
    void solve();

    inline void setAlpha( double a ){ alpha = a; }
    inline double getAlpha( ){ return alpha ; }

  protected:

    virtual void optimize()=0;

  private:
    
    //Checks the bounds of chomp, and smoothly pushes the trajectory
    //  back into the bounds.
    void checkBounds();
    
    //called for every interation. It returns true, if 
    //  the optimization is not finished.
    bool iterate();

    //check if chomp is finished.
    bool checkFinished(EventType event);
    
    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective );

};


}//namespace

#endif
