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

#ifndef _HMC_H_
#define _HMC_H_

#include "../utils/utils.h"
#include "../containers/Metric.h"
#include "../containers/Trajectory.h"

namespace mopt{

class HMC{

  private:
    
    // hmc_lambda : the parameter that determines the frequency, and
    //              magnitude of random resampling
    double lambda, previous_energy;

    // doNotReject : if true, the energy of the system will not
    //               be evaluated, and low energy systems will not be
    //               rejected.
    bool doNotReject;

    // hmc_resample_iter : the next iteration that resampling will be
    //                     done on
    size_t resample_iter;
    
    //old_data : the previous xi, saved from the last resample iteration,
    //          if rejection is on, this will be restored if the, 
    //          energy of the trajectory does not increase.
    //old_momentum : like old_xi, this is a saved previous state,
    //               for use if the current trajectory is rejected.
    MatX old_momentum;
    double * old_data;
    
    static const std::string TAG;

  //PUBLIC MEMBER FUNCIONS
  public:


    HMC( double lambda=0.02, bool doNotReject=true );

    ~HMC();    
    

    //SHould be called each iteration of the optimizer,
    //  if it is on a resample iteration,
    //  it will resample the momentum.
    void iterate( size_t current_iteration,
                  double lastObjective,
                  const Metric & metric,
                  Trajectory & traj,
                  MatX & momentum);
    
    //setup the random seed for HMC
    void setSeed(unsigned long seed=0);

  private: 
    
    
    //checks the current HMC iteration, and rejects it if the 
    //  energy of the system is too low.
    bool checkForRejection( Trajectory & traj,
                            MatX & momentum,
                            double lastObjective);
    
    //samples a random momentum from the probability dist given by:
    //  exp( -0.5*xAx )
    void getRandomMomentum( const Metric & metric,
                            size_t current_iteration,
                            MatX & momentum );
        
};

}//namespace

#endif
