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


#include "HMC.h"
#include "mzcommon/mersenne.h"
#include "mzcommon/gauss.h"


namespace mopt {

const std::string HMC::TAG = "HMC";


HMC::HMC( double lambda, bool doNotReject):
    lambda( lambda ),
    previous_energy( 0 ),
    doNotReject( doNotReject ),
    resample_iter ( - log( mt_genrand_real1()) / lambda ),
    old_data( NULL )
{
}

HMC::~HMC()
{
    if(old_data){delete old_data; } 
}

void HMC::setSeed( unsigned long seed )
{
    mt_init_genrand( seed );
}

void HMC::iterate(size_t current_iteration,
                  double lastObjective,
                  const Metric & metric,
                  Trajectory & trajectory,
                  MatX & momentum )
{
    
    debug_status( TAG, "iterate", "start" );

    //return if there is something to do for this iteration
    if(current_iteration == resample_iter){
        
        debug_status( TAG, "iterate", "resampling" );
        
        if( doNotReject || 
            !checkForRejection(trajectory, momentum, lastObjective ) )
        {
            getRandomMomentum( metric, current_iteration, momentum );
        }
        
        resample_iter =  current_iteration + 1 
                       - log( mt_genrand_real1()) / lambda;
    }
    
    debug_status( TAG, "iterate", "end" );
}

void HMC::getRandomMomentum(const Metric & metric,
                            size_t current_iteration,
                            MatX & momentum )
{
    debug_status( TAG, "getRandomMomentum", "start" );
    
    const double hmc_alpha = 2/lambda *
                             exp( lambda * current_iteration );

    //this is the standard deviation of the gaussian distribution.
    const double sigma = 1.0/sqrt(hmc_alpha); 
    
    metric.sampleNormalDistribution( sigma, momentum );
    
    debug_status( TAG, "getRandomMomentum", "end" );
    
}


bool HMC::checkForRejection( Trajectory & traj,
                             MatX & momentum,
                             double lastObjective )
{

    debug_status( TAG, "checkForRejection", "start" );
    
    //test the probability of the current state against that of the
    //  old state.
    //start by calculating the current probability
    //the energy of the momentum vector.
    double kinetic_energy = momentum.squaredNorm() * 0.5;

    //the potential energy is the value of the last objective function.
    //  the total energy is below.
    double current_energy = exp(-kinetic_energy) * exp(-lastObjective);

    if ( current_energy < previous_energy ){
        double probability = current_energy / previous_energy;

        //if the probability is too low, 
        //  revert to the previous trajectory
        if ( mt_genrand_real1() > probability ){

            assert( momentum.cols() == old_momentum.cols());
            assert( momentum.rows() == old_momentum.rows());
            
            //perform a swap of the new and old data.
            double * temporary = old_data;
            old_data = traj.getData();
            traj.setData( temporary );
            momentum = old_momentum;
            
            debug_status( TAG, "checkForRejection", "end" );
            
            return true;
        }
    }
    
    //if old_data has not been initialized, initialize it.
    if ( !old_data ){ old_data = new double[ traj.size() ]; }
        
    traj.copyDataTo( old_data );

    previous_energy = current_energy;
    old_momentum = momentum;
    
    debug_status( TAG, "checkForRejection", "end" );
    
    return false;

}

}//namespace
