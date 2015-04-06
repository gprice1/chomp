
#include "HMC.h"
#include "../containers/Trajectory.h"
#include "mzcommon/mersenne.h"
#include "mzcommon/gauss.h"



namespace chomp {


void HMC::setSeed( unsigned long seed ){
    mt_init_genrand( seed );
}

void HMC::iteration(size_t current_iteration, Trajectory & traj,
                    MatX & momentum,
                    const Metric & metric, 
                    double lastObjective ){
    
    //return if there is nothing to do for this iteration
    if(current_iteration != resample_iter){ return; }


    std::cout << "Resampling momentum" << std::endl;

    if( doNotReject || !checkForRejection(traj, momentum, lastObjective ) ){
      
        getRandomMomentum( momentum, current_iteration );
        
        //TODO should this be here?
        metric.multiplyLowerInverseTranspose( momentum );

        traj.update( momentum );
    }

    resample_iter =  current_iteration + 1 
                   - log( mt_genrand_real1()) / lambda;
    std::cout << "Resampled momentum, next iteration: " 
              << resample_iter <<std::endl;
}

void HMC::getRandomMomentum(MatX & momentum, size_t current_iteration)
{
    const double hmc_alpha = 2/lambda *
                             exp( lambda * current_iteration );

    //this is the standard deviation of the gaussian distribution.
    const double sigma = 1.0/sqrt(hmc_alpha); 
    
    //std::cout << "\n\nPrevious Momentum: \n" << momentum;

    //get random univariate gaussians to start.
    for ( int i = 0; i < momentum.size(); i ++ ){
        //since we are using the leapfrog method, we divide by 2
        momentum(i) = gauss_ziggurat( sigma );
    }
}


bool HMC::checkForRejection( Trajectory & traj,
                             MatX & momentum,
                             double lastObjective ) {

    //test the probability of the current state against that of the
    //  old state.
    //start by calculating the current probability
    //the energy of the momentum vector.
    double kinetic_energy = momentum.squaredNorm() * 0.5;

    //the potential energy is the value of the last objective function.
    //  the total energy is below.
    double current_energy = exp(-kinetic_energy) * exp(-lastObjective);

    std::cout << "Current Energy: " << current_energy << "\n"; 
    std::cout << "Previous Energy: " << previous_energy << "\n";

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
            
            return true;
        }
    }
    
    //if old_data has not been initialized, initialize it.
    if ( !old_data ){ old_data = new double[ traj.size() ]; }
        
    traj.copyDataTo( old_data );

    previous_energy = current_energy;
    old_momentum = momentum;
    return false;

}

void HMC::setupRun(){

    previous_energy = 0.0;
    resample_iter = - log( mt_genrand_real1()) / lambda;
}

void HMC::setupHMC( ChompObjectiveType objective_type, double chomp_alpha ){
    
    this->alpha = chomp_alpha * 0.5;
    
    //apply the smoothness operator K to each of the columns.
    if (objective_type == MINIMIZE_VELOCITY){
        n_operator = 2;
        smoothing_operator[0] = -1;
        smoothing_operator[1] = 1;
        //fill the last one with a garbage value to make sure that it is
        // never used
        smoothing_operator[2] = HUGE_VAL;

    }else if (objective_type == MINIMIZE_ACCELERATION){
        n_operator = 3;
        smoothing_operator[0] = 1; 
        smoothing_operator[1] = -2; 
        smoothing_operator[2] = 1; 
    }
}


}//namespace
