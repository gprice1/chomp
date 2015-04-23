
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
    
    //get random univariate gaussians to start.
    for ( int i = 0; i < momentum.size(); i ++ ){
        momentum(i) = gauss_ziggurat( sigma );
    }
    
    //transform the univariate gaussians into a multivariate
    //  gaussian via multiplying by the inverse transpose of the
    //  metric.
    metric.multiplyLowerInverseTranspose( momentum );
    
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
