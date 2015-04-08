
#include "ChompOptimizerBase.h"
#include "../containers/ChompGradient.h"
#include "../containers/Trajectory.h"
#include "HMC.h"

namespace chomp {

const char* ChompOptimizerBase::TAG = "ChompOptimizerBase";

ChompOptimizerBase::ChompOptimizerBase( ProblemDescription & problem,
                                         ChompObserver * observer,
                                         double obstol,
                                         double timeout_seconds,
                                         size_t max_iter) : 
    OptimizerBase( problem, observer,
                   obstol, timeout_seconds, max_iter),
    alpha( 0.1 ),
    min_iter( 0 ),
    use_momentum( false ),
    hmc( NULL )
{
    if ( timeout_seconds <= 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }

    debug_status( TAG, "construction", "end" );
}

void ChompOptimizerBase::solve(){

    debug_status( TAG, "solve", "start" );
    
    
    if ( timeout_seconds <= 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    
    if ( hmc ){ 
        use_momentum = true;
        hmc->setupHMC( problem.getTrajectory().getObjectiveType(), alpha );
        hmc->setupRun();
    }

    if ( use_momentum ){
        momentum.resize( problem.N(), problem.M() );
        momentum.setZero();
    }
    
    g.resize( problem.N(), problem.M() );
    
    bool not_finished = true;

    debug_status( TAG, "solve", "before evaluateObjective" );

    //Do the optimization
    //TODO find out why prepareIter is not needed anymore.
    //prepareIter( trajectory );

    debug_status( TAG, "solve", "evaluated_objective" );
    
    while (not_finished) { not_finished = iterate(); }
    
    debug_status( TAG, "solve", "end" );
} 

bool ChompOptimizerBase::iterate(){
    debug_status( TAG, "iterate", "start" );
    
    last_objective = current_objective;
    current_objective = problem.evaluateGradient( g );
    
    //perform optimization
    optimize();

    //check and correct the bounds 
    checkBounds();
    
    //increment the iteration.
    current_iteration ++;

    //TODO add back in the HMC once everything is figured out.
    /*
    if ( hmc ) {
        //TODO make sure that we are doing Global Chomp,
        //  or add momentum into local chomp
        hmc->iteration( curr_iter, problemtrajectory, momentum,
                        problem.getLMatrix(),
                        last_objective );
    }
    */

    //check whether optimization is completed.

    debug_status( TAG, "iterate", "end" );
    
    return checkFinished( event );
}

bool ChompOptimizerBase::checkFinished(ChompEventType event)
{
    
    //get the new value of the objective


    if ( canTimeout && stop_time < TimeStamp::now() ) {
        notify(CHOMP_TIMEOUT);
        return false;
    }


    //test for termination conditions
    bool greater_than_min = current_iteration > min_iter;
    bool greater_than_max = current_iteration > max_iter;
    bool converged = goodEnough( last_objective, current_objective );
    bool observer_flag = notify(event);

    if (greater_than_max || ( greater_than_min && converged ) ||
        observer_flag )
    {
        return false;
    } 

    return true;
}

// returns true if performance has converged
bool ChompOptimizerBase::goodEnough(double oldObjective,
                                    double newObjective )
{
    return (fabs((oldObjective-newObjective)/newObjective)<obstol);
}

//since this is a templated function 
void ChompOptimizerBase::checkBounds()
{

    debug_status( TAG, "checkBounds", "start" );
    
    const MatX & lower_bound = problem.getLowerBounds();
    const MatX & upper_bound = problem.getUpperBounds();

    const bool check_lower = (lower_bound.size() == problem.M());
    const bool check_upper = (upper_bound.size() == problem.M());

    //if there are bounds to check, check for the violations.
    if ( !check_upper || !check_lower ){ return; }

    bool violation = true;

    //Terminate if it does more than 10 iterations of bounds 
    //  correction (an arbitrary number), or if there are no
    //  more bounds violations
    for ( int count = 0; count < 10 && violation; count ++ ){
        
        double max_violation = 0.0;
        std::pair< int, int > max_index;

        //iterate over all of the points, and make a matrix which
        //  has the magnitude of all the bounds violations
        //  for each variable in the trajectory
        for ( int j = 0; j < problem.N(); j ++ ) {
            const double upper = (check_upper ? upper_bound(j) : HUGE_VAL);
            const double lower = (check_lower ? lower_bound(j) :-HUGE_VAL);
            
            for( int i = 0; i < problem.N(); i ++ ){

                //is there a violation of a lower bound?
                const double value = problem.getTrajectory()(i,j);
                if ( check_lower && value < lower )
                {
                    const double magnitude = value - lower;
                    g(i,j) = magnitude; 
                    
                    //save the max magnitude, and its index.
                    if ( -magnitude > max_violation ){
                        max_violation = -magnitude;
                        max_index.first = i;
                        max_index.second = j;
                    }
                //is there a violation of an upper bound?
                }else if ( check_upper && value > upper ){
                    const double magnitude = value - upper; 
                    g(i,j) = magnitude; 

                    //save the max magnitude, and its index.
                    if ( magnitude > max_violation ){
                        max_violation = magnitude;
                        max_index.first = i;
                        max_index.second = j;
                    }
                }else {
                    g(i,j) = 0;
                }
            }
        }
        

        //There are bounds violations if the largest
        //violation does not have a magnitude of zero.
        violation = ( max_violation > 0.0 );
        
        //TODO, maybe only check bounds violations on
        //  non-subsampled trajectories
        if( violation ){
            problem.getMetric().solve(g);

            //scale the bounds_violation matrix so that it sets the
            //  largest violation to zero.
            double current_mag = g( max_index.first, max_index.second);
            const double scale = (current_mag > 0 ?
                                  max_violation/current_mag :
                                 -max_violation/current_mag  );
            problem.updateTrajectory( g * scale );
        }
    }


    debug_status( TAG, "checkBounds", "end" );
}

}//namespace
