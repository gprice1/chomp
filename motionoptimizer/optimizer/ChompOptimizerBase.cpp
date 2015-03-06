
#include "ChompOptimizerBase.h"
#include "ChompGradient.h"
#include "Trajectory.h"
#include "HMC.h"

namespace chomp {

    
ChompOptimizerBase::ChompOptimizerBase( Trajectory & traj,
                                        ConstraintFactory * factory,
                                        ChompGradient * gradient,
                                        ChompObserver * observer,
                                        double obstol,
                                        double timeout_seconds,
                                        size_t max_iter,
                                        const MatX & lower_bounds,
                                        const MatX & upper_bounds) : 
    OptimizerBase( traj, factory, gradient, observer,
                   obstol, timeout_seconds, max_iter,
                   lower_bounds, upper_bounds ),
    alpha( 0.1 ),
    hmag( 0 ), 
    last_objective( HUGE_VAL ),
    current_objective( HUGE_VAL ),
    curr_iter( 0 ),
    min_iter( 0 ),
    use_momentum( false ),
    hmc( NULL )
{

    if ( timeout_seconds < 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }

}

void ChompOptimizerBase::solve(){

    if (notify(CHOMP_INIT, 0, last_objective, -1, hmag)) { return; }

    if ( timeout_seconds < 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    
    if ( hmc ){ 
        use_momentum = true;
        hmc->setupHMC( gradient->objective_type, alpha );
        hmc->setupRun();
    }

    if ( use_momentum ){
        momentum.resize( trajectory.N(), trajectory.M() );
        momentum.setZero();
    }
    
    bool not_finished = true;

    //Do the optimization
    //TODO find out why prepareIter is not needed anymore.
    //prepareIter( trajectory );
    last_objective = gradient->evaluateObjective( trajectory );
    while (not_finished) { not_finished = iterate(); }
    
    notify(CHOMP_FINISH, 0, last_objective, -1, hmag);
} 

bool ChompOptimizerBase::iterate(){
    
    //perform optimization
    optimize();

    //check and correct the bounds 
    checkBounds();
    
    //increment the iteration.
    curr_iter ++;
    
    //prepare for the next iteration.
    //TODO find out why prepareIter is not needed anymore.
    //prepareIter( trajectory );

    if ( hmc ) {
        //TODO make sure that we are doing Global Chomp,
        //  or add momentum into local chomp
        hmc->iteration( curr_iter, trajectory, momentum,
                        gradient->getInvAMatrix(),
                        last_objective );
    }
    
    gradient->getGradient( trajectory );

    //check whether optimization is completed.
    return checkFinished( event );
}

bool ChompOptimizerBase::checkFinished(ChompEventType event)
{
    
    //get the new value of the objective
    current_objective = gradient->evaluateObjective( trajectory );


    if ( canTimeout && stop_time < TimeStamp::now() ) {
        notify(CHOMP_TIMEOUT, curr_iter, current_objective,
                                         last_objective, hmag);
        return false;
    }


    //test for termination conditions
    bool greater_than_min = curr_iter > min_iter;
    bool greater_than_max = curr_iter > max_iter;
    bool converged = goodEnough( last_objective, current_objective );
    bool observer_flag = notify(event, curr_iter,
                                current_objective,
                                last_objective,
                                hmag);

    if (greater_than_max || ( greater_than_min && converged ) ||
        observer_flag )
    {
        return false;
    } 

    last_objective = current_objective;

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

    const bool check_upper = (upper_bounds.size() == trajectory.M());
    const bool check_lower = (lower_bounds.size() == trajectory.M());

    //if there are bounds to check, check for the violations.

    if ( check_upper || check_lower ){
        
        bool violation;
        bounds_violations.resize( trajectory.rows(), trajectory.cols() );

        const int max_checks = 10;
        int count = 0;

        do{
            count ++;
            double max_violation = 0.0;
            std::pair< int, int > max_index;

            for ( int j = 0; j < trajectory.cols(); j ++ ) {
                const double upper = (check_upper ? upper_bounds(j) : 0 );
                const double lower = (check_upper ? lower_bounds(j) : 0 );
                
                for( int i = 0; i < trajectory.rows(); i ++ ){

                    //is there a violation of a lower bound?
                    if ( check_lower && trajectory(i,j) < lower ){
                        const double magnitude = trajectory(i,j) - lower;
                        bounds_violations(i,j) = magnitude; 
                        
                        //save the max magnitude, and its index.
                        if ( -magnitude > max_violation ){
                            max_violation = -magnitude;
                            max_index.first = i;
                            max_index.second = j;
                        }
                    //is there a violation of an upper bound?
                    }else if ( check_upper && trajectory(i,j) > upper ){
                        const double magnitude = trajectory(i,j) - upper; 
                        bounds_violations(i,j) = magnitude; 

                        //save the max magnitude, and its index.
                        if ( magnitude > max_violation ){
                            max_violation = magnitude;
                            max_index.first = i;
                            max_index.second = j;
                        }
                    }else {
                        bounds_violations(i,j) = 0;
                    }
                }
            }
            

            //There are violations in the largest violation does not
            //  have a magnitude of zero.
            violation = ( max_violation > 0.0 );
            
            if( violation ){
                //smooth out the bounds violation matrix. With
                //  the appropriate smoothing matrix.
                if ( bounds_violations.rows() == trajectory.rows() ){
                    skylineCholSolve( gradient->L, bounds_violations );
                }else {
                    assert( gradient->L_sub.rows() == trajectory.rows() );
                    skylineCholSolve( gradient->L_sub, bounds_violations );
                }

                //scale the bounds_violation matrix so that it sets the
                //  largest violation to zero.
                double current_mag = bounds_violations( max_index.first,
                                                        max_index.second );
                const double scale = (current_mag > 0 ?
                                      max_violation/current_mag :
                                     -max_violation/current_mag  );
                trajectory.update( bounds_violations * scale );
            }
        //continue to check and fix limit violations as long as they exist.
        }while( violation && count < max_checks );
    }
}

}//namespace
