
#include "ChompOptimizerBase.h"
#include "ChompGradient.h"
#include "Constraint.h"

namespace chomp {

ChompOptimizerBase::ChompOptimizerBase( ConstraintFactory * f,
                                        const MatX & xinit,
                                        const MatX & pinit,
                                        const MatX & pgoal,
                                        const MatX & lower_bounds,
                                        const MatX & upper_bounds,
                                        ChompObjectiveType object_type,
                                        double total_time ) :
    factory( f ), 
    observer( NULL ),
    objective_type( object_type ),
    N(xinit.rows()), M(xinit.cols()),
    xi( xinit ),
    lower_bounds( lower_bounds ), upper_bounds( upper_bounds )
{

    if ( timeout_seconds < 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }

}


bool ChompOptimizerBase::checkFinished(ChompEventType event)
{
    
    //get the new value of the objective
    current_objective = gradient->evaluateObjective();

    //test for termination conditions
    bool greater_than_min = curr_iter > min_iter;
    bool greater_than_max = curr_iter > max_iter;
    bool converged = goodEnough( last_objective, current_objective );
    bool observer_flag = notify(event, curr_iter,
                                current_objective,
                                last_objective,
                                hmag);

    if ( canTimeout && stop_time < TimeStamp::now() ) {
        notify(CHOMP_TIMEOUT, curr_iter, current_objective,
                                         last_objective, hmag);
        return false;
    }

    if (greater_than_max || ( greater_than_min && converged ) ||
        observer_flag )
    {
        return false;
    } 

    lastObjective = curObjective;

    debug << "Ending Iteration" << std::endl;
    return true;
}

// returns true if performance has converged
bool ChompOptimizerBase::goodEnough(double oldObjective,
                                    double newObjective )
{
    return (fabs((oldObjective-newObjective)/newObjective)<objRelErrTol);
}

//since this is a templated function 
void ChompOptimizerBase::checkBounds( )
{

    const bool check_upper = (upper_bounds.size() == M);
    const bool check_lower = (lower_bounds.size() == M);

    //if there are bounds to check, check for the violations.

    if ( check_upper || check_lower ){
        
        bool violation;
        bounds_violations.resize( xi.rows(), xi.cols() );

        const int max_checks = 10;
        int count = 0;

        do{
            count ++;
            double max_violation = 0.0;
            std::pair< int, int > max_index;

            for ( int j = 0; j < xi.cols(); j ++ ) {
                const double upper = (check_upper ? upper_bounds(j) : 0 );
                const double lower = (check_upper ? lower_bounds(j) : 0 );
                
                for( int i = 0; i < xi.rows(); i ++ ){

                    //is there a violation of a lower bound?
                    if ( check_lower && xi(i,j) < lower ){
                        const double magnitude = xi(i,j) - lower;
                        bounds_violations(i,j) = magnitude; 
                        
                        //save the max magnitude, and its index.
                        if ( -magnitude > max_violation ){
                            max_violation = -magnitude;
                            max_index.first = i;
                            max_index.second = j;
                        }
                    //is there a violation of an upper bound?
                    }else if ( check_upper && xi(i,j) > upper ){
                        const double magnitude = xi(i,j) - upper; 
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
                if ( bounds_violations.rows() == xi.rows() ){
                    skylineCholSolve( gradient->L, bounds_violations );
                }else {
                    assert( gradient->L_sub.rows() == xi.rows() );
                    skylineCholSolve( gradient->L_sub, bounds_violations );
                }

                //scale the bounds_violation matrix so that it sets the
                //  largest violation to zero.
                double current_mag = bounds_violations( max_index.first,
                                                        max_index.second );
                const double scale = (current_mag > 0 ?
                                      max_violation/current_mag :
                                     -max_violation/current_mag  );
                const_cast<Eigen::MatrixBase<Derived>&>(xi) -= 
                                             bounds_violations * scale;
            }
        //continue to check and fix limit violations as long as they exist.
        }while( violation && count < max_checks );
    }
}

}//namespace
