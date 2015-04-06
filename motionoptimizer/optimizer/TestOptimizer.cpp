
#include "TestOptimizer.h"
#include "../containers/ChompGradient.h"
#include "../containers/Trajectory.h"
#include "HMC.h"

namespace chomp {

const char* TestOptimizer::TAG = "TestOptimizer";

TestOptimizer::TestOptimizer( ProblemDescription & problem,
                                         ChompObserver * observer,
                                         double obstol,
                                         double timeout_seconds,
                                         size_t max_iter) : 
    OptimizerBase( problem, observer,
                   obstol, timeout_seconds, max_iter),
    alpha( 0.1 ),
    min_iter( 0 ),
    x_data( NULL ),
    g_data( NULL ),
    H_data( NULL ),
    h_data( NULL ),
    x(NULL, 0, 0),
    g(NULL, 0, 0),
    H(NULL, 0, 0),
    h(NULL, 0, 0)
{
    if ( timeout_seconds <= 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }

    debug_status( TAG, "construction", "end" );
}

void TestOptimizer::solve(){

    debug_status( TAG, "solve", "start" );
    
    
    if (notify(CHOMP_INIT)) { return; }

    if ( timeout_seconds <= 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    x_data = new double [ problem.size() ];
    new (&x) MatMap( g_data, problem.N(), problem.M() );
    problem.copyTrajectoryTo( x_data );

    g_data = new double [ problem.size() ];
    new (&g) MatMap( g_data, problem.N(), problem.M() );

    if (problem.isConstrained() ){
        h_data = new double [ problem.getConstraintDims() ];
        new (&h) MatMap( h_data, problem.getConstraintDims(), 1 );
        
        H_data = new double [ problem.size() * problem.getConstraintDims()];
        new (&H) MatMap( h_data, problem.size(), 
                         problem.getConstraintDims() );
    }
        
    bool not_finished = true;

    debug_status( TAG, "solve", "before evaluateObjective" );

    //Do the optimization
    //TODO find out why prepareIter is not needed anymore.
    //prepareIter( trajectory );

    debug_status( TAG, "solve", "evaluated_objective" );
    
    while (not_finished) { not_finished = iterate(); }
    
    notify(CHOMP_FINISH);

    problem.copyToTrajectory( x_data );

    debug_status( TAG, "solve", "end" );
} 

void TestOptimizer::optimize()
{

  
    debug_status( TAG, "optimize", "start" );
    const Metric & metric = problem.getMetric();
    
    
    debug_status( TAG, "optimize", "after getting constraints" );

    //If there are no constraints,
    //  run the update without constraints.
    if (!problem.isConstrained()) {
        debug_status( TAG, "optimize" , "start unconstrained" );
        
        if (!problem.isCovariant() ){
            metric.solve( g );
        }
        
        x -= g * alpha;        

        debug_status( TAG, "optimize" , "end unconstrained" );
        
    //chomp update with constraints
    } else {
        debug_status( TAG, "optimize", "start constrained" );
      
        const int M = problem.M();
        const int N = problem.N();
        
        if (problem.isCovariant() ){
            
            P = ( H.transpose()*H ).inverse();
            
            debug_status( TAG, "optimize", "got P_inv" );

            const int n_by_m = problem.size();
            
            delta = (alpha*(MatX::Identity(n_by_m, n_by_m)
                     - H*P*H.transpose() )
                    * MatMap( g_data, n_by_m , 1 )
                    + H*P*h).transpose();

        }else {
            P = H;
            
            metric.solve( MatMap( P.data(), N, M * P.cols() ) );

            //debug << "H = \n" << H << "\n";
            //debug << "P = \n" << P << "\n";
          
            HP = H.transpose()*P;

            cholSolver.compute(HP);
            Y = cholSolver.solve(P.transpose());

            //debug << "HP*Y = \n" << HP*Y << "\n";
            //debug << "P.transpose() = \n" << P.transpose() << "\n";
            debug_assert( P.transpose().isApprox( HP*Y ));

            int newsize = H.rows();
            
            assert(newsize == N * M);
            assert(g.rows() == N && g.cols() == M);
            
            ConstMatMap g_flat(g.data(), newsize, 1);
            W = (MatX::Identity(newsize,newsize) - H * Y)
                * g_flat * alpha;

            metric.solve( MatMap( W.data(), N, M * W.cols() ) );

            Y = cholSolver.solve(h);

            delta = W + P * Y;

            assert(delta.rows() == newsize && delta.cols() == 1);
            
        }
        x -= MatMap(delta.data(), N, M);
        debug_status( TAG, "optimize", "end constraint" );
        
    }

    debug_status( TAG, "optimize", "end" );
}


bool TestOptimizer::iterate()
{
    debug_status( TAG, "iterate", "start" );
    

    last_objective = current_objective;
    current_objective = problem.evaluateGradient( x_data, g_data );
    if (problem.isConstrained() ){
        constraint_magnitude = problem.evaluateConstraint( x_data,
                                                           h_data,
                                                           H_data);
    }
    
    //perform optimization
    optimize();

    //check and correct the bounds 
    checkBounds();
    
    //increment the iteration.
    current_iteration ++;

    //check whether optimization is completed.

    debug_status( TAG, "iterate", "end" );
    
    return checkFinished( event );
}

bool TestOptimizer::checkFinished(ChompEventType event)
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
bool TestOptimizer::goodEnough(double oldObjective,
                                    double newObjective )
{
    return (fabs((oldObjective-newObjective)/newObjective)<obstol);
}

//since this is a templated function 
void TestOptimizer::checkBounds()
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
        violations.resize( problem.N(), problem.M());

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
                    violations(i,j) = magnitude; 
                    
                    //save the max magnitude, and its index.
                    if ( -magnitude > max_violation ){
                        max_violation = -magnitude;
                        max_index.first = i;
                        max_index.second = j;
                    }
                //is there a violation of an upper bound?
                }else if ( check_upper && value > upper ){
                    const double magnitude = value - upper; 
                    violations(i,j) = magnitude; 

                    //save the max magnitude, and its index.
                    if ( magnitude > max_violation ){
                        max_violation = magnitude;
                        max_index.first = i;
                        max_index.second = j;
                    }
                }else {
                    violations(i,j) = 0;
                }
            }
        }
        

        //There are bounds violations if the largest
        //violation does not have a magnitude of zero.
        violation = ( max_violation > 0.0 );
        
        //TODO, maybe only check bounds violations on
        //  non-subsampled trajectories
        if( violation ){
            problem.getMetric().solve(violations);

            //scale the bounds_violation matrix so that it sets the
            //  largest violation to zero.
            double current_mag = violations( max_index.first,
                                             max_index.second);
            const double scale = (current_mag > 0 ?
                                  max_violation/current_mag :
                                 -max_violation/current_mag  );
            problem.updateTrajectory( violations * scale );
        }
    }


    debug_status( TAG, "checkBounds", "end" );
}

}//namespace
