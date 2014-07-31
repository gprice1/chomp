
#include "ChompNLopt.h"



namespace chomp {

ChompNLopt::ChompNLopt(
              const MatX& xi_init,
              const MatX& pinit,
              const MatX& pgoal,
              double obstol,
              int max_iter,
              int N_max,
              const std::vector<double> & lower,
              const std::vector<double> & upper) :
    N( xi_init.rows() ), M(xi_init.cols()),
    N_max( N_max ),
    max_iter( max_iter ),
    xi( xi_init ),
    lower( lower ), upper( upper )

{
    assert( M == pinit.size());
    assert( M == pgoal.size());
    
    //currently using the slsqp algorithm for optimization.
    algorithm = nlopt::LD_SLSQP;

    gradient = new ChompGradient( pinit, pgoal );

}



ChompNLopt::~ChompNLopt()
{
    if ( optimizer ){delete optimizer; }
    if ( gradient ){delete gradient; }
}

void ChompNLopt::solve()
{
    
    //optimize at the current resolution
    optimize();

    //if the current resolution is not the final resolution,
    //  upsample and then optimize. Repeat until the final resolution
    //  is reached or exceeded.
    while ( N < N_max ){

        //upsample the trajectory and prepare for the next
        //  stage of optimization.
        MatX xi_up;
        upsampleTrajectory( xi, gradient->q0,
                                gradient->q1,
                                gradient->dt,
                                gradient->objective_type, xi_up );
        xi = xi_up;
        N = xi.rows();

        optimize();
    }
}

double ChompNLopt::optimize(){
    
    //create the optimizer
    optimizer = new nlopt::opt( algorithm, xi.size() );
    
    //If the lower and upper bounds vectors are of the correct size,
    //  set the upper and lower bounds.
    if (lower.size() == size_t(M)){ setLowerBounds( lower ); }
    if (upper.size() == size_t(M)){ setUpperBounds( upper ); }
    
    //set the objective function and the termination conditions.
    optimizer->set_min_objective( ChompGradient::NLoptFunction, gradient );
    if ( obstol != 0 ){ optimizer->set_ftol_rel( obstol ); }
    if ( max_iter != 0 ){ optimizer->set_maxeval( max_iter ); }
    
    //prepare the gradient for the run.
    gradient->prepareRun( N );
    
    //call the optimization routine, get the result and the value
    //  of the objective function.
    std::vector<double> trajectory;
    matToVec( xi, trajectory );
    result = optimizer->optimize( trajectory , objective_value );
    vecToMat( trajectory, xi );

    //clean up by deleting the optimizer.
    delete optimizer;
    
    //return the final value of the objective function.
    return objective_value;
}



void ChompNLopt::copyNRows( const std::vector<double> & original, 
                            std::vector<double> & result)
{
    result.resize(N*M);
    
    //eigen matrices are stored in column major format
    for ( int i = 0; i < M; i ++ ){
        for ( int  j = 0; j < N; j ++ ){
            result[i*N + j] = original[i];
        }
    }
}

void ChompNLopt::setLowerBounds( const std::vector<double> & lower)
{
    assert( lower.size() == size_t(M) );

    std::vector<double> lower_bounds;
    copyNRows( lower, lower_bounds );
    
    optimizer->set_lower_bounds( lower_bounds );
}

void ChompNLopt::setUpperBounds( const std::vector<double> & upper)
{
    assert( upper.size() == size_t(M) );

    std::vector<double> upper_bounds;
    copyNRows( upper, upper_bounds );
    
    optimizer->set_upper_bounds( upper_bounds );
}


}//namespace


