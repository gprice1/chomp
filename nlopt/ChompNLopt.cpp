
#include "ChompNLopt.h"
#include "ChompGradient.h"
#include "ConstraintFactory.h"


namespace chomp {


const std::string getNLoptReturnString( nlopt::result & result ){
    if ( result == nlopt::SUCCESS ){ return "SUCCESS"; }
    if ( result == nlopt::STOPVAL_REACHED ){ return "STOPVAL_REACHED"; }
    if ( result == nlopt:: FTOL_REACHED ){ return "FTOL_REACHED"; }
    if ( result == nlopt::XTOL_REACHED  ){ return "XTOL_REACHED "; }
    if ( result == nlopt::MAXEVAL_REACHED ){ return "MAXEVAL_REACHED"; }
    if ( result == nlopt::MAXTIME_REACHED ){ return "MAXTIME_REACHED"; }
    return "UNKNOWN_RESULT";
}


ChompNLopt::ChompNLopt(
              ConstraintFactory * factory,
              const MatX& xi_init,
              const MatX& pinit,
              const MatX& pgoal,
              double obstol,
              int max_iter,
              int N_max,
              ChompObjectiveType obj_t,
              const MatX & lower_bounds,
              const MatX & upper_bounds) :
    ChompOptimizerBase( factory, xi_init, pinit, pgoal,
                        lower_bounds, upper_bounds, obj_t),
    N_max( N_max ),
    max_iter( max_iter ),
    optimizer( NULL )
{

    //currently using the slsqp algorithm for optimization.
    //LD_VAR2
    //SLSQP
    //MMA
    algorithm = nlopt::LD_MMA;
}


ChompNLopt::~ChompNLopt()
{
    if ( optimizer ){delete optimizer; }
}

void ChompNLopt::solve(bool global, bool local)
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


void ChompNLopt::prepareNLoptConstraints(){

    //if there is no factory, do nothing
    if ( !factory ){ return; }

    //clear the old constraints and then get the new constraints.
    factory->getAll( N );

    //Get the dimensionality of the constraint, and fill the
    //  constraint_tolerances vector with the appropriate number of
    //  elements.
    int constraint_dims = factory->numOutput();

    if ( constraint_dims != 0 ){ 
        constraint_tolerances.resize( constraint_dims, 1e-5 );
        
        //the algorithm must be one that can handle equality constraints,
        //  if it is not one that can handle equality constraints,
        //  use the AUGLAG algorithm, with the specified local optimizer.
        if ( algorithm != nlopt::LD_SLSQP ){

            //create the new optimizer, and set the local optimizer.
            nlopt::opt * new_optimizer = new nlopt::opt( nlopt::AUGLAG, 
                                                         xi.size() );
            new_optimizer->set_local_optimizer( *optimizer );

            //delete the old optimizer, and set the optimizer variable
            //  to the new optimizer.
            delete optimizer; 
            optimizer = new_optimizer;
            if ( obstol != 0 ){ optimizer->set_ftol_rel( obstol ); }
            if ( max_iter != 0 ){ optimizer->set_maxeval( max_iter ); }

        }

        //pass the constraint function into nlopt.
        optimizer->add_equality_mconstraint(
                                ConstraintFactory::NLoptConstraint,
                                factory, constraint_tolerances);
    }

}


double ChompNLopt::optimize(){
    
    double previous_objective_value = objective_value;
    notify(CHOMP_INIT, 0, objective_value, -1, 0);
    
    //create the optimizer
    assert( xi.size() == N * M );
    optimizer = new nlopt::opt( algorithm, xi.size() );
    if ( obstol != 0 ){ optimizer->set_ftol_rel( obstol ); }
    if ( max_iter != 0 ){ optimizer->set_maxeval( max_iter ); }

    //prepare the gradient for the run.
    gradient->prepareRun( N );

    //prepare the constraints for optimization.
    //  This MUST be called before set_min_objective and giveBoundsToNLopt,
    //  because prepareNLoptConstraints can change the optimization
    //  routine.
    if ( factory ){  prepareNLoptConstraints(); }
    giveBoundsToNLopt();

    //set the objective function and the termination conditions.
    optimizer->set_min_objective( ChompGradient::NLoptFunction, gradient );
    
        
    //call the optimization routine, get the result and the value
    //  of the objective function.
    std::vector<double> trajectory;
    matToVec( xi, trajectory );

    try{
        result = optimizer->optimize( trajectory , objective_value );
    }catch( std::exception & e ){
        std::cout << "Caught exception: " << e.what() << std::endl;
    }

    vecToMat( trajectory, xi );

    //clean up by deleting the optimizer.
    delete optimizer;
    optimizer = NULL;
    
  
    //notify the observer of the happenings.
    notify( CHOMP_FINISH, 0, objective_value, 
            previous_objective_value, 0);
    std::cout << "Finished with exit code: "
              << getNLoptReturnString(result) << "\n";

    //return the final value of the objective function.
    return objective_value;
}



void ChompNLopt::copyNRows( const MatX & original_bounds, 
                            std::vector<double> & result)
{
    result.resize(N*M);
    
    //eigen matrices are stored in column major format
    for ( int i = 0; i < M; i ++ ){
        for ( int  j = 0; j < N; j ++ ){
            result[i*N + j] = original_bounds(i);
        }
    }
}

void ChompNLopt::giveBoundsToNLopt()
{
    
    assert( optimizer != NULL );

    //set the lower bounds if the lower vector is 
    //  of the correct size.
    if ( lower_bounds.size() == size_t( M ) ){
        std::vector<double> nlopt_lower_bounds;
        copyNRows( lower_bounds, nlopt_lower_bounds );
        
        optimizer->set_lower_bounds( nlopt_lower_bounds );

    }

    //set the upper bounds if the upper matrix is of the 
    //  correct size.
    if ( upper_bounds.size() == size_t( M ) ){

        std::vector<double> nlopt_upper_bounds;
        copyNRows( upper_bounds , nlopt_upper_bounds );
        
        optimizer->set_upper_bounds( nlopt_upper_bounds );
    }
}


}//namespace


