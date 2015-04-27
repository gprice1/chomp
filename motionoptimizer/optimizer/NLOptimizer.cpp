
#include "NLOptimizer.h"

namespace mopt {

const char * NLOptimizer::TAG = "NLOptimizer";

const std::string getNLoptReturnString( nlopt::result & result ){
    if ( result == nlopt::SUCCESS ){ return "SUCCESS"; }
    if ( result == nlopt::STOPVAL_REACHED ){ return "STOPVAL_REACHED"; }
    if ( result == nlopt:: FTOL_REACHED ){ return "FTOL_REACHED"; }
    if ( result == nlopt::XTOL_REACHED  ){ return "XTOL_REACHED "; }
    if ( result == nlopt::MAXEVAL_REACHED ){ return "MAXEVAL_REACHED"; }
    if ( result == nlopt::MAXTIME_REACHED ){ return "MAXTIME_REACHED"; }
    return "UNKNOWN_RESULT";
}


NLOptimizer::NLOptimizer( ProblemDescription & problem,
                          Observer * observer,
                          double obstol,
                          double timeout_seconds,
                          size_t max_iter) : 

    OptimizerBase( problem, observer,
                   obstol, timeout_seconds,
                   max_iter),
    algorithm( nlopt::LD_MMA )
{
}


NLOptimizer::~NLOptimizer(){}

void NLOptimizer::solve()
{
    debug_status( TAG, "solve", "start" );

    //create the optimizer
    nlopt::opt optimizer( algorithm, problem.size() );

    //set termination conditions.
    if ( timeout_seconds > 0){optimizer.set_maxtime(timeout_seconds);}
    if ( obstol > 0 ){optimizer.set_ftol_rel( obstol ); }
    if ( max_iter > 0 ){ optimizer.set_maxeval( max_iter ); }

    //prepare the constraints for optimization.
    //  This MUST be called before set_min_objective and giveBoundsToNLopt,
    //  because prepareNLoptConstraints can change the optimization
    //  routine.
    if ( problem.isConstrained() ){ 
        prepareConstraints( optimizer );
    }

    if ( problem.isCollisionConstraint() ){
        prepareCollisionConstraint( optimizer );
    }

    prepareBounds( optimizer );
    
    //set the objective function and the termination conditions.
    optimizer.set_min_objective( objectiveFunction, this );
    
    debug_status( TAG, "solve", "before copying data" );
    
    //call the optimization routine, get the result and the value
    //  of the objective function.
    std::vector<double> optimization_vector;
    problem.copyTrajectoryTo( optimization_vector );
    
    debug_status( TAG, "solve", "pre-optimize" );
    
    //many of the optimizer algorithms fail, but still return very good
    //  trajectories, so catch the failures and return the trajectory,
    //  even if failure occurs.
    try{
        double objective_value;
        result = optimizer.optimize(optimization_vector, objective_value);
        current_objective = objective_value;
    }catch( std::exception & e ){
        std::cout << "Caught exception: " << e.what() << std::endl;
    }
    
    debug_status( TAG, "solve", "post-optimize" );
    
    problem.copyToTrajectory( optimization_vector );

    //notify the observer of the happenings.
    debug << "Finished with exit code: "
              << getNLoptReturnString(result) << "\n";

    debug_status( TAG, "solve", "end" );

}

void NLOptimizer::prepareCollisionConstraint( nlopt::opt & optimizer )
{
    
    const nlopt::algorithm current_algorithm = optimizer.get_algorithm();

    if ( current_algorithm != nlopt::LD_MMA || 
         current_algorithm != nlopt::AUGLAG ||
         current_algorithm != nlopt::AUGLAG_EQ ||
         current_algorithm != nlopt::LD_SLSQP ||
         current_algorithm != nlopt::LN_COBYLA )
    {
        nlopt::opt new_optimizer( nlopt::AUGLAG, problem.size());
        new_optimizer.set_local_optimizer( optimizer );

        optimizer = new_optimizer;
    }
    
    optimizer.add_inequality_constraint( collisionConstraintFunction,                                                 this, 1e-4 );
}

void NLOptimizer::prepareConstraints( nlopt::opt & optimizer )
{

    //Get the dimensionality of the constraint, and fill the
    //  constraint_tolerances vector with the appropriate number of
    //  elements.
    int constraint_dims = problem.getConstraintDims() ;

    if ( constraint_dims != 0 ){ 
        
        //the algorithm must be one that can handle equality constraints,
        //  if it is not one that can handle equality constraints,
        //  use the AUGLAG algorithm, with the specified local optimizer.
        if ( algorithm != nlopt::LD_SLSQP && 
             algorithm != nlopt::LN_COBYLA ){

            //create the new optimizer, and set the local optimizer.
            nlopt::algorithm algorithm2;

            if ( algorithm == nlopt::LD_MMA  ){
                algorithm2 = nlopt::AUGLAG_EQ;
            }else {
                algorithm2 = nlopt::AUGLAG;
            }
            
            nlopt::opt new_optimizer( algorithm2, problem.size());
            new_optimizer.set_local_optimizer( optimizer );

            optimizer = new_optimizer;
        }
        
        std::vector<double> constraint_tolerances;
        constraint_tolerances.resize( constraint_dims, 1e-5 );
        
        //pass the constraint function into nlopt.
        optimizer.add_equality_mconstraint( constraintFunction,
                                            this, constraint_tolerances);
    }
}


void NLOptimizer::prepareBounds( nlopt::opt & optimizer )
{
    
    //TODO throw error because the bounds are the wrong size.
    
    //set the lower bounds if the lower vector is 
    //  of the correct size.
    if ( problem.isBounded() ){
        std::vector<double> lower_bounds;
        std::vector<double> upper_bounds;
        problem.getFullBounds( lower_bounds, upper_bounds );
        
        if ( problem.getLowerBounds().size() == problem.M() )
        {
            optimizer.set_lower_bounds( lower_bounds );
        }
        
        //set the upper bounds if the upper matrix is of the 
        //  correct size.
        if ( problem.getUpperBounds().size() == problem.M() )
        {
            optimizer.set_upper_bounds( upper_bounds );
        }
    }
}


//a wrapper function for passing the ChompGradient to NLopt.
double NLOptimizer::objectiveFunction(unsigned n,
                                      const double * x,
                                      double* grad,
                                      void *data)
{
    NLOptimizer * opt = reinterpret_cast<NLOptimizer*>(data);

    opt->last_objective = opt->current_objective;
    opt->current_objective = opt->problem.evaluateObjective( x, grad );
    
    //TODO - report the constraint violations
    opt->notify( event );
    
    opt->current_iteration ++;

    return opt->current_objective;
}

//a wrapper function for passing the ChompGradient to NLopt.
void NLOptimizer::constraintFunction(unsigned constraint_dim,
                                       double * h,
                                       unsigned n_by_m,
                                       const double * x,
                                       double* H,
                                       void *data)
{
    NLOptimizer * opt = reinterpret_cast<NLOptimizer*>(data);

    opt->constraint_magnitude = opt->problem.evaluateConstraint( x, h, H);
}

//a wrapper function for passing the ChompGradient to NLopt.
double NLOptimizer::collisionConstraintFunction(unsigned n_by_m,
                                              const double * x,
                                              double* h,
                                              void *data)
{
    NLOptimizer * opt = reinterpret_cast<NLOptimizer*>(data);

     return opt->problem.evaluateCollisionFunction( x, h);
}



}//namespace


