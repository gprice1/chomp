
#include "MotionOptimizer.h"

namespace chomp {

nlopt::algorithm getNLoptAlgorithm( OptimizationAlgorithm alg ){

    switch (alg){
    case MMA_NLOPT: return nlopt::LD_MMA;
    case CCSAQ_NLOPT: return nlopt::LD_CCSAQ;
    case SLSQP_NLOPT: return nlopt::LD_SLSQP;
    case LBFGS_NLOPT: return nlopt::LD_LBFGS;
    case TNEWTON_PRECOND_RESTART_NLOPT:
        return nlopt::LD_TNEWTON_PRECOND_RESTART;
    case TNEWTON_RESTART_NLOPT: return nlopt::LD_TNEWTON_RESTART;
    case TNEWTON_NLOPT: return nlopt::LD_TNEWTON;
    case VAR1_NLOPT: return nlopt::LD_VAR1;
    case VAR2_NLOPT: return nlopt::LD_VAR2;

    //TODO Throw error.
    default: return nlopt::LD_MMA;
    }
}

const char* MotionOptimizer::TAG = "MotionOptimizer";

//constructor.
MotionOptimizer::MotionOptimizer(ChompObserver * observer,
                                 double obstol,
                                 double timeout_seconds,
                                 size_t max_iter,
                                 const MatX & lower_bounds,
                                 const MatX & upper_bounds,
                                 OptimizationAlgorithm alg1,
                                 OptimizationAlgorithm alg2,
                                 int N_max) :
    observer( observer ),
    N_max( N_max ),
    full_global_at_final( false ),
    do_subsample( true ),
    obstol( obstol ),
    timeout_seconds( timeout_seconds ),
    alpha( -1 ),
    max_iterations( max_iter ),
    algorithm1( alg1 ),
    algorithm2( alg2 )
{
    debug << "MotionOptimizer initialized" << std::endl;
}

void MotionOptimizer::solve(){
    
    debug_status( TAG, "solve", "start");
    
    N_min = problem.N();

    //optimize at the current resolution
    optimize( getOptimizer(algorithm1) );

    //if the current resolution is not the final resolution,
    //  upsample and then optimize. Repeat until the final resolution
    //  is reached or exceeded.
    while ( problem.N() < N_max ){

        //upsample the trajectory and prepare for the next
        //  stage of optimization.
        problem.upsample();

        if (do_subsample ){
            optimize( getOptimizer( algorithm1 ), true );

            OptimizationAlgorithm alg = (algorithm2 == NONE ?
                                         algorithm1 :
                                         algorithm2);
            optimize( getOptimizer( alg ) );
            
        }else  {
            optimize( getOptimizer( algorithm1 ));
        }
            
    }

    //If full_global_at_final
    if ( full_global_at_final && do_subsample && N_min < problem.N() ){
        optimize( getOptimizer( algorithm1 ));
    }
    
    debug_status( TAG, "solve", "end");
}


void MotionOptimizer::optimize( OptimizerBase * optimizer, 
                                bool subsample ){
    
    //if the optimizer is NULL, do not evaluate it.
    if ( !optimizer ) { return; }

    debug_status( TAG, "optimize", "start");

    if( subsample ){ problem.subsample(); }

    problem.prepareRun();
    
    optimizer->solve();

    delete optimizer;

    if( subsample ){ problem.stopSubsample(); }

    debug_status( TAG, "optimize", "end");

}

OptimizerBase * MotionOptimizer::getOptimizer(OptimizationAlgorithm alg )
{
    //create the optimizer

    //If the algorith is CovariantChomp, or if the Algorithm is a
    //  chomp variant with a covaraint trajectory, do CovariantChomp
    if ( ((alg == GLOBAL_CHOMP || alg == LOCAL_CHOMP) &&
         problem.isCovariantOptimization() ) &&
         !problem.isSubsampled() )
    {
        alg = COVARIANT_CHOMP;
    }

    if ( alg == COVARIANT_CHOMP ){
        problem.doCovariantOptimization();
        ChompCovariantOptimizer * opt = new ChompCovariantOptimizer(
                                  problem,
                                  observer, 
                                  obstol, timeout_seconds,
                                  max_iterations);
        if (alpha > 0){ opt->setAlpha( alpha ); }
        return opt;

    } else if ( alg == GLOBAL_CHOMP ){
        ChompOptimizer * opt = new ChompOptimizer(
                                  problem,
                                  observer, 
                                  obstol, timeout_seconds,
                                  max_iterations);
        if (alpha > 0){ opt->setAlpha( alpha ); }
        return opt;
        
    } else if ( alg == LOCAL_CHOMP ){
        ChompLocalOptimizer * opt = new ChompLocalOptimizer(
                                      problem, observer, 
                                      obstol, timeout_seconds,
                                      max_iterations);
        
        if (alpha > 0){ opt->setAlpha( alpha ); }
        return opt;
        
    } else if ( alg > GLOBAL_CHOMP && alg < NONE){
#ifdef NLOPT_FOUND
        NLOptimizer * opt = new NLOptimizer(
                                  problem, observer, 
                                  obstol, timeout_seconds,
                                  max_iterations);
        
        opt->setAlgorithm( getNLoptAlgorithm( alg ) );
        return opt;
#else 
        std::cerr << "NLopt optimization libraries are"
                  << " not available, please use a different"
                  << " optimization algorithm" << std::endl;
        return NULL; 
#endif

    }
    
    return NULL;
}

void MotionOptimizer::setLowerBounds( const MatX & lower )
{
    problem.setLowerBounds( lower );
}
void MotionOptimizer::setLowerBounds( const std::vector<double> & lower)
{
    problem.setLowerBounds( ConstMatMap(lower.data(), lower.size(), 1 ));
}
void MotionOptimizer::setLowerBounds( const double * lower, int M)
{
    problem.setLowerBounds( ConstMatMap(lower, M, 1) );
}


void MotionOptimizer::setUpperBounds(const MatX & upper )
{
    problem.setUpperBounds( upper );    
}
void MotionOptimizer::setUpperBounds( const std::vector<double> & upper)
{
    problem.setUpperBounds( ConstMatMap(upper.data(), upper.size(), 1 ));
}

void MotionOptimizer::setUpperBounds( const double * upper, int M)
{
    problem.setUpperBounds( ConstMatMap(upper, M, 1) );
}

void MotionOptimizer::setBounds( const MatX & lower, const MatX & upper )
{
    problem.setLowerBounds( lower );
    problem.setUpperBounds( upper );
}
void MotionOptimizer::setBounds( const std::vector<double> & lower, 
                                 const std::vector<double> & upper)
{
    setLowerBounds( lower );
    setUpperBounds( upper );
}
void MotionOptimizer::setBounds( const double* lower, 
                                 const double* upper,
                                 int M)
{
    setLowerBounds( lower, M );
    setUpperBounds( upper, M );
}

void MotionOptimizer::addConstraint( Constraint * c,
                                     double start_time,
                                     double end_time)
{
    problem.factory.addConstraint( c, start_time, end_time );
}


}//namespace
