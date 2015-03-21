
#include "MotionOptimizer.h"

namespace chomp {

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
    gradient( trajectory ),
    factory ( trajectory ),
    observer( observer ),
    N_max( N_max ),
    use_goalset( false ),
    full_global_at_final( false ),
    do_subsample( true ),
    obstol( obstol ),
    timeout_seconds( timeout_seconds ),
    max_iterations( max_iter ),
    lower_bounds( lower_bounds ),
    upper_bounds( upper_bounds ),
    algorithm1( alg1 ),
    algorithm2( alg2 )
    
{
    debug << "MotionOptimizer initialized" << std::endl;
}

void MotionOptimizer::solve(){
    
    debug_status( TAG, "solve", "start");
    
    N_min = trajectory.N();

    //optimize at the current resolution
    optimize( getOptimizer(algorithm1) );

    //if the current resolution is not the final resolution,
    //  upsample and then optimize. Repeat until the final resolution
    //  is reached or exceeded.
    while ( trajectory.rows() < N_max ){

        //upsample the trajectory and prepare for the next
        //  stage of optimization.
        trajectory.upsample();

        if (do_subsample) { optimize( getOptimizer( algorithm1 ), true ); }
        optimize( getOptimizer( algorithm2 ) );
    }

    debug_status( TAG, "solve", "end");
}


void MotionOptimizer::optimize( OptimizerBase * optimizer, 
                                bool subsample ){
    
    debug_status( TAG, "optimize", "start");

    if( subsample ){ trajectory.subsample(); }
    
    if ( !factory.empty() ){ factory.getAll( trajectory.N() ); }
    if ( use_goalset ){ prepareGoalSet(); }

    gradient.prepareRun( trajectory, use_goalset );

    optimizer->solve();

    delete optimizer;

    if (use_goalset ){ finishGoalSet(); }
    if( subsample ){ trajectory.endSubsample(); }

    debug_status( TAG, "optimize", "end");

}

OptimizerBase * MotionOptimizer::getOptimizer(OptimizationAlgorithm alg )
{
    //create the optimizer
    //TODO include other optimization schemes

    switch ( alg ){
    case GLOBAL_CHOMP:
        return new ChompOptimizer(trajectory, &factory, 
                                  &gradient, observer, 
                                  obstol, timeout_seconds,
                                  max_iterations, 
                                  lower_bounds, upper_bounds );
    case LOCAL_CHOMP:
        return new ChompLocalOptimizer(trajectory, &factory, 
                                  &gradient, observer, 
                                  obstol, timeout_seconds,
                                  max_iterations, 
                                  lower_bounds, upper_bounds );

    case THE_OTHER:
#if NLOPT_FOUND
        //TODO add nlopt optimization schemes
        std::cerr << "NLopt optimization schemes are unimplemented"
                  << std::endl;
#else 
        std::cerr << "NLopt optimization libraries are"
                  << " not available, please use a different"
                  << " optimization algorithm" << std::endl;
        return NULL; 
#endif
        return NULL;
    default:
        std::cerr << "Unrecognized optimization scheme";
        return NULL;
    }
    return NULL;
}
void MotionOptimizer::setGoalset( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;
}

void MotionOptimizer::prepareGoalSet(){
    
    //do not subsample if doing goalset run.
    trajectory.startGoalSet();
    
    //add the goal constraint to the constraints vector.
    factory.constraints.push_back( goalset );
}

void MotionOptimizer::finishGoalSet(){
    
    use_goalset = false;
    
    trajectory.endGoalSet();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory.constraints.pop_back();

}


void MotionOptimizer::setLowerBounds( const MatX & lower )
{
    assert( lower.size() == trajectory.M() );
    lower_bounds = lower;
}
void MotionOptimizer::setLowerBounds( const std::vector<double> & lower)
{
    assert( lower.size() == size_t(trajectory.M())  );
    lower_bounds = ConstMatMap(lower.data(), trajectory.M() , 1 );
}
void MotionOptimizer::setLowerBounds( const double * lower)
{
    lower_bounds = ConstMatMap(lower, trajectory.M() , 1 );
}


void MotionOptimizer::setUpperBounds(const MatX & upper )
{
    assert( upper.size() == trajectory.M()  );
    upper_bounds = upper;
}
void MotionOptimizer::setUpperBounds( const std::vector<double> & upper)
{
    assert( upper.size() == size_t(trajectory.M() ) );
    upper_bounds = ConstMatMap(upper.data(), trajectory.M(), 1 );
}
void MotionOptimizer::setUpperBounds( const double * upper)
{
    upper_bounds = ConstMatMap(upper, trajectory.M() , 1 );
}

void MotionOptimizer::setBounds( const MatX & lower, const MatX & upper )
{
    setLowerBounds( lower );
    setUpperBounds( upper );
}
void MotionOptimizer::setBounds( const std::vector<double> & lower, 
                                 const std::vector<double> & upper){
    setLowerBounds( lower );
    setUpperBounds( upper );
}
void MotionOptimizer::setBounds( const double* lower, 
                                    const double* upper){
    setLowerBounds( lower );
    setUpperBounds( upper );
}






}//namespace
