
#include "MotionOptimizer.h"

namespace mopt {


const char* MotionOptimizer::TAG = "MotionOptimizer";

//constructor.
MotionOptimizer::MotionOptimizer(Observer * observer,
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
}

void MotionOptimizer::solve()
{
    
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
        
        //if we are subsampling, optimize at the current
        //  resolution with subsamling,
        //  then optimize without subsampling.
        if (do_subsample ){
            optimize( getOptimizer( algorithm1 ), true );

            OptimizationAlgorithm alg = (algorithm2 == NONE ?
                                         algorithm1 :
                                         algorithm2);
            optimize( getOptimizer( alg ) );
            
        //if we are not subsampling, optimize without subsampling
        }else  {
            optimize( getOptimizer( algorithm1 ));
        }
            
    }

    //If full_global_at_final
    if ( full_global_at_final && do_subsample && N_min < problem.N()){
        optimize( getOptimizer( algorithm1 ));
    }
    
    debug_status( TAG, "solve", "end");
}



void MotionOptimizer::optimize( OptimizerBase * optimizer, 
                                bool subsample )
{
    
    //if the optimizer is NULL, do not evaluate it.
    if ( !optimizer ) { return; }

    debug_status( TAG, "optimize", "start");
    
    //prepare the problem to be run at the current resolution
    problem.prepareRun( subsample );
    
    //do not optimize if the observer throws an error.
    //TODO either throw error or say what exactly happened
    if ( !optimizer->notify( INIT ) ){
        optimizer->solve();
    } else {
        debug << "Observer threw error on INIT, stopping optimization\n";
    }
    
    //the run is over, so tell the problem to clean up stuff
    //  pertaning to the previous run.
    problem.endRun();
    
    //notify the the finish
    optimizer->notify( FINISH );

    //delete the used optimizer.
    delete optimizer;


    debug_status( TAG, "optimize", "end");

}

OptimizerBase * MotionOptimizer::getOptimizer(OptimizationAlgorithm alg )
{
    //create the optimizer

    //If the algorith is CovariantChomp, or if the Algorithm is a
    //  chomp variant with a covaraint trajectory, do CovariantChomp
    if ( alg == LOCAL_CHOMP &&
         problem.isCovariantOptimization() &&
         !problem.isSubsampled() )
    {
        alg = CHOMP;
    }
    
    //Collision constraints should not be run with 'chomp-like'
    //  algorithms because 
    //  1) they have not been hooked up to it
    //  2) they were designed with collision as a soft constraint.
    if ( (alg == LOCAL_CHOMP || alg == CHOMP || alg == TEST ) &&
          problem.collision_constraint)
    {
        //TODO throw error
        problem.collision_constraint = false;
    }

    if ( alg == CHOMP ){
        ChompOptimizer * opt = new ChompOptimizer(
                                  problem,
                                  observer, 
                                  obstol, timeout_seconds,
                                  max_iterations);
        if (alpha >= 0){ opt->setAlpha( alpha ); }
        return opt;
    } else if ( alg == TEST ){
        TestOptimizer * opt = new TestOptimizer(
                                      problem, observer, 
                                      obstol, timeout_seconds,
                                      max_iterations);
        
        if (alpha >= 0){ opt->setAlpha( alpha ); }
        return opt;

    } else if ( alg == LOCAL_CHOMP ){
        ChompLocalOptimizer * opt = new ChompLocalOptimizer(
                                      problem, observer, 
                                      obstol, timeout_seconds,
                                      max_iterations);
        
        if (alpha >= 0){ opt->setAlpha( alpha ); }
        return opt;
        
    } else if ( alg > TEST && alg < NONE){
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

#ifdef NLOPT_FOUND
nlopt::algorithm getNLoptAlgorithm( OptimizationAlgorithm alg )
{

    switch (alg){
        case MMA_NLOPT:   return nlopt::LD_MMA;
        case CCSAQ_NLOPT: return nlopt::LD_CCSAQ;
        case SLSQP_NLOPT: return nlopt::LD_SLSQP;
        case LBFGS_NLOPT: return nlopt::LD_LBFGS;
        case VAR1_NLOPT:  return nlopt::LD_VAR1;
        case VAR2_NLOPT:  return nlopt::LD_VAR2;
        case TNEWTON_PRECOND_RESTART_NLOPT:
            return nlopt::LD_TNEWTON_PRECOND_RESTART;
        case TNEWTON_RESTART_NLOPT: return nlopt::LD_TNEWTON_RESTART;
        case TNEWTON_PRECOND_NLOPT: return nlopt::LD_TNEWTON_PRECOND;
        case TNEWTON_NLOPT: return nlopt::LD_TNEWTON;
        //TODO Throw error.
        default: return nlopt::LD_MMA;
    }
}
#endif

std::string algorithmToString( OptimizationAlgorithm alg )
{

    switch (alg){
    case LOCAL_CHOMP:   return "LOCAL_CHOMP";
    case CHOMP:         return "CHOMP";
    case TEST:          return "TEST";
#ifdef NLOPT_FOUND
    case MMA_NLOPT:     return "MMA";
    case CCSAQ_NLOPT:   return "CCSAQ";
    case SLSQP_NLOPT:   return "SLSQP";
    case LBFGS_NLOPT:   return "LBFGS";
    case TNEWTON_NLOPT: return "NEWTON";
    case VAR1_NLOPT:    return "VAR1";
    case VAR2_NLOPT:    return "VAR2";
    case TNEWTON_PRECOND_RESTART_NLOPT: return "TNEWTON_PRECOND_RESTART";
    case TNEWTON_RESTART_NLOPT:         return "TNEWTON_RESTART";
    case TNEWTON_PRECOND_NLOPT:         return "TNEWTON_PRECOND";
#endif
    //TODO Throw error.
    case NONE:          return "NONE";
    default: 
        return "ERROR";
    }
}

OptimizationAlgorithm algorithmFromString( const std::string & str )
{
    if ( str == "LOCAL_CHOMP" ){return LOCAL_CHOMP; }
    else if ( str == "CHOMP" ){ return CHOMP;}
    else if ( str == "TEST" ){  return TEST;}

#ifdef NLOPT_FOUND
    else if ( str == "MMA" ){   return MMA_NLOPT;}
    else if ( str == "CCSAQ" ){ return CCSAQ_NLOPT;}
    else if ( str == "SLSQP" ){ return SLSQP_NLOPT;}
    else if ( str == "LBFGS" ){ return LBFGS_NLOPT;}
    else if ( str == "NEWTON" ){return TNEWTON_NLOPT;}
    else if ( str == "VAR1" ){  return VAR1_NLOPT;}
    else if ( str == "VAR2" ){  return VAR2_NLOPT;}
    else if ( str == "TNEWTON_RESTART" ){ return TNEWTON_RESTART_NLOPT;}
    else if ( str == "TNEWTON_PRECOND" ){ return TNEWTON_PRECOND_NLOPT;}
    else if ( str == "TNEWTON_PRECOND_RESTART" ){
       return TNEWTON_PRECOND_RESTART_NLOPT;
    }
#endif
    
    else { return NONE; }
}

//simple getters and setters
void MotionOptimizer::setNMax( int n_max ){ N_max = n_max; }
int  MotionOptimizer::getNMax( ) const { return N_max; }

void MotionOptimizer::setGoalset( Constraint * goal)
{
    problem.setGoalset(goal);
}
const Constraint * MotionOptimizer::getGoalset() const 
{
    return problem.getGoalset();
}

void MotionOptimizer::setTimeoutSeconds( double s )
{ 
    timeout_seconds = s;
}
double MotionOptimizer::getTimeoutSeconds() const
{ 
    return timeout_seconds;
}

void MotionOptimizer::setMaxIterations( size_t max )
{ 
    max_iterations = max; 
}
size_t MotionOptimizer::getMaxIterations() const
{ 
    return max_iterations;
}

void MotionOptimizer::setFunctionTolerance( double tol )
{ 
    obstol = tol;
}
double MotionOptimizer::getFunctionTolerance() const
{ 
    return obstol;
}

void MotionOptimizer::setAlgorithm(OptimizationAlgorithm a1,
                                   OptimizationAlgorithm a2)
{ 
    algorithm1 = a1;
    algorithm2 = a2;
}



void MotionOptimizer::setAlgorithm1(OptimizationAlgorithm a1)
{ 
    algorithm1 = a1;
}
void MotionOptimizer::setAlgorithm2(OptimizationAlgorithm a2)
{ 
    algorithm2 = a2;
}

void MotionOptimizer::setAlgorithm( const std::string & a1, 
                                    const std::string & a2)
{
    algorithm1 = algorithmFromString( a1 ); 
    algorithm2 = algorithmFromString( a2 ); 
}
void MotionOptimizer::setAlgorithm1( const std::string & a1)
{
    algorithm1 = algorithmFromString( a1 ); 
}

void MotionOptimizer::setAlgorithm2( const std::string & a2) 
{
    algorithm2 = algorithmFromString( a2 ); 
}

OptimizationAlgorithm MotionOptimizer::getAlgorithm() const
{ 
    return algorithm1;
}
OptimizationAlgorithm MotionOptimizer::getAlgorithm1() const
{ 
    return algorithm1;
}
OptimizationAlgorithm MotionOptimizer::getAlgorithm2() const
{ 
    return algorithm2;
}

void MotionOptimizer::dontSubsample()
{ 
    do_subsample = false;
}
void MotionOptimizer::doSubsample()
{ 
    do_subsample = true;
}
void MotionOptimizer::setSubsample( bool subsample )
{ 
    do_subsample = subsample;
}

void MotionOptimizer::setAlpha( double a )
{ 
    alpha = a;
}
double MotionOptimizer::getAlpha() const
{ 
    return alpha;
}

void MotionOptimizer::doFullGlobalAtFinal()
{ 
    full_global_at_final = true;
}
void MotionOptimizer::dontFullGlobalAtFinal()
{ 
    full_global_at_final = false;
}
bool MotionOptimizer::getFullGlobalAtFinal() const
{ 
    return full_global_at_final;
}

Trajectory & MotionOptimizer::getTrajectory()
{ 
    return problem.trajectory;
}

const Trajectory & MotionOptimizer::getTrajectory() const 
{ 
    return problem.trajectory;
}
void MotionOptimizer::setTrajectory( const Trajectory & trajectory )
{ 
    problem.trajectory = trajectory;
}

void MotionOptimizer::setCollisionFunction(CollisionFunction * coll_func)
{
    assert( coll_func != NULL );
    
    problem.collision_function = coll_func;
}

const CollisionFunction * MotionOptimizer::getCollisionFunction()
const
{
    return problem.collision_function;
}

void MotionOptimizer::setObserver( Observer * obs )
{ 
    observer = obs;
}

Observer * MotionOptimizer::getObserver()
{ 
    return observer;
}
const Observer * MotionOptimizer::getObserver() const 
{ 
    return observer;
}

void MotionOptimizer::doCovariantOptimization()
{ 
    problem.is_covariant = true;
}
void MotionOptimizer::dontCovariantOptimization()
{ 
    problem.is_covariant = false;
}
void MotionOptimizer::setCovariantOptimization( bool covariant)
{ 
    problem.is_covariant = covariant;
}
bool MotionOptimizer::isCovariantOptimization() const
{ 
    return problem.isCovariantOptimization();
}


void MotionOptimizer::doCollisionConstraint()
{
    problem.collision_constraint = true;
}
void MotionOptimizer::dontCollisionConstraint()
{
    problem.collision_constraint = false;
}
void MotionOptimizer::setCollisionConstraint( 
                                bool do_collision_constraint )
{
    problem.collision_constraint = do_collision_constraint;
}
bool MotionOptimizer::isCollisionConstraint() const
{
    return problem.collision_constraint;
}


}//namespace
