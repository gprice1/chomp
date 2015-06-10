

#include "ProblemDescription.h"


namespace mopt {

    
const char* ProblemDescription::TAG = "ProblemDescription";


ProblemDescription::ProblemDescription() :
    goalset( NULL ),
    use_goalset( false ),
    is_covariant( false ),
    doing_covariant( false )
{
    TIMER_START( "total" );
}
ProblemDescription::~ProblemDescription(){}


void ProblemDescription::copyTrajectoryTo( double * data )
{
    if (doing_covariant){ covariant_trajectory.copyDataTo( data );}
    else { trajectory.copyDataTo( data ); }
}

void ProblemDescription::copyTrajectoryTo( std::vector<double> & data )
{
    if (doing_covariant){ covariant_trajectory.copyDataTo( data );}
    else { trajectory.copyDataTo( data ); }
}


void ProblemDescription::prepareRun( bool subsample)
{
    
    debug_status( "ProblemDescription", "prepareRun", "start");
    
    if ( subsample ){ trajectory.subsample(); }
    
    //if the factory is not empty of constraints,
    //  get all of the constraints for the current resolution
    if ( !factory.empty() ){ factory.getAll( trajectory.N() ); }

    //are we going to use the goalset on this iteration?
    //  If we are subsampling, then even if there is a goalset
    //  to use, do not use it.
    bool use_goalset = goalset && !subsample;
    if ( use_goalset ){ 
        //do not subsample if doing goalset run.
        trajectory.startGoalset();

        //add the goal constraint to the constraints vector.
        factory.addGoalset( goalset );
    }
    
    //are we doing covariant on the current 
    //  optimization? Do not do covariant optimization
    //  if we are subsampling 
    doing_covariant = is_covariant && !subsample;

    metric.initialize( trajectory.fullN(),
                       trajectory.getObjectiveType(),
                       false,
                       use_goalset );
    if ( subsample ){
        subsampled_metric.initialize( trajectory.N(),
                                      trajectory.getObjectiveType(),
                                      true);
    }


    //prepare the gradient for a run at this resolution
    smoothness_function.prepareRun( trajectory, metric);
    
    //are we doing covariant optimization at this stage of optimiation?
    //  do not do covariant optimization if there is subsampling
    if ( doing_covariant ){
        trajectory.getCovariantTrajectory(metric, covariant_trajectory);
    }

    debug_status( "ProblemDescription", "prepareRun", "end");
}

void ProblemDescription::endRun()
{
    debug_status( "ProblemDescription", "endRun", "start");
    
    if ( doing_covariant ){
        covariant_trajectory.getNonCovariantTrajectory( metric,
                                                        trajectory );
    }
    
    if( use_goalset ){     
        //Restore the trajectory to the non-goalset type
        trajectory.endGoalset();

        //remove the goal constraint, so that it is not deleted along
        //  with the other constraints.
        factory.removeGoalset();
    }
    
    if ( trajectory.isSubsampled() ){ 
        trajectory.endSubsample();
    }
    
    debug_status( "ProblemDescription", "endRun", "end");
}


void ProblemDescription::upsample()
{
    trajectory.upsample();
}

double ProblemDescription::evaluateCollisionFunction( const double * xi,
                                                            double * g)
{
    //TODO throw error.
    assert( collision_function );

    prepareData( xi );
    
    double value;
    if ( g ){
        MatMap g_map( g, N(), M() );
        value = collision_function->evaluate( trajectory, g_map );
        if ( doing_covariant ){ metric.multiplyLowerInverse( g_map ); } 
    }else {
        value = collision_function->evaluate( trajectory );
    }
    
    return value;
}


double ProblemDescription::evaluateObjective ( const double * xi,
                                               double * g )
{
    TIMER_START( "gradient" );
    
    if ( xi ){ prepareData( xi ); }
    else     { prepareData();     }

    double value;
    
    if ( g ) {
        value = computeObjective( MatMap(g,
                                         trajectory.N(),
                                         trajectory.M() ) );
    } else {
        value = computeObjective( MatX(0,0) );
    }
    
    TIMER_STOP( "gradient" );

    return value;
}


double ProblemDescription::evaluateConstraint( MatX & h )
{
    if ( factory.empty() ) {return 0; }

    TIMER_START( "constraint" );
    
    prepareData();
    
    h.resize( factory.numOutput(), 1 );

    const double value = factory.evaluate( trajectory, h);

    TIMER_STOP( "constraint" );

    return value;
}

double ProblemDescription::evaluateConstraint( MatX & h, MatX & H )
{
    if ( factory.empty() ) {return 0; }

    TIMER_START( "constraint" );

    prepareData();
    
    h.resize( factory.numOutput(), 1 );
    H.resize( size(), factory.numOutput() );
    
    double magnitude = factory.evaluate( trajectory, h, H );

    if ( doing_covariant ) {
        metric.multiplyLowerInverse( 
                    MatMap ( H.data(), trajectory.N(),
                             trajectory.M()*factory.numOutput() ) );
    }
    
    TIMER_STOP( "constraint" );

    return magnitude;
    
}

double ProblemDescription::evaluateConstraint( const double * xi, 
                                                     double * h,
                                                     double * H )
{
    if ( factory.empty() ) {return 0; }

    TIMER_START( "constraint" );
    
    assert( h ); // make sure that h is not NULL
    
    prepareData( xi );
    
    MatMap h_map( h, factory.numOutput(), 1 );

    if ( !H ){ 
        const double val = factory.evaluate( trajectory, h_map );
        TIMER_STOP( "constraint" );
        return val;
    }
        
    MatMap H_map( H, trajectory.size(), factory.numOutput() );

    double magnitude = factory.evaluate(trajectory, h_map, H_map );

    if( doing_covariant ){
        //TODO find out if this is correct
        MatMap H_map2( H, trajectory.N(),
                       trajectory.M()*factory.numOutput() );
        metric.multiplyLowerInverse( H_map2 );
    }
    
    TIMER_STOP( "constraint" );
    return magnitude;
}

bool ProblemDescription::evaluateConstraint( MatX & h_t,
                                             MatX & H_t,
                                             int t )
{
    TIMER_START( "constraint" );
    
    //TODO make this throw an error
    //  A covariant constraint jacobian H should not be local to
    //  the timestep, so we cannot use this method
    //  of constraint detection for covariant optimization
    assert( !is_covariant );
    
    Constraint * c = factory.getConstraint( t );
    if ( c == NULL || c->numOutputs() == 0 ) {
        TIMER_STOP( "constraint" );
        return false;
    }

    c->evaluateConstraints( trajectory.row( t ), h_t, H_t );

    TIMER_STOP( "constraint" );
    
    return true;
}



int ProblemDescription::getConstraintDims()
{
    return factory.numOutput();
}
    
void ProblemDescription::copyToTrajectory( const double * data )
{
    if( doing_covariant ){
        covariant_trajectory.copyToData( data );
    }else {
        trajectory.copyToData( data );
    }
}

void ProblemDescription::copyToTrajectory( const std::vector<double> data )
{
    if( doing_covariant ){
        covariant_trajectory.copyToData( data );
    }else {
        trajectory.copyToData( data );
    }
}

void ProblemDescription::prepareData(const double * xi )
{
    if ( doing_covariant ){ 
        if (xi){ covariant_trajectory.setData( xi ); }
        
        covariant_trajectory.getNonCovariantTrajectory(
                                     metric, 
                                     trajectory );
    }else if ( xi ) {
        trajectory.setData( xi );
    }
}
   
void ProblemDescription::getFullBounds( std::vector< double > & lower,
                                        std::vector< double > & upper )
{
    //resize the bounds
    lower.resize( size() );
    upper.resize( size() );
    
    if ( doing_covariant ){
        //get the covariant bounds
        metric.solveCovariantBounds( this->lower_bounds, 
                                     this->upper_bounds,
                                     MatMap( lower.data(), N(), M() ),
                                     MatMap( upper.data(), N(), M() ));
        
    } else {
        if ( lower_bounds.size() == M() ){
            MatMap lower_map( lower.data(), N(), M() );
            for ( int i = 0; i < N(); i++ ){
                lower_map.row( i ) = lower_bounds;
            }
        }
        if ( upper_bounds.size() == M() ){
            MatMap upper_map( upper.data(), N(), M() );
            for ( int i = 0; i < N(); ++i ){
                upper_map.row( i ) = upper_bounds;
            }
        }
    }
}

void ProblemDescription::getTimes( 
        std::vector< std::pair<std::string, double> > & times ) const
{
#ifdef DO_TIMING
    timer.getAllTotal( times );
#else 
    std::cout << "Timing information is not available"
              << " because it was disabled at compile time.\n";
#endif
}

void ProblemDescription::printTimes( bool verbose ) const
{
#ifdef DO_TIMING
    std::vector< std::pair<std::string, double> > times;
    timer.getAllTotal( times );

    if (times.size() > 0 ){
        for ( size_t i = 0; i < times.size() ; i ++ )
        {
            if (verbose){ std::cout << times[i].first << ": "; }
            std::cout << times[i].second;
            if ( i != times.size() - 1 ){ std::cout << ", "; }
        }
    }else {
        std::cout << "No timing information available.";
    }

#else 
    std::cout << "Timing information is not available"
              << " because it was disabled at compile time.\n";
#endif

}
std::string ProblemDescription::getTimesString( bool verbose )
{
#ifdef DO_TIMING
    std::ostringstream stream;
    std::vector< std::pair<std::string, double> > times;
    timer.stop("total");
    timer.start("total");
    timer.getAllTotal( times );

    if (times.size() > 0 ){
        for ( size_t i = 0; i < times.size() ; i ++ )
        {
            if (verbose){ stream << times[i].first << ":"; }
            stream << times[i].second;
            if ( i != times.size() - 1 ){ stream << ", "; }
        }
    }else {
        stream << "No timing information available.";
    }

    return stream.str();

#else 
    std::cout << "Timing information is not available"
              << " because it was disabled at compile time.\n";
    return std::string();
#endif

}

//undefine the timer macros
#ifdef TIMER_START
    #undef TIMER_START
#endif

#ifdef TIMER_STOP        
    #undef TIMER_STOP
#endif

#ifdef TIMER_GET_TOTAL
    #undef TIMER_GET_TOTAL
#endif


}//namespace
