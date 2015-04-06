

#include "ProblemDescription.h"

namespace chomp {
ProblemDescription::ProblemDescription() :
    goalset( NULL ),
    ok_to_sample( true ),
    use_goalset( false ),
    is_covariant( false ),
    doing_covariant( false )
{
    TIMER_START( "total" );
}
ProblemDescription::~ProblemDescription(){}


void ProblemDescription::prepareSample()
{
    ok_to_sample = true;
    if( use_goalset ){ stopGoalset(); }
}

void ProblemDescription::prepareRun()
{
    
    debug_status( "ProblemDescription", "prepareRun", "start");
    
    if( use_goalset ){ stopGoalset(); }
    
    //if the factory is not empty of constraints,
    //  get all of the constraints for the current resolution
    if ( !factory.empty() ){ factory.getAll( trajectory.N() ); }

    //are we going to use the goalset on this iteration?
    //  If we are subsampling, then even if there is a goalset
    //  to use, do not use it.
    bool use_goalset = goalset && !trajectory.isSubsampled();
    if ( use_goalset ){ startGoalset(); }
    
    doing_covariant = ( is_covariant && !trajectory.isSubsampled() );

    //prepare the gradient for a run at this resolution
    gradient.prepareRun( trajectory, use_goalset, doing_covariant );
    
    //are we doing covariant optimization at this stage of optimiation?
    //  do not do covariant optimization if there is subsampling
    if ( doing_covariant ){
        trajectory.getCovariantTrajectory(gradient.getMetric(),
                                          covariant_trajectory);
    }

    debug_status( "ProblemDescription", "prepareRun", "end");
}


void ProblemDescription::upsample()
{
    if( !ok_to_sample ){ prepareSample(); }
    trajectory.upsample();
}

void ProblemDescription::subsample()
{
    if( !ok_to_sample ){ prepareSample(); }
    trajectory.subsample();
}

void ProblemDescription::stopSubsample()
{
    trajectory.endSubsample();
}

double ProblemDescription::evaluateGradient( MatX & g )
{

    TIMER_START( "gradient" );
    
    if ( doing_covariant ){ 
        prepareCovariant();
        gradient.evaluate( trajectory, g, &covariant_trajectory );
        return gradient.evaluateObjective( covariant_trajectory, true );
    }
    
    gradient.evaluate( trajectory, g );

    const double val = gradient.evaluateObjective(trajectory);

    TIMER_STOP( "gradient" );
    
    return val;
}

double ProblemDescription::evaluateGradient( const double * xi,
                                                   double * g )
{
    TIMER_START( "gradient" );
    if ( doing_covariant ){ prepareCovariant( xi ); }
    else { trajectory.setData( xi ); }
    
    if ( !g ){ 
        const double val = gradient.evaluateObjective( trajectory );
        TIMER_STOP( "gradient" );
        return val;
    }
            
    MatMap g_map( g, trajectory.N(), trajectory.M() );

    if ( doing_covariant ){
        gradient.evaluate( trajectory, g_map, &covariant_trajectory );
        const double val = 
            gradient.evaluateObjective( covariant_trajectory, true );
        TIMER_STOP( "gradient" );
        return val;
    }

    gradient.evaluate( trajectory, g_map );
    const double val = gradient.evaluateObjective( trajectory );

    TIMER_STOP( "gradient" );

    return val;
}



double ProblemDescription::evaluateConstraint( MatX & h )
{
    if ( factory.empty() ) {return 0; }

    TIMER_START( "constraint" );
    if ( doing_covariant ){ prepareCovariant(); }
    
    h.resize( factory.numOutput(), 1 );

    const double value = factory.evaluate( trajectory, h);

    TIMER_STOP( "constraint" );

    return value;
}

double ProblemDescription::evaluateConstraint( MatX & h, MatX & H )
{
    if ( factory.empty() ) {return 0; }

    TIMER_START( "constraint" );

    if ( doing_covariant ){ prepareCovariant(); }

    h.resize( factory.numOutput(), 1 );
    H.resize( size(), factory.numOutput() );
    
    double magnitude = factory.evaluate( trajectory, h, H );

    if ( doing_covariant ){ 
        //TODO find out if this is correct
        MatMap H_map( H.data(), trajectory.N(),
                      trajectory.M()*factory.numOutput() );
        gradient.getMetric().multiplyLowerInverse( H );
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
    
    if ( doing_covariant ){ prepareCovariant( xi );}
    else { trajectory.setData( xi ); }

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
        gradient.getMetric().multiplyLowerInverse( H_map2 );
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

double ProblemDescription::evaluateObjective()
{
    if ( doing_covariant ){ prepareCovariant(); }

    //TODO include the fextra term
    return gradient.evaluateObjective( trajectory );
}

double ProblemDescription::evaluateObjective( const double * xi )
{

    if ( doing_covariant ){ 
        covariant_trajectory.setData(xi);
        return gradient.evaluateObjective( covariant_trajectory, true );
    }

    trajectory.setData( xi );
    //TODO include the fextra term.
    return gradient.evaluateObjective( trajectory );
}

int ProblemDescription::getConstraintDims()
{
    return factory.numOutput();
}

void ProblemDescription::startGoalset()
{ 
    //startGoalSet can effect the data in trajectory,
    //  so we must check for that.
    if( !ok_to_sample ){ prepareSample(); }

    //do not subsample if doing goalset run.
    trajectory.startGoalset();

    //add the goal constraint to the constraints vector.
    factory.addGoalset( goalset );
    
    ok_to_sample = false;
}
    
void ProblemDescription::stopGoalset()
{ 

    use_goalset = false;
    
    //Restore the trajectory to the non-goalset type
    trajectory.endGoalset();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory.removeGoalset();
    
}
    
void ProblemDescription::copyToTrajectory( const double * data )
{
    if( doing_covariant ){
        covariant_trajectory.copyToData( data );
        covariant_trajectory.getNonCovariantTrajectory( 
                                    gradient.getMetric(), 
                                    trajectory );
    }else {
        trajectory.copyToData( data );
    }

    if( !ok_to_sample ){ prepareSample(); }
    
}

void ProblemDescription::copyToTrajectory( const std::vector<double> data )
{

    if( doing_covariant ){
        covariant_trajectory.copyToData( data );
        covariant_trajectory.getNonCovariantTrajectory( 
                                    gradient.getMetric(), 
                                    trajectory );
    }else {
        trajectory.copyToData( data );
    }

    if( !ok_to_sample ){ prepareSample(); }
    
}

void ProblemDescription::prepareCovariant(const double * xi )
{
    if (xi){ 
        covariant_trajectory.setData( xi );
    }
    
    covariant_trajectory.getNonCovariantTrajectory( gradient.getMetric(), 
                                                    trajectory );
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
