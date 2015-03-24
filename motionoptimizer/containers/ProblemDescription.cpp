

#include "ProblemDescription.h"

namespace chomp {
ProblemDescription::ProblemDescription() :
    goalset( NULL ),
    ok_to_run( false ),
    ok_to_sample( true ),
    use_goalset( false ),
    is_covariant( false )
{
}
ProblemDescription::~ProblemDescription(){}


void ProblemDescription::prepareSample()
{
    ok_to_sample = true;
    if( use_goalset ){ stopGoalset(); }
    
    if ( is_covariant ){
        covariant_trajectory.restoreData();
    } else {
        trajectory.restoreData();
    }
}

void ProblemDescription::prepareRun()
{
    
    ok_to_run = true;
    
    if( use_goalset ){ stopGoalset(); }

    if ( !factory.empty() ){ factory.getAll( trajectory.N() ); }

    //are we going to use the goalset on this iteration?
    //  If we are subsampling, then even if there is a goalset
    //  to use, do not use it.
    bool use_goalset = goalset && !trajectory.isSubsampled();
    
    if ( use_goalset ){ startGoalset(); }

    gradient.prepareRun( trajectory, use_goalset );

    if ( is_covariant ){
        trajectory.getCovariantTrajectory(gradient.getLMatrix(),
                                          covariant_trajectory);
    }
}


void ProblemDescription::upsample()
{
    
    if( !ok_to_sample ){ prepareSample(); }
    
    ok_to_run = false;

    trajectory.upsample();
}

void ProblemDescription::subsample()
{

    if( !ok_to_sample ){ prepareSample(); }
    
    ok_to_run = false;

    trajectory.subsample();
}

void ProblemDescription::stopSubsample()
{

    ok_to_run = false;
    
    trajectory.endSubsample();
}

double ProblemDescription::evaluateGradient( MatX & g )
{
    if ( !ok_to_run ){ prepareRun(); }
    if ( is_covariant ){ 
        prepareCovariant();
        gradient.evaluate( trajectory, g, &covariant_trajectory );
        return gradient.evaluateObjective( covariant_trajectory, true );
    }
    
    gradient.evaluate( trajectory, g );
    return gradient.evaluateObjective(trajectory);
}

double ProblemDescription::evaluateGradient( const double * xi,
                                                   double * g )
{
    if ( !ok_to_run ){ prepareRun(); }
    if ( is_covariant ){
        prepareCovariant( xi );
    } else {
        trajectory.setData( xi );
    }
    
    //setData() is a dangerous call that will throw segfaults
    //  if sampling occurs on top of it.
    ok_to_sample = false;

    if ( !g ){ return gradient.evaluateObjective( trajectory ); }
            
    MatMap g_map( g, trajectory.N(), trajectory.M() );

    if ( is_covariant ){
        gradient.evaluate( trajectory, g_map, &covariant_trajectory );
        return gradient.evaluateObjective( covariant_trajectory, true );
    }

    gradient.evaluate( trajectory, g_map );
    return gradient.evaluateObjective( trajectory );
}



double ProblemDescription::evaluateConstraint( MatX & h )
{
    if ( factory.empty() ) {return 0; }
    if ( !ok_to_run ){ prepareRun(); }
    if ( is_covariant ){ prepareCovariant(); }
    
    h.resize( factory.numOutput(), 1 );

    return factory.evaluate( trajectory, h);
}

double ProblemDescription::evaluateConstraint( MatX & h, MatX & H )
{
    if ( factory.empty() ) {return 0; }

    if ( !ok_to_run ){ prepareRun(); }
    
    if ( is_covariant ){ prepareCovariant(); }

    h.resize( factory.numOutput(), 1 );
    H.resize( size(), factory.numOutput() );
    
    double magnitude = factory.evaluate( trajectory, h, H );

    skylineCholMultiplyInverse( getLMatrix(), H );

    return magnitude;
    
}

double ProblemDescription::evaluateConstraint( const double * xi, 
                                                     double * h,
                                                     double * H )
{
    if ( !ok_to_run ){ prepareRun(); }
    if ( factory.empty() ) {return 0; }
    assert( h ); // make sure that h is not NULL
    
    if ( is_covariant ){
        prepareCovariant( xi );
    }else {
        trajectory.setData( xi );
    }

    //setData() is a dangerous call that will throw segfaults
    //  if sampling occurs on top of it.
    ok_to_sample = false;

    MatMap h_map( h, factory.numOutput(), 1 );

    if ( !H ){ return factory.evaluate( trajectory, h_map ); }
        
    MatMap H_map( H, trajectory.size(), factory.numOutput() );

    return factory.evaluate(trajectory, h_map, H_map );
}

bool ProblemDescription::evaluateConstraint( MatX & h_t,
                                             MatX & H_t,
                                             int t )
{
    //TODO make this throw an error
    //  A covariant constraint jacobian H should not be local to
    //  the timestep, so we cannot use this method
    //  of constraint detection for covariant optimization
    assert( !is_covariant );
    
    Constraint * c = factory.getConstraint( t );
    if ( c == NULL || c->numOutputs() == 0 ) { return false; }

    c->evaluateConstraints( trajectory.row( t ), h_t, H_t );
    return true;
}

double ProblemDescription::evaluateObjective()
{
    if ( !ok_to_run ){ prepareRun(); }
    if ( is_covariant ){ prepareCovariant(); }

    //TODO include the fextra term
    return gradient.evaluateObjective( trajectory );
}

double ProblemDescription::evaluateObjective( const double * xi )
{
    if ( !ok_to_run ){ prepareRun(); }
    
    if ( is_covariant ){ 
        prepareCovariant();
    } else {
        trajectory.setData( xi );
    }
    
    //setData() is a dangerous call that will throw segfaults
    //  if sampling occurs on top of it.
    ok_to_sample = false;

    //TODO include the fextra term.
    return gradient.evaluateObjective( trajectory );
}

int ProblemDescription::getConstraintDims()
{
    if( !ok_to_run ){ prepareRun(); }
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
    
    use_goalset = true;

    ok_to_sample = false;
}
    
void ProblemDescription::stopGoalset()
{ 

    use_goalset = false;
    if( !ok_to_sample ){ prepareSample(); }
    
    trajectory.endGoalset();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory.removeGoalset();
    
}
    
void ProblemDescription::copyToTrajectory( const double * data )
{
    if( !ok_to_sample ){ prepareSample(); }
    trajectory.copyToData( data );
}
void ProblemDescription::copyToTrajectory( const std::vector<double> data )
{
    if( !ok_to_sample ){ prepareSample(); }
    trajectory.copyToData( data );
}

void ProblemDescription::prepareCovariant(const double * xi )
{
    if (xi){ 
        covariant_trajectory.setData( xi );
    }
    
    covariant_trajectory.getNonCovariantTrajectory( gradient.getLMatrix(), 
                                                    trajectory );
    
}
    


}//namespace
