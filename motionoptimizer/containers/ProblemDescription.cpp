

#include "ProblemDescription.h"

namespace chomp {
ProblemDescription::ProblemDescription() :
    goalset( NULL ),
    ok_to_run( false ),
    use_goalset( false )
{
}


void ProblemDescription::prepareRun(){
    
    ok_to_run = true;
    
    if( use_goalset ){ stopGoalset(); }

    if ( !factory.empty() ){ factory.getAll( trajectory.N() ); }

    //are we going to use the goalset on this iteration?
    //  If we are subsampling, then even if there is a goalset
    //  to use, do not use it.
    bool use_goalset = goalset && !trajectory.isSubsampled();
    
    if ( use_goalset ){ startGoalset(); }

    gradient.prepareRun( trajectory, use_goalset );
}


void ProblemDescription::upsample(){
    
    ok_to_run = false;
    if( use_goalset ){ stopGoalset(); }

    trajectory.upsample();
}

void ProblemDescription::subsample(){

    ok_to_run = false;
    if( use_goalset ){ stopGoalset(); }

    trajectory.subsample();
}

void ProblemDescription::stopSubsample(){

    ok_to_run = false;
    trajectory.endSubsample();
}

double ProblemDescription::evaluateGradient( MatX & g )
{
    if ( !ok_to_run ){ prepareRun(); }
    gradient.evaluate( trajectory, g );

    return gradient.evaluateObjective(trajectory);
}

double ProblemDescription::evaluateGradient( const double * xi,
                                                   double * g )
{
    if ( !ok_to_run ){ prepareRun(); }
    
    trajectory.setData( xi );

    if ( !g ){ return gradient.evaluateObjective( trajectory ); }
            
    MatMap g_map( g, trajectory.N(), trajectory.M() );
    gradient.evaluate( trajectory, g_map );
    
    return gradient.evaluateObjective( trajectory );
}



double ProblemDescription::evaluateConstraint( MatX & h )
{
    if ( factory.empty() ) {return 0; }

    if ( !ok_to_run ){ prepareRun(); }
    
    h.resize( factory.numOutput(), 1 );

    return factory.evaluate( trajectory, h);
}

double ProblemDescription::evaluateConstraint( MatX & h, MatX & H )
{
    if ( factory.empty() ) {return 0; }

    if ( !ok_to_run ){ prepareRun(); }

    h.resize( factory.numOutput(), 1 );
    H.resize( size(), factory.numOutput() );
    
    return factory.evaluate( trajectory, h, H );
}

double ProblemDescription::evaluateConstraint( const double * xi, 
                                                     double * h,
                                                     double * H )
{
    if ( !ok_to_run ){ prepareRun(); }
    if ( factory.empty() ) {return 0; }
    assert( h ); // make sure that h is not NULL
    
    trajectory.setData( xi );
    MatMap h_map( h, factory.numOutput(), 1 );

    if ( !H ){ return factory.evaluate( trajectory, h_map ); }
        
    MatMap H_map( H, trajectory.size(), factory.numOutput() );

    return factory.evaluate(trajectory, h_map, H_map );
}

bool ProblemDescription::evaluateConstraint( MatX & h_t,
                                             MatX & H_t,
                                             int t )
{
    
    Constraint * c = factory.getConstraint( t );
    if ( c == NULL || c->numOutputs() == 0 ) { return false; }

    c->evaluateConstraints( trajectory.row( t ), h_t, H_t );
    return true;
}

double ProblemDescription::evaluateObjective()
{
    if ( !ok_to_run ){ prepareRun(); }

    //TODO include the fextra term
    return gradient.evaluateObjective( trajectory );
}

double ProblemDescription::evaluateObjective( const double * xi )
{
    if ( !ok_to_run ){ prepareRun(); }

    trajectory.setData( xi );
    
    //TODO include the fextra term.
    return gradient.evaluateObjective( trajectory );
}

int ProblemDescription::getConstraintDims()
{
    if( !ok_to_run ){ prepareRun(); }
    return factory.numOutput();
}

void ProblemDescription::startGoalset(){ 

    //do not subsample if doing goalset run.
    trajectory.startGoalset();

    //add the goal constraint to the constraints vector.
    factory.addGoalset( goalset );
    
    use_goalset = true;
}
    
void ProblemDescription::stopGoalset(){ 

    use_goalset = false;
    
    trajectory.endGoalset();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory.removeGoalset();
    
}
    
void ProblemDescription::copyToTrajectory( const double * data )
{
    trajectory.restoreData();
    trajectory.copyToData( data );
}
void ProblemDescription::copyToTrajectory( const std::vector<double> data )
{
    trajectory.restoreData();
    trajectory.copyToData( data );
}



}//namespace
