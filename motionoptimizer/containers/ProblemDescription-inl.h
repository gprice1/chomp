
template <class Derived>
double ProblemDescription::evaluateObjective( 
                           const Eigen::MatrixBase<Derived> & g )
{
    TIMER_START( "gradient" );
    
    prepareData();
    
    double value = computeObjective( g );
    
    TIMER_STOP( "gradient" );
    
    return value;

}

template <class Derived>
inline double ProblemDescription::computeObjective( 
                           const Eigen::MatrixBase<Derived> & g )
{
    
    double value;
    if (g.size() == 0 ){
            value = smoothness_function.evaluate(trajectory, metric);
    
        if ( collision_function && !collision_constraint ){ 
            value += collision_function->evaluate( trajectory );
        }
    }else {
        value = smoothness_function.evaluate( trajectory, metric, g );
        
        if ( collision_function && !collision_constraint ){
            value += collision_function->evaluate( trajectory, g );
        }
            
        if ( doing_covariant ){ metric.multiplyLowerInverse( g ); } 
    }

    return value;
}

inline void ProblemDescription::doCollisionConstraint()
{ 
    collision_constraint = true;
}
inline void ProblemDescription::dontCollisionConstraint()
{
    collision_constraint = false;
}

inline bool ProblemDescription::isCollisionConstraint() const
{
    return collision_constraint;
}

inline bool ProblemDescription::isConstrained() const 
{ 
    return !factory.empty();
}
inline bool ProblemDescription::isCovariant() const
{ 
    return doing_covariant;
}

inline void ProblemDescription::setGoalset( Constraint * goal )
{
    goalset = goal;
}

inline const Constraint * ProblemDescription::getGoalset() const
{ 
    return goalset;
}

inline int ProblemDescription::N() const { return trajectory.N(); }
inline int ProblemDescription::M() const { return trajectory.M(); }
inline int ProblemDescription::size() const { return trajectory.size(); }

inline void ProblemDescription::setUpperBounds(const MatX & upper)
{ 
    upper_bounds = upper;
}

inline void ProblemDescription::setLowerBounds(const MatX & lower)
{ 
    lower_bounds = lower;
}

inline const MatX & ProblemDescription::getUpperBounds() const 
{ 
    return upper_bounds;
}
inline const MatX & ProblemDescription::getLowerBounds() const
{
    return lower_bounds;
}
inline bool ProblemDescription::isBounded() const 
{
    if ( lower_bounds.size() == M() ) { return true; }
    if ( upper_bounds.size() == M() ) { return true; }
    return false;
}

inline const Trajectory & ProblemDescription::getTrajectory() const 
{
    return trajectory;
}
inline Trajectory & ProblemDescription::getTrajectory() 
{
    return trajectory;
}

inline const ConstraintFactory & ProblemDescription::getFactory() const
{
    return factory;
}


inline void ProblemDescription::doCovariantOptimization()
{ 
    is_covariant = true;
}
inline void ProblemDescription::dontCovariantOptimization()
{ 
    is_covariant = false;
}
inline bool ProblemDescription::isCovariantOptimization() const
{ 
    return is_covariant;
}

inline bool ProblemDescription::isSubsampled() const
{
    return trajectory.isSubsampled();
}

inline const Metric & ProblemDescription::getMetric()
{
    if ( trajectory.isSubsampled() ){
        return subsampled_metric;
    }
    return metric;
}

template <class Derived>
inline void ProblemDescription::updateTrajectory( 
                                const Eigen::MatrixBase<Derived> & delta) 
{
    if ( doing_covariant ){
        covariant_trajectory.update( delta );
    }
    trajectory.update( delta );
}

inline void ProblemDescription::updateTrajectory( const double * delta )
{
    if ( doing_covariant ){
        covariant_trajectory.update( ConstMatMap(delta, 
                                                 trajectory.N(),
                                                 trajectory.M() ) );
    }

    trajectory.update( ConstMatMap(delta, trajectory.N(), trajectory.M()));
}

//this is a local optimization specific update function, and
//  you cannot do local optimization while doing covariant optimization,
//  so the assertion is necessary
template <class Derived>
inline void ProblemDescription::updateTrajectory( 
                                const Eigen::MatrixBase<Derived> & delta,
                                int t) 
{
    debug_assert( !doing_covariant );
    trajectory.update( delta, t );
}

//this is a local optimization specific update function, and
//  you cannot do local optimization while doing covariant optimization,
//  so the assertion is necessary
inline void ProblemDescription::updateTrajectory( const double * delta,
                                                  int t )
{
    debug_assert( !doing_covariant );
    trajectory.update( 
                ConstMatMap(delta, trajectory.N(), trajectory.M()),
                t);
}



