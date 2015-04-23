
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

inline const Gradient & ProblemDescription::getGradient() const 
{
    return gradient;
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
        return gradient.getSubsampledMetric();
    }
    return gradient.getMetric();
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



