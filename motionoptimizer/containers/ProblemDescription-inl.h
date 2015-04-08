

inline const Metric & ProblemDescription::getMetric(){
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



