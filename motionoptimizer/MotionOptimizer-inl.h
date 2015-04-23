//simple getters and setters, all of which are inline
inline void MotionOptimizer::setNMax( int n_max ){ N_max = n_max; }
inline int  MotionOptimizer::getNMax( ) const { return N_max; }

inline void MotionOptimizer::setGoalset( Constraint * goal)
{
    problem.setGoalset(goal);
}
inline const Constraint * MotionOptimizer::getGoalset() const 
{
    return problem.getGoalset();
}

inline void MotionOptimizer::setTimeoutSeconds( double s )
{ 
    timeout_seconds = s;
}
inline double MotionOptimizer::getTimeoutSeconds() const
{ 
    return timeout_seconds;
}

inline void MotionOptimizer::setMaxIterations( size_t max )
{ 
    max_iterations = max; 
}
inline size_t MotionOptimizer::getMaxIterations() const
{ 
    return max_iterations;
}

inline void MotionOptimizer::setFunctionTolerance( double tol )
{ 
    obstol = tol;
}
inline double MotionOptimizer::getFunctionTolerance() const
{ 
    return obstol;
}

inline void MotionOptimizer::setAlgorithm(OptimizationAlgorithm a)
{ 
    algorithm1 = a;
}

inline OptimizationAlgorithm MotionOptimizer::getAlgorithm() const
{ 
    return algorithm1;
}

inline void MotionOptimizer::setAlgorithm1(OptimizationAlgorithm a1)
{ 
    algorithm1 = a1;
}
inline void MotionOptimizer::setAlgorithm2(OptimizationAlgorithm a2)
{ 
    algorithm2 = a2;
}

inline OptimizationAlgorithm MotionOptimizer::getAlgorithm1() const
{ 
    return algorithm1;
}
inline OptimizationAlgorithm MotionOptimizer::getAlgorithm2() const
{ 
    return algorithm2;
}

inline void MotionOptimizer::dontSubsample()
{ 
    do_subsample = false;
}
inline void MotionOptimizer::doSubsample()
{ 
    do_subsample = true;
}
inline void MotionOptimizer::setSubsample( bool subsample )
{ 
    do_subsample = subsample;
}

inline void MotionOptimizer::setAlpha( double a )
{ 
    alpha = a;
}
inline double MotionOptimizer::getAlpha() const
{ 
    return alpha;
}

inline void MotionOptimizer::doFullGlobalAtFinal()
{ 
    full_global_at_final = true;
}
inline void MotionOptimizer::dontFullGlobalAtFinal()
{ 
    full_global_at_final = false;
}
inline bool MotionOptimizer::getFullGlobalAtFinal() const
{ 
    return full_global_at_final;
}

inline Trajectory & MotionOptimizer::getTrajectory()
{ 
    return problem.trajectory;
}

inline const Trajectory & MotionOptimizer::getTrajectory() const 
{ 
    return problem.getTrajectory();
}
inline void MotionOptimizer::setTrajectory( const Trajectory & trajectory )
{ 
    problem.trajectory = trajectory;
}

inline void MotionOptimizer::setCollisionFunction( 
                                CollisionFunction * coll_func)
{
    problem.collision_function = coll_func;
}
inline const CollisionFunction * MotionOptimizer::getCollisionFunction()
const
{
    return problem.collision_function;
}

inline void MotionOptimizer::setObserver( Observer * obs )
{ 
    observer = obs;
}

inline Observer * MotionOptimizer::getObserver()
{ 
    return observer;
}
inline const Observer * MotionOptimizer::getObserver() const 
{ 
    return observer;
}

inline void MotionOptimizer::doCovariantOptimization()
{ 
    problem.doCovariantOptimization();
}
inline void MotionOptimizer::dontCovariantOptimization()
{ 
    problem.dontCovariantOptimization();
}
inline bool MotionOptimizer::isCovariantOptimization() const
{ 
    return problem.isCovariantOptimization();
}

