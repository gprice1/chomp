#ifndef _PROBLEM_DESCRIPTION_H_
#define _PROBLEM_DESCRIPTION_H_

#include "Trajectory.h"
#include "ChompGradient.h"
#include "ConstraintFactory.h"
#include "Metric.h"

#ifdef DO_TIMING
    #include "../utils/timer.h"
    #include <sstream>
    #define TIMER_START( x )        timer.start( x )   
    #define TIMER_STOP( x )         timer.stop( x )    
    #define TIMER_GET_TOTAL( x )    timer.getTotal( x )
#else 
    #define TIMER_START( x )     (void)0 
    #define TIMER_STOP( x )      (void)0
    #define TIMER_GET_TOTAL( x ) (void)0
#endif 

namespace chomp {
    
class ProblemDescription {
    
    //this is used to give the MotionOptimizer 
    //  class the necessary machinery to pass the
    //  Trajectory back and forht.
    friend class MotionOptimizer;
    
private:
    Trajectory trajectory;
    ChompGradient gradient;
    ConstraintFactory factory;
    
    //in the case that we are doing covariant optimization,
    //  this trajectory holds the covariant state
    Trajectory covariant_trajectory;

    MatX lower_bounds, upper_bounds;
    
    Constraint * goalset;

    bool ok_to_sample, use_goalset;
    bool is_covariant, doing_covariant;

    MatX g_full;
    
    //Defines for timing things
#ifdef DO_TIMING
    Timer timer;
    static const std::string gradient_timer;
    static const std::string copy_timer;
    static const std::string constraint_timer;
#endif
    
public:
    
    ProblemDescription();
    ~ProblemDescription();

    void upsample();
    void subsample();
    void stopSubsample();
 
    const Metric & getMetric();

    inline bool isConstrained() const { return !factory.empty(); }
    inline bool isCovariant() const { return doing_covariant; }
    
    inline void setGoalset( Constraint * goal ){goalset = goal;}
    inline const Constraint * getGoalset() const { return goalset; }
    
    inline int N(){ return trajectory.N(); }
    inline int M(){ return trajectory.M(); }
    inline int size(){ return trajectory.size(); }
    
    inline void setUpperBounds(const MatX & upper){ upper_bounds = upper;}
    inline void setLowerBounds(const MatX & lower){ lower_bounds = lower;}
    inline const MatX & getUpperBounds(){ return upper_bounds;}
    inline const MatX & getLowerBounds(){ return lower_bounds;}

    inline void copyTrajectoryTo( double * data )
                { trajectory.copyDataTo( data ); }
    inline void copyTrajectoryTo( std::vector<double> & data )
                { trajectory.copyDataTo( data ); }

    void copyToTrajectory( const double * data );
    void copyToTrajectory( const std::vector<double> data );

    inline const Trajectory & getTrajectory() const {return trajectory;}
    inline const ChompGradient & getGradient() const {return gradient;}
    inline const ConstraintFactory & getFactory() const {return factory;}

    double evaluateGradient( MatX & g );
    double evaluateGradient( const double * xi, double * g );

    double evaluateConstraint( MatX & h );
    double evaluateConstraint( MatX & h, MatX & H );
    double evaluateConstraint( const double * xi,
                                     double * h,
                                     double * H = NULL);
    
    int getConstraintDims();
    
    bool evaluateConstraint( MatX & h_t, MatX & H_t, int t );

    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta);
    void updateTrajectory( const double * delta );
    
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta, int t);
    void updateTrajectory( const double * delta, int t );

    double evaluateObjective();
    double evaluateObjective( const double * xi );
    
    inline void doCovariantOptimization(){ is_covariant = true; }
    inline void dontCovariantOptimization(){ is_covariant = false; }
    inline bool isCovariantOptimization(){ return is_covariant; }
    
    inline bool isSubsampled(){ return trajectory.isSubsampled(); }


    void getTimes( 
        std::vector< std::pair<std::string, double> > & times ) const ;
    
    void printTimes( bool verbose=true ) const ;

    std::string getTimesString( bool verbose=true ) ;
    
private:

    void prepareRun();
    void prepareSample();
    void prepareCovariant( const double * xi = NULL );
    
    void startGoalset();
    void stopGoalset();
    
};//Class ProblemDescription


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
        prepareCovariant();
    }
    trajectory.update( delta );
}

inline void ProblemDescription::updateTrajectory( const double * delta )
{
    if ( doing_covariant ){
        covariant_trajectory.update( ConstMatMap(delta, 
                                                 trajectory.N(),
                                                 trajectory.M() ) );
        prepareCovariant();
    }

    trajectory.update( ConstMatMap(delta, trajectory.N(), trajectory.M()));
}

template <class Derived>
inline void ProblemDescription::updateTrajectory( 
                                const Eigen::MatrixBase<Derived> & delta,
                                int t) 
{
    debug_assert( !doing_covariant );
    trajectory.update( delta, t );
}

inline void ProblemDescription::updateTrajectory( const double * delta,
                                                  int t )
{
    debug_assert( !doing_covariant );
    trajectory.update( 
                ConstMatMap(delta, trajectory.N(), trajectory.M()),
                t);
}


}//namespace


#endif 
