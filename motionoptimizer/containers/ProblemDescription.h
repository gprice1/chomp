#ifndef _PROBLEM_DESCRIPTION_H_
#define _PROBLEM_DESCRIPTION_H_

#include "Trajectory.h"
#include "SmoothnessFunction.h"
#include "CollisionFunction.h"
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

namespace mopt {
    
class ProblemDescription {
    
    //this is used to give the MotionOptimizer 
    //  class the necessary machinery to pass the
    //  Trajectory back and forht.
    friend class MotionOptimizer;
    
private:
    Trajectory trajectory;
    SmoothnessFunction smoothness_function;
    CollisionFunction * collision_function;
    ConstraintFactory factory;

    Metric metric, subsampled_metric;
    
    //in the case that we are doing covariant optimization,
    //  this trajectory holds the covariant state
    Trajectory covariant_trajectory;

    MatX lower_bounds, upper_bounds;
    
    Constraint * goalset;

    bool use_goalset;
    bool is_covariant, doing_covariant;
    bool collision_constraint;

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
     
    template <class Derived>
    double evaluateObjective( const Eigen::MatrixBase<Derived> & g );
    double evaluateObjective( const double * xi=NULL, double * g=NULL );
    
    double evaluateConstraint( MatX & h );
    double evaluateConstraint( MatX & h,
                               MatX & H);
    double evaluateConstraint( const double * xi,
                                     double * h,
                                     double * H = NULL);
    bool evaluateConstraint( MatX & h_t, MatX & H_t, int t );


    //Functions with the trajectory
    int N() const;
    int M() const;
    int size() const;

    void copyTrajectoryTo( double * data );
    void copyTrajectoryTo( std::vector<double> & data );

    void copyToTrajectory( const double * data );
    void copyToTrajectory( const std::vector<double> data );

    const Trajectory & getTrajectory() const;
    Trajectory & getTrajectory();

    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta);
    void updateTrajectory( const double * delta );
    
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta, int t);
    void updateTrajectory( const double * delta, int t );

    
    //Functions with constraints:
    bool isConstrained() const; 
    bool isCovariant() const;

    int getConstraintDims();

    void doCollisionConstraint();
    void dontCollisionConstraint();
    bool isCollisionConstraint() const ;
    
    const ConstraintFactory & getFactory() const;
    
    void setGoalset( Constraint * goal );
    const Constraint * getGoalset() const;


    const Metric & getMetric();

    //Set stuff with bounds    
    void setUpperBounds(const MatX & upper);
    void setLowerBounds(const MatX & lower);
    const MatX & getUpperBounds() const ;
    const MatX & getLowerBounds() const ;

    bool isBounded() const;
    
    void doCovariantOptimization();
    void dontCovariantOptimization();
    bool isCovariantOptimization() const; 
    
    bool isSubsampled() const;
    
    void getTimes( 
        std::vector< std::pair<std::string, double> > & times ) const ;
    
    void printTimes( bool verbose=true ) const ;

    std::string getTimesString( bool verbose=true ) ;

    void getFullBounds( std::vector< double > & lower,
                        std::vector< double > & upper );
    
private:

    void prepareRun( bool subsample);
    void endRun();
    
    //prepareData must have been called before this function.
    template <class Derived>
    double computeObjective( const Eigen::MatrixBase<Derived> & g);

    void prepareData( const double * xi = NULL );
    
    
};//Class ProblemDescription



#include "ProblemDescription-inl.h"

}//namespace


#endif 
