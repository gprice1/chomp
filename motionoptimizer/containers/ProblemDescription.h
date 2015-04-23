#ifndef _PROBLEM_DESCRIPTION_H_
#define _PROBLEM_DESCRIPTION_H_

#include "Trajectory.h"
#include "Gradient.h"
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
    Gradient gradient;
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
     
    const Metric & getMetric();

    bool isConstrained() const; 
    bool isCovariant() const;
    
    void setGoalset( Constraint * goal );
    const Constraint * getGoalset() const;
    
    int N() const;
    int M() const;
    int size() const;
    
    void setUpperBounds(const MatX & upper);
    void setLowerBounds(const MatX & lower);
    const MatX & getUpperBounds() const ;
    const MatX & getLowerBounds() const ;

    bool isBounded() const;

    void copyTrajectoryTo( double * data );
    void copyTrajectoryTo( std::vector<double> & data );

    void copyToTrajectory( const double * data );
    void copyToTrajectory( const std::vector<double> data );

    const Trajectory & getTrajectory() const;
    Trajectory & getTrajectory();
    
    const Gradient & getGradient() const;
    const ConstraintFactory & getFactory() const;

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
    
    void prepareData( const double * xi = NULL );
    
    
};//Class ProblemDescription



#include "ProblemDescription-inl.h"

}//namespace


#endif 
