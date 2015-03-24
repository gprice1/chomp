#ifndef _PROBLEM_DESCRIPTION_H_
#define _PROBLEM_DESCRIPTION_H_

#include "Trajectory.h"
#include "ChompGradient.h"
#include "ConstraintFactory.h"

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

    MatX lower_bounds, upper_bounds;
    
    Constraint * goalset;

    //This is true, if prepareRun has been called since the
    //  last change to the trajectory, or goalset has been made.
    bool ok_to_run, use_goalset;

public:
    
    ProblemDescription();

    void upsample();
    void subsample();
    void stopSubsample();
 
    const MatX & getLMatrix();

    inline bool isConstrained() const { return !factory.empty(); }

    inline void setGoalset( Constraint * goal )
            {goalset = goal; ok_to_run = false;}
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
    /*
    bool evaluateConstraint( const double * x_t, 
                                   double * h,
                                   double * H = NULL);
    */

    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta);
    void updateTrajectory( const double * delta );
    
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta, int t);
    void updateTrajectory( const double * delta, int t );

    double evaluateObjective();
    double evaluateObjective( const double * xi );
    

private:

    void prepareRun();
    void startGoalset();
    void stopGoalset();
    


};//Class ProblemDescription


inline const MatX & ProblemDescription::getLMatrix(){
    if ( trajectory.isSubsampled() ){
        return gradient.getLsubMatrix();
    }
    return gradient.getLMatrix();
}

template <class Derived>
inline void ProblemDescription::updateTrajectory( 
                                const Eigen::MatrixBase<Derived> & delta) 
{
    trajectory.update( delta );
}

inline void ProblemDescription::updateTrajectory( const double * delta )
{
    trajectory.update( ConstMatMap(delta, trajectory.N(), trajectory.M()));
}

template <class Derived>
inline void ProblemDescription::updateTrajectory( 
                                const Eigen::MatrixBase<Derived> & delta,
                                int t) 
{
    trajectory.update( delta, t );
}

inline void ProblemDescription::updateTrajectory( const double * delta,
                                                  int t )
{
    trajectory.update( 
                ConstMatMap(delta, trajectory.N(), trajectory.M()),
                t);
}


}//namespace


#endif 
