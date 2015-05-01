

#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "../utils/utils.h"
#include "Metric.h"

namespace mopt {

//forward declaration of Constraint factory
class ConstraintFactory;

class Trajectory {
  

  private:
    double * data, * cached_data;
    DynamicMatMap xi;
    MatMap full_xi;

    MatX q0, q1;

    ObjectiveType objective_type;
    
    static const char* TAG;
    double dt, total_time;
    bool is_subsampled;

    //this is a scratch variable for efficiently returning
    //  ticks with border padding. Look at the function getTick. 
    mutable MatX border_tick; 
    

  public:
    Trajectory();
    
    Trajectory( const MatX & q0, const MatX & q1, int N,
                ObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & q0,
                const std::vector<double> & q1,
                int N,
                ObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const MatX & xi,
               const MatX & q0,
               const MatX & q1,
               ObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);

    //Initialization from a whole prexisting trajectory.
    Trajectory( const MatX & traj,
                ObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & traj,
                int M,
                ObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const std::vector<std::vector<double> > & traj,
               ObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);
    Trajectory(const double * traj,
               int M, int N,
               ObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);

    ~Trajectory();    
    Trajectory & operator= (const Trajectory & other);
    
    void startGoalset();
    void endGoalset();

    //restores the data to the original after 
    //  it had been replaced in a call to setData.
    void restoreData();
    
    void setData( const double * new_data );
    void setData( double * new_data );

    double * getData();
    const double * getData() const;

    void copyToData( const std::vector<double> & vec );
    void copyToData( const double * new_data );
    
    void copyDataTo( double * copy ) const ;
    void copyDataTo( std::vector<double> & vec ) const;

    const DynamicMatMap & getTraj() const;    
    void getTraj( std::vector<double> & vec ) const ;
    void getTraj( double * array ) const ;

    void getState( std::vector<double> & vec, int i ) const ;
    void getState( double * array, int i ) const ;
    
    double getDt() const;

    void subsample();
    void endSubsample();

    bool isSubsampled() const;

    double & operator ()( int i, int j );
    double const & operator()(int i, int j)const;
    
    int N() const;
    int M() const;

    int fullN() const;
    int fullM() const; 

    int rows() const;
    int cols() const; 

    int size() const;

    const MatX & getQ0() const;
    const MatX & getStartPoint() const;

    const MatX & getQ1() const;
    const MatX & getEndPoint() const;

    //get individual matrix blocks or rows
    Row row( int i );
    ConstRow row( int i ) const;

    Col col( int i );
    ConstCol col( int i ) const;

    const DynamicMatMap & getXi() const; 

    const MatMap & getFullXi() const;

    void setObjectiveType( ObjectiveType otype);
    ObjectiveType getObjectiveType() const;

    //upsample the trajectory by two times.
    void upsample();
    //upsample the trajectory until it is greater than Nmax.
    void upsampleTo( int Nmax );
    
    //upsample, respecting the given constraints.
    void constrainedUpsample( ConstraintFactory & factory,
                              double htol,
                              double hstep);
    void constrainedUpsampleTo( ConstraintFactory & factory,
                                double htol,
                                double hstep,
                                int Nmax );

    //initialize from a matrix.
    void initialize(const MatX & traj);
    void initialize(const std::vector<double> & traj, int rows);
    void initialize(const std::vector<std::vector<double> > & traj);
    void initialize(const double * traj, int rows, int cols);
    
    //initialize from two endpoints.
    void initialize( const MatX & pinit,
                     const MatX & pgoal,
                     int rows );
    void initialize( const std::vector<double> & pinit,
                     const std::vector<double> & pgoal,
                     int rows);
    void initialize( const double * pinit,
                     const double * pgoal,
                     int rows, int cols);
    
    
    //updates the trajectory via a matrix delta. Delta
    // must be the same size and shape as the trajectory,
    //  or the subsampled trajectory
    template <class Derived>
    void update( const Eigen::MatrixBase<Derived> & delta );

    //updates the trajectory at the given row, by the vector delta,
    //  delta should have the same number of columns as the trajectory.
    template <class Derived>
    void update( const Eigen::MatrixBase<Derived> & delta, int row);
    
    //a small helper function to get states that smoothly fall
    //  off the positive and negative edges of the trajectory
    MatX getTick(int tick) const;
    
    //utility functions for turning the trajectory into a 
    //  string, and then printing it.
    std::string toString() const;
    void print() const;
    
    void getCovariantTrajectory( const Metric & metric,
                                 Trajectory & other ) const;
    void getNonCovariantTrajectory( const Metric & metric,
                                    Trajectory & other) const; 
    
    void copyDataToOther(Trajectory & other ) const;

  private:
    void createInitialTrajectory(); 

  private:
    //The below functions are all inlines that are implemented
    //  inside of Trajectory-inl.h

    //copy over information 
    template <class Derived>
    void initializeData( const Eigen::MatrixBase<Derived> & traj );
    
    template <class Derived>
    void initializeData( const Eigen::MatrixBase<Derived> & pinit,
                         const Eigen::MatrixBase<Derived> & pgoal,
                         int rows );

    void remapXi( int n, int full_n, int m );
    void resizeOther( Trajectory & other ) const ;

};

#include "Trajectory-inl.h"

}//namespace

#endif 
