

#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "../utils/utils.h"
#include "Metric.h"

namespace chomp{

class Trajectory {
  public:

    double * data, * cached_data;
    DynamicMatMap xi;
    MatMap full_xi;

    MatX q0, q1;

    ChompObjectiveType objective_type;
    
    static const char* TAG;
  private:
    double dt, total_time;
    bool is_subsampled;

  public:
    Trajectory();
    
    Trajectory( const MatX & q0, const MatX & q1, int N,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & q0,
                const std::vector<double> & q1,
                int N,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const MatX & xi,
               const MatX & q0,
               const MatX & q1,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);

    //Initialization from a whole prexisting trajectory.
    Trajectory( const MatX & traj,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & traj,
                int M,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const std::vector<std::vector<double> > & traj,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);
    Trajectory(const double * traj,
               int M, int N,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
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

    double * getData(){ return data; }
    const double * getData() const { return data; } 
    
    void copyToData( const std::vector<double> & vec );
    void copyToData( const double * new_data );
    
    void copyDataTo( double * copy ) const ;
    void copyDataTo( std::vector<double> & vec ) const;

    const DynamicMatMap & getTraj()const { return xi; }
    //ConstDynamicMatMap  & getTraj()      { return xi; }
    
    void getTraj( std::vector<double> & vec ) const ;
    void getTraj( double * array ) const ;

    void getState( std::vector<double> & vec, int i ) const ;
    void getState( double * array, int i ) const ;
    
    inline double getDt() const { return dt; }

    void subsample();
    void endSubsample();

    bool isSubsampled() const { return is_subsampled; }
    
    inline double & operator ()( int i, int j ){
        return xi( i, j );
    }
    inline double const & operator()(int i, int j)const {
        return xi( i, j );
    }

    
    inline int N() const { return xi.rows(); }
    inline int M() const { return xi.cols(); }

    inline int fullN() const { return full_xi.rows(); }
    inline int fullM() const { return full_xi.cols(); }

    inline int rows() const { return xi.rows();}
    inline int cols() const { return xi.cols();}

    inline int size() const { return xi.size();}

    inline const MatX & getQ0() const { return q0; }
    inline const MatX & getStartPoint() const { return q0; }

    inline const MatX & getQ1() const { return q1; };
    inline const MatX & getEndPoint() const { return q1; }

    //get individual matrix blocks or rows
    inline Row row( int i ) { return xi.row(i); }
    inline ConstRow row( int i ) const{ return xi.row(i); }

    inline Col col( int i ){  return xi.col(i); }
    inline ConstCol col( int i ) const { return xi.col(i); }

    inline const DynamicMatMap & getXi() const { return xi; } 
    //inline DynamicMatMap & getXi() { return xi; } 

    inline const MatMap & getFullXi() const { return full_xi; } 
    //inline MatMap & getFullXi() { return full_xi; } 

    inline void setObjectiveType( ChompObjectiveType otype){
        objective_type = otype;
    }
    inline ChompObjectiveType getObjectiveType() const {
        return objective_type;
    }

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
    inline MatX getTick(int tick) const;
    
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

    //copy over information 
    template <class Derived>
    void initializeData( const Eigen::MatrixBase<Derived> & traj );
    
    template <class Derived>
    void initializeData( const Eigen::MatrixBase<Derived> & pinit,
                         const Eigen::MatrixBase<Derived> & pgoal,
                         int rows );

    void remapXi( int n, int full_n, int m );

};


inline void Trajectory::remapXi( int n, int full_n, int m ){
    
    //find out which inner stride to use.
    //  if full_n is not equal to n, then we are subsampling,
    //  and the correct stride is 2, otherwise, the correct stride is 1.
    const int inner_stride = full_n == n ? 1 : 2;
    
    new (&xi) DynamicMatMap( data, n, m, 
                             DynamicStride(full_n, inner_stride) );
    new (&full_xi) MatMap( data, full_n, m);
}

//definitions of the inline functions.

//initialize the data from a preexistant matrix
template <class Derived>
inline void Trajectory::initializeData( 
                               const Eigen::MatrixBase<Derived> & traj )
{

    q0 = traj.row(0);
    q1 = traj.row( traj.rows() - 1 );

    const int n = traj.rows() - 2;
    const int m = traj.cols();

    //if there is previously existant data,
    //  delete it.
    if (cached_data){ 
        delete cached_data;
        cached_data = NULL;
    } else if (data){
        delete data;
    }

    //allocate the data
    data = new double[ n * m ];

    //remap the data
    remapXi( n, n, m );

    //copy over the data from the trajectory
    full_xi = traj.block( 1, 0, n, m );

    dt = total_time / xi.rows()+1;
    
}

//definitions of the inline functions.
template <class Derived>
inline void Trajectory::initializeData( 
                               const Eigen::MatrixBase<Derived> & pinit,
                               const Eigen::MatrixBase<Derived> & pgoal,
                               int n )
{
    q0 = pinit;
    q1 = pgoal;

    const int m = q0.size();

    //if there is previously existant data,
    //  delete it.
    if (cached_data){ 
        delete cached_data;
        cached_data = NULL;
    } else if (data){
        delete data;
    }
    
    //allocate the data
    data = new double[ n * m ];

    //remap the data
    remapXi( n, n, m );

    //copy over the data from the trajectory
    createInitialTrajectory();

    dt = total_time / xi.rows()+1;
    
}


inline MatX Trajectory::getTick(int tick) const
{
    //if the tick is negative, get a state that falls off the
    //  the edge of the trajectory
    if (tick < 0) { return getPos(q0, (tick+1)*dt); }
    
    //if the tick is larger than the number of states,
    //  get a state that falls off the positive edge of the
    //  trajectory
    if (tick >= N() ) { return getPos(q1, (tick-N())*dt);}

    //if the tick corresponds to a state in the trajectory,
    //  return the corresponding state.
    return xi.row( tick );
} 


template <class Derived>
inline void Trajectory::update(const Eigen::MatrixBase<Derived> & delta)
{
    xi -= delta;
}

template <class Derived>
inline void Trajectory::update(const Eigen::MatrixBase<Derived> & delta,
                               int index)
{
    xi.row( index ) -= delta;
}

}//namespace

#endif 
