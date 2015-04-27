#ifndef _METRIC_H_
#define _METRIC_H_

#include <Eigen/Dense>
#include "../utils/utils.h"

namespace mopt {

class Metric {

private:
    MatX L;
    Eigen::VectorXd coefficients;
    Eigen::VectorXd goalset_coefficients;
    

public:

    Metric( int n, 
            ObjectiveType otype,
            bool do_subsample=false,
            bool do_goalset=false);
    
    //a constructor that creates an empty metric
    Metric(){}

    ~Metric(){}
    void initialize( int n, 
                     ObjectiveType otype,
                     bool do_subsample=false,
                     bool do_goalset=false);

    void resize(int n ); 

    void doGoalset();
    void stopGoalset();
    
    void printL() const;

    bool isGoalset() const;
    
    int size() const; 
    int width() const;
    
    bool empty() const;
    
    //the old diagmul call
    template <class Derived>        
    void sampleNormalDistribution(
            double standard_deviation,
            Eigen::MatrixBase<Derived> const & result) const;

    //the old diagmul call
    template <class Derived1, class Derived2 >        
    void multiply( const Eigen::MatrixBase<Derived1>& x,
                   const Eigen::MatrixBase<Derived2>& Ax_const ) const ;
    
    //the old diagmul call
    template <class Derived >        
    void multiply( const Eigen::MatrixBase<Derived>& result) const ;

    //Perform a multiplication operation with the L matrix. 
    //  L is a lower triangular matrix stored in a skyline format.
    //  x_const = L * x_const;
    template <class Derived>
    void multiplyLower( const Eigen::MatrixBase<Derived>& result ) const ;

    template <class Derived1, class Derived2>
    void multiplyLower( const Eigen::MatrixBase<Derived1>& original,
                        const Eigen::MatrixBase<Derived2>& result) const;


    //Perform a multiplication operation with the L matrix. 
    //  L is a lower triangular matrix stored in a skyline format.
    //  x_const = L * x_const;
    template <class Derived >
    void multiplyLowerTranspose( 
         const Eigen::MatrixBase<Derived>& result) const;
    
    template <class Derived1, class Derived2 >
    void multiplyLowerTranspose( 
         const Eigen::MatrixBase<Derived1>& original,
         const Eigen::MatrixBase<Derived2>& result) const;


    //Perform an inverse transpose multiplication operation
    //with the L matrix. 
    //  L is a lower triangular matrix stored in a skyline format.
    //  x_const = L^-T * x_const;
    //This is used by the HMC class to generate random smooth momenta.
    template <class Derived>
    void multiplyLowerInverseTranspose(
         const Eigen::MatrixBase<Derived>& result) const ;
    
    template <class Derived1, class Derived2>
    void multiplyLowerInverseTranspose(
         const Eigen::MatrixBase<Derived1>& original,
         const Eigen::MatrixBase<Derived2>& result) const;


    //Perform an inverse multiplication operation with the L matrix. 
    //  x_const = L^-1 * x_const;
    template <class Derived>
    void multiplyLowerInverse(
         const Eigen::MatrixBase<Derived>& result) const;
    
    template <class Derived1, class Derived2>
    void multiplyLowerInverse(
         const Eigen::MatrixBase<Derived1>& original,
         const Eigen::MatrixBase<Derived2>& result) const;


    //Solve for the inverse multiplication of the A matrix using
    //  the cholesky decomposition.  
    //  L is a lower triangular matrix stored in a skyline format.
    //  x_const = L^-T * L^-1 * x_const, which is equivalent to:
    //      x_const = A^-1 * x_const
    template <class Derived>
    void solve(const Eigen::MatrixBase<Derived>& x_const) const ;

    template <class Derived1, class Derived2>
    void solve( const Eigen::MatrixBase<Derived1>& original,
                const Eigen::MatrixBase<Derived2>& result) const ;

    //Creates the b matrix which factors in the contribution of 
    //  the endpoints to both the objective function and its gradients.
    template <class Derived1, class Derived2 >
    double createBMatrix(const Eigen::MatrixBase<Derived1>& x0,
                         const Eigen::MatrixBase<Derived1>& x1,
                         const Eigen::MatrixBase<Derived2>& b_const,
                         double dt) const ;

    
    template <class Derived>
    void solveCovariantBounds( const MatX & lower, const MatX & upper,
                               const Eigen::MatrixBase<Derived> & covariant_lower,
                               const Eigen::MatrixBase<Derived> & covariant_upper ) const ;



private:
    
    /////////////////////////////////////////////////////////////////////
    //////////////////////Matrix Methods for Goalset-CHOMP///////////////
    ////These functions are called internally from their non-goalset 
    //  counterpars
   
    //this is a diag mul for goal set chomp
    template <class Derived1, class Derived2>
    void multiplyGoalset(const Eigen::MatrixBase<Derived1>& original,
                        const Eigen::MatrixBase<Derived2>& result ) const;

    template <class Derived>
    void multiplyGoalset( const Eigen::MatrixBase<Derived>& result ) const;
    
    //initializer for creating the L matrix
    void initialize( int n );

    //this is a skyline chol for goal set chomp
    void createGoalsetLMatrix( int n );
    void createLMatrix( int n );


    //This version of createBMatrix is used for goal set chomp.
    template <class Derived1, class Derived2>
    double createGoalsetBMatrix(
                         const Eigen::MatrixBase<Derived1>& x0,
                         const Eigen::MatrixBase<Derived2>& b_const,
                         double dt) const ;
    
    //methods for getting the value of the L matrix at a given 
    //  row, column index. 
    //  Since L is stored in skyline format, this is a non-trivial 
    //      operation (granted though, it is not difficult.)
    double getLValue( int row_index, int col_index ) const;
    double getCoefficientValue( int row_index, int col_index ) const;

};

#include "Metric-inl.h"

}//namespace

#endif
