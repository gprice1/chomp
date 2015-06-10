/*
* Copyright (c) 2008-2015, Matt Zucker and Temple Price
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "mzcommon/mersenne.h"
#include "mzcommon/gauss.h"


//the old diagmul call
template <class Derived>        
void Metric::sampleNormalDistribution( 
                double standard_deviation,
                Eigen::MatrixBase<Derived> const & result_const) const
{

    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived>& result = 
        const_cast<Eigen::MatrixBase<Derived>&>(result_const);
    
 
    for (int i = size()-1; i>=0 ; --i) {
        const int j1 = std::min( i+width(), size() );

        for ( int j = 0; j < result.cols(); ++j ){
            result(i,j) = gauss_ziggurat( standard_deviation );
        }
        
        for (int j=i+1 ; j<j1 ; ++j) {
            result.row(i) -= getLValue(j, i) * result.row(j);
        }
        
        result.row(i) /= getLValue( i, i );
    }
}

template <class Derived1, class Derived2>
void Metric::multiply(const Eigen::MatrixBase<Derived1>& original,
                      const Eigen::MatrixBase<Derived2>& result_const)const 
{
    
    if ( isGoalset() ){
        multiplyGoalset( original, result_const );
        return;
    }
    
    assert( result_const.rows() == original.rows() &&
            result_const.cols() == original.cols() );

    Eigen::MatrixBase<Derived2>& result = 
      const_cast<Eigen::MatrixBase<Derived2>&>(result_const);
    

    for (int i=0; i < size(); ++i){

        const int j0 = std::max(i - width() + 1, 0);
        const int j1 = std::min(i + width()    , size() );

        result.row(i) = original.row(j0) * getCoefficientValue(i, j0);

        for (int j=j0+1; j<=i; ++j) {
            result.row(i) += original.row(j) * getCoefficientValue(i, j);
        }
    
        for (int j = i+1; j<j1; ++j) {
            result.row(i) += original.row(j) * getCoefficientValue(j, i);
        }
    }
}

template <class Derived >        
void Metric::multiply( const Eigen::MatrixBase<Derived>& result) const 
{
    //TODO make this more intelligent.
    Eigen::MatrixXd original = result;
    multiply( original, result );

}

//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived>
void Metric::multiplyLower(
     const Eigen::MatrixBase<Derived>& result_const) const 
{
    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived>& result = 
        const_cast<Eigen::MatrixBase<Derived>&>(result_const);

    
    for (int i = size()-1; i >= 0; --i) {
        
        result.row(i) *= getLValue( i, i );
        
        const int j1 = std::max( 0, i - width() );
        for (int j = i-1; j > j1; --j) {
            result.row(i) += getLValue(i, j) * result.row( j ); //
        }
    }
}

//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived1, class Derived2>
void Metric::multiplyLower(
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const) const
{
    assert(original.rows() == size() );
    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived2>& result = 
        const_cast<Eigen::MatrixBase<Derived2>&>(result_const);

    
    for (int i = size()-1; i >= 0; --i) {
        
        result.row(i) = original.row(i) * getLValue( i, i );
        
        const int j1 = std::max( 0, i - width());
        for (int j = i-1; j > j1; --j) {
            result.row(i) += getLValue(i, j) * result.row( j );
        }
    }
}


//Perform a multiplication operation with the L matrix. 
//  L is a lower triangular matrix stored in a skyline format.
//  x_const = L * x_const;
template <class Derived>
void Metric::multiplyLowerTranspose(
     const Eigen::MatrixBase<Derived>& result_const) const
{
    assert( result_const.rows() == size() );

    Eigen::MatrixBase<Derived>& result = 
        const_cast<Eigen::MatrixBase<Derived>&>(result_const);

     
    for (int i=0 ; i < size() ; ++i )
    {
        result.row(i) *= getLValue(i, i);

        const int j1 = std::min( size(), i+width() );
        for (int j = i+1; j < j1; ++j) {
            result.row(i) += getLValue( j, i) * result.row( j ); //
        }
    }
}

template <class Derived1, class Derived2>
void Metric::multiplyLowerTranspose(
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const) const
{
    assert( result_const.rows() == size() );

    Eigen::MatrixBase<Derived2>& result = 
        const_cast<Eigen::MatrixBase<Derived2>&>(result_const);

     
    // i is the column index 
    // j indexes into the correct row of xi
    for (int i=0 ; i < size() ; ++i )
    {
        result.row(i) = original.row(i) * getLValue(i, i);

        const int j1 = std::min( size(), i+width());
        for (int j = i+1; j < j1; ++j) {
            result.row(i) += getLValue(j, i) * original.row( j ); //
        }
    }
}

//////////////////////////////////////////////////////////////////////
//This is used by the HMC class to generate random smooth momenta.
template <class Derived>
void Metric::multiplyLowerInverseTranspose(
     const Eigen::MatrixBase<Derived>& result_const) const 
{ 

    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived>& result = 
        const_cast<Eigen::MatrixBase<Derived>&>(result_const);
    

    for (int i = size()-1; i>=0 ; --i) {
        const int j1 = std::min( i+width(), size() );
        
        for (int j=i+1 ; j<j1 ; ++j) {
            result.row(i) -= getLValue(j, i) * result.row(j);
        }
        
        result.row(i) /= getLValue( i, i );
    }
}


template <class Derived1, class Derived2>
void Metric::multiplyLowerInverseTranspose(
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const) const
{
    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived2>& result = 
        const_cast<Eigen::MatrixBase<Derived2>&>(result_const);
    
    
    for (int i = size()-1; i >= 0; --i) {
        const int j1 = std::min( i+width(), size());
        
        result.row(i) = original.row(i);

        for (int j=i+1; j<j1; ++j) {
            result.row(i) -= getLValue(j, i) * result.row(j);
        }
        
        result.row(i) /= getLValue( i, i );
    }
}



template <class Derived>
void Metric::multiplyLowerInverse( 
     const Eigen::MatrixBase<Derived>& result_const) const
{
    assert(result_const.rows() == size() );

    Eigen::MatrixBase<Derived>& result = 
      const_cast<Eigen::MatrixBase<Derived>&>(result_const);


    for (int i=0; i < size(); ++i) {
        int j0 = std::max(0, i - width() + 1);
        
        for (int j=j0; j<i; ++j) {
            result.row(i) -= getLValue(i, j) * result.row(j);
        }
        
        result.row(i) /= getLValue(i,i);
    }
}

template <class Derived1, class Derived2>
void Metric::multiplyLowerInverse( 
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const) const
{
    assert(result_const.rows() == size() );
    assert(original.rows()     == size() );

    Eigen::MatrixBase<Derived2>& result = 
      const_cast<Eigen::MatrixBase<Derived2>&>(result_const);


    for (int i = 0; i < size(); ++i) {
        int j0 = std::max(0, i - width() + 1);
        
        result.row(i) = original.row(i);
        
        for (int j=j0; j<i; ++j) {
            result.row(i) -= getLValue(i, j) * result.row(j);
        }
        
        result.row(i) /= getLValue(i,i);
    }
}

template <class Derived>
void Metric::solve( const Eigen::MatrixBase<Derived>& result ) const
{
    multiplyLowerInverse( result );
    multiplyLowerInverseTranspose( result );
}

template <class Derived1, class Derived2>
void Metric::solve(
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const) const
{
    multiplyLowerInverse( original, result_const );
    multiplyLowerInverseTranspose( original, result_const );
}


template <class Derived1, class Derived2>
double Metric::createBMatrix(
                     const Eigen::MatrixBase<Derived1>& x0,
                     const Eigen::MatrixBase<Derived1>& x1,
                     const Eigen::MatrixBase<Derived2>& b_const,
                     double dt) const
{
    
    if ( isGoalset() ){ return createGoalsetBMatrix( x0, b_const, dt ); }
    
    assert(b_const.rows() == size());
    assert(b_const.cols() == x0.cols());
    assert(b_const.cols() == x1.cols());


    Eigen::MatrixBase<Derived2>& b = 
      const_cast<Eigen::MatrixBase<Derived2>&>(b_const);

    b.setZero();
    double c = 0;

    for (int i0=0; i0 < width()-1; ++i0) {
        const int nn = width() - i0 - 1;
        const int i1 = size() - i0 - 1;

        for (int j=0; j<nn; ++j) {
            // so when j=o, we want t0=0, and when j<o we want t0<0
            //    when j=o, we want t1=0, and when j<o we want t1>0
            int t0 = j - width() + 1;
            int t1 = -t0;
        
            b.row(i0) += coefficients(j) * getPos(x0, t0*dt);
            b.row(i1) += coefficients(j) * getPos(x1, t1*dt);
        }
        c += mydot(b.row(i0), b.row(i0));

        if (i0 != i1) { c += mydot(b.row(i1), b.row(i1)); }
    }

    return 0.5*c;
}
  
  
//this is a diag mul for goal set chomp
template <class Derived1, class Derived2 >
void Metric::multiplyGoalset(
     const Eigen::MatrixBase<Derived1>& original,
     const Eigen::MatrixBase<Derived2>& result_const ) const
{
    assert( result_const.rows() == original.rows() && 
            result_const.cols() == original.cols() );
    
    Eigen::MatrixBase<Derived2>& result = 
      const_cast<Eigen::MatrixBase<Derived2>&>(result_const);
    
    const int n = original.rows();
    
    const int start_gs = n - goalset_coefficients.rows();

    for (int i=0; i < n; ++i) {

        int j0 = std::max(i - width() + 1, 0 );
        int j1 = std::min(i + width()  , n );
      
        result.row(i) = original.row(j0) * getCoefficientValue(i, j0);

        for (int j=j0+1; j<=i; ++j) {
            
            if ( i >= start_gs && j >= start_gs ){ 
                result.row(i) += original.row(j) *
                    goalset_coefficients( i - start_gs, j - start_gs );
            }else { 
                result.row(i) += original.row(j) * 
                                 getCoefficientValue(i,j); 
            }
        }
    
        for (int j=i+1; j<j1; ++j) {
            
            if ( i >= start_gs && j >= start_gs ){ 
                result.row(i) += original.row(j) *
                    goalset_coefficients( i - start_gs, j - start_gs );
            }else { 
                result.row(i) += original.row(j) *
                                 getCoefficientValue(j,i);
            }
        }
    }

}


template <class Derived>
void Metric::multiplyGoalset( 
     const Eigen::MatrixBase<Derived>& result ) const
{
    //TODO maybe make this more efficient
    Eigen::MatrixXd original = result;
    multiplyGoalset( original, result );
}

//This version of createBMatrix is used for goal set chomp.
template <class Derived1, class Derived2>
double Metric::createGoalsetBMatrix(
                     const Eigen::MatrixBase<Derived1>& x0,
                     const Eigen::MatrixBase<Derived2>& b_const,
                     double dt) const
{

    
    assert(b_const.rows() == size() );
    assert(b_const.cols() == x0.cols());

    Eigen::MatrixBase<Derived2>& b = 
      const_cast<Eigen::MatrixBase<Derived2>&>(b_const);

    b.setZero();

    double c = 0;
    
    const int nc = width();
    const int o = nc-1;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        b.row(i0) += coefficients(j)*getPos(x0, t0*dt);
      }
      c += mydot(b.row(i0), b.row(i0));
    }

    return 0.5*c;
}

template <class Derived>
void Metric::solveCovariantBounds( const MatX & lower, const MatX & upper,
                const Eigen::MatrixBase<Derived> & covariant_lower_const,
                const Eigen::MatrixBase<Derived> & covariant_upper_const )
const
{
    
    assert( upper.size() > 0 && lower.size() > 0 );
    assert( upper.size() == covariant_lower_const.cols() );
    assert( upper.size() == covariant_upper_const.cols() );
    assert( covariant_lower_const.rows() == L.rows() );
    assert( covariant_upper_const.rows() == L.rows() );
     
    Eigen::MatrixBase<Derived> & covariant_upper = 
         const_cast<Eigen::MatrixBase<Derived>&>(covariant_upper_const);
    Eigen::MatrixBase<Derived> & covariant_lower = 
         const_cast<Eigen::MatrixBase<Derived>&>(covariant_lower_const);

    covariant_upper.setZero();
    covariant_lower.setZero();

    for (int i=0 ; i < size() ; ++i )
    {

        const int j1 = std::min( size(), i+width() );
        for (int j = i+1; j < j1; ++j) {
            const double value = getLValue( j, i);
            if (value > 0 ){
                //if value is positive, we want the smallest possible product
                //  in lower, which is a small number (the bounds in lower)
                //  times the positive value 
                //  Likewise, we want to add the largest possible value into
                //  covariant_upper, so we multiply the positive value by the 
                //  large number in the upper bounds vector.
                covariant_lower.row( i ) += value * lower;
                covariant_upper.row( i ) += value * upper;
            } else if (value < 0 ){
                covariant_lower.row( i ) += value * upper;
                covariant_upper.row( i ) += value * lower;
            }
        }

        //TODO find out why this is necessary
        if ( covariant_lower.row( i ) == covariant_upper.row( i ) ){
            covariant_lower.row(i) = lower;
            covariant_upper.row(i) = upper;
        }
    }
}

