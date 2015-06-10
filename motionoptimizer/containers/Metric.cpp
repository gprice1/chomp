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


#include "Metric.h"

namespace mopt {

Metric::Metric( int n ,
                ObjectiveType otype,
                bool do_subsample,
                bool do_goalset)
{
    initialize( n, otype, do_subsample, do_goalset );
}

bool Metric::isGoalset() const
{
    return goalset_coefficients.size() > 0;
};
    
int Metric::size() const { return L.rows(); }

int Metric::width()const { return L.cols(); }

bool Metric::empty() const { return L.size() > 0; }

double Metric::getLValue( int row_index, int col_index ) const 
{
    const int offset = coefficients.size() - 1;
    return L(row_index, offset + col_index - row_index );
}


double Metric::getCoefficientValue( int row_index, int col_index ) const 
{
    const int offset = coefficients.size() - 1;
    return coefficients( offset + col_index - row_index );
}

void Metric::initialize( int n ,
                         ObjectiveType otype,
                         bool do_subsample,
                         bool do_goalset )
{
    if (otype == MINIMIZE_VELOCITY) {
        
        if ( do_subsample ){
            coefficients.resize(1);
            coefficients << 2;
        }else {
            coefficients.resize(2);
            coefficients << -1, 2;
        }
        
        if (do_goalset ){
            goalset_coefficients.resize(1,1);
            goalset_coefficients << 1;
        }
    } else {
        if (do_subsample ){
            coefficients.resize(2);
            coefficients<< 1, 6;
        }else {
            coefficients.resize(3);
            coefficients << 1, -4, 6;
        }
        
        if (do_goalset){ 
            goalset_coefficients.resize(2,2);
            goalset_coefficients << 6, -3,
                                   -3,  2 ;
        }
    }


    if ( do_goalset ){ createGoalsetLMatrix( n ); }
    else { createLMatrix( n ); }
}

void Metric::resize( int n )
{
    
    assert( coefficients.size() > 0 );
    
    if (goalset_coefficients.size() > 0 ){
        createGoalsetLMatrix( n );
    }else {
        createLMatrix( n );
    }

}

void Metric::doGoalset()
{
   
    if (coefficients.size() == 3 ){ 
        goalset_coefficients.resize(2,2);
        goalset_coefficients << 6, -3,
                               -3,  2 ;
    } else if ( coefficients.size() == 2 ){
        goalset_coefficients.resize(1,1);
        goalset_coefficients << 1;
    }
    
    createGoalsetLMatrix( size() );

}

void Metric::stopGoalset()
{
    //do nothing if we are not in goalset configuration already.
    if ( goalset_coefficients.size() == 0 ){ return; }
    
    goalset_coefficients.resize(0,0);

    createLMatrix( size() );
}

void Metric::createLMatrix( int n )
{
    L.resize(n, coefficients.size());
    
    const int nc = width();
    const int o = nc-1;
  
    for (int j=0; j<n; ++j) {
        
        const int i1 = std::min(j+nc, n);
        for (int i=j; i<i1; ++i) {

            double sum = 0;
            const int k0 = std::max(0,i-o);

            for (int k=k0; k<j; ++k) {
                sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
            }

            if (i == j) {
                L(j,o) = sqrt(coefficients(o) - sum);
            } else {
                L(i,j-i+o) = (coefficients(j-i+o) - sum) / L(j,o);
            }
        }
    }

}

//this is a skyline chol for goal set chomp
void Metric::createGoalsetLMatrix(int n ) 
{

    L.resize(n, coefficients.size());

    int nc = coefficients.size();
    int o = nc-1;

    const int start_gs = n - goalset_coefficients.rows();

  
    for (int j=0; j<n; ++j) {
      
        //i1 is the current row, forwarded by the amount of coefficients.
        int i1 = std::min( j+width() , n);
      
        for (int i=j; i<i1; ++i) {

            double sum = 0;

            int k0 = std::max(0, i-o);

            for (int k=k0; k<j; ++k) {
                sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
            }
            
            double coeff;
            if ( i >= start_gs && j >= start_gs){
                coeff = goalset_coefficients( i - start_gs,
                                              j - start_gs);
            } else {
                coeff = getCoefficientValue(i,j);
            }

            if ( i == j ){
                L(j,o) = sqrt(coeff - sum);
            }else {
                L(i,j-i+o) = (coeff - sum) / L(j,o);
            }
        }
    }
}

void Metric::printL() const
{
    std::cout << "L matrix: \n" << L << std::endl; 
}






}// namespace
