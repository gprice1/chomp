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

template <class Derived>
inline void SmoothnessFunction::subsampleGradient(int N_sub, 
        const Eigen::MatrixBase<Derived> & g_sub_const)
{   
    debug_status( TAG, "getSubsampledSmoothnessFunction", "start" );
    
    //cast away the const to edit the g_sub matrix
    Eigen::MatrixBase<Derived>& g_sub = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_sub_const);

    //Subsample the g_full matrix. This uses the matrix map method
    //  to subsample the matrix, basically, we make a Matrix Map that
    //  make g_full look like the subsampled matrix that we want,
    //  and then we set g_sub to be equal to that image.
    g_sub = DynamicMatMap( g_full.data(),
                           g_sub.rows(), g_sub.cols(), 
                           DynamicStride(g_full.rows(), 2) );

    debug_status( TAG, "getSubsampledSmoothnessFunction", "end" );
    
}


template <class Derived>
inline double SmoothnessFunction::evaluate(
                    const Trajectory & trajectory,
                    const Metric & metric,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    debug_status( TAG, "evaluateGradient", "start" );
    
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
            
    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    metric.multiply( trajectory.getFullXi(), Ax);
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    if (trajectory.isSubsampled() ){
        g_full = Ax + b; 
        subsampleGradient(trajectory.N(), g );
    }else {
        g = Ax + b; 
    }
    
    debug_status( TAG, "evaluateGradient", "end" );
    
    //return the computed smoothness value of the trajectory
    return computeValue( trajectory );
    
}


