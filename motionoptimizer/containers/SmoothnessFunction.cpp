/*
* Copyright (c) 2008-2014, Matt Zucker
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


#include "SmoothnessFunction.h"


namespace mopt {

const char* SmoothnessFunction::TAG = "SmoothnessFunction";

void SmoothnessFunction::prepareRun(const Trajectory & trajectory,
                                    const Metric & metric)
{
    debug_status( TAG, "prepareRun", "start" );
    
    Ax.resize( trajectory.fullN(), trajectory.M() );
    b.resize( trajectory.fullN(), trajectory.M() );
    
    //get the b matrix, and get its contribution to the
    //  objective function.
    c = metric.createBMatrix(trajectory.getStartPoint(),
                             trajectory.getEndPoint(),
                             b,
                             trajectory.getDt());
        
    debug_status( TAG, "prepareRun", "end" );
}

double SmoothnessFunction::evaluate( const Trajectory & trajectory,
                                     const Metric & metric )
{
    metric.multiply( trajectory.getFullXi(), Ax);
    return computeValue( trajectory );
}


double SmoothnessFunction::computeValue(const Trajectory & trajectory) const
{
    return 0.5 * mydot(trajectory.getFullXi(), Ax) + 
           mydot(trajectory.getFullXi(), b) +
           c;
}

}// namespace

