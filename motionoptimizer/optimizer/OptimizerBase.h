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

#ifndef _OPTIMIZER_BASE_H_
#define _OPTIMIZER_BASE_H_

#include "../utils/utils.h"
#include "../utils/Observer.h"
#include "../containers/ProblemDescription.h"

namespace mopt {

class OptimizerBase{
   
    //this allows the observer to peek at the 
    //  problem member variables
    friend class Observer;
    
  public:
    
    ProblemDescription & problem;
    Observer * observer;

    double obstol, timeout_seconds;
    size_t max_iter;
    
    double last_objective, current_objective, constraint_magnitude;
    size_t current_iteration;

    
    OptimizerBase( ProblemDescription & problem,
                   Observer * observer = NULL,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1));

    virtual ~OptimizerBase(){}

    virtual void solve()=0;

    //notify the observer
    int notify(EventType event) const;
    
};

}//namespace

#endif
