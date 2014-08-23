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

#ifndef _CONSTRAINT_FACTORY_H_
#define _CONSTRAINT_FACTORY_H_

#include "chomputil.h"
#include <vector>

namespace chomp {

  class Constraint;
  

  class ConstraintFactory {
  public:
    
    std::vector<Constraint*> constraints;
    
    //base constructor, if constraints is non-empty, it deletes all
    //  of the constraints which are non-NULL.
    virtual ~ConstraintFactory(); 

    void clearConstraints();

    virtual Constraint* getConstraint(size_t t, size_t total) =0;
    
    void getAll(size_t total);
    
    size_t numOutput();

    virtual void evaluate(const MatX& xi, 
                          MatX& h_tot, 
                          MatX& H_tot, 
                          int step=1);

    virtual void evaluate(ConstMatMap& xi, 
                          MatMap& h_tot, 
                          MatMap& H_tot);
    virtual void evaluate(ConstMatMap& xi, 
                          MatMap& h_tot);
    
    virtual void evaluate( unsigned constraint_dim,
                           double* result,
                           unsigned n_by_m,
                           const double * x,
                           double* grad );

    static void NLoptConstraint(unsigned constraint_dim,
                                double* result,
                                unsigned n_by_m,
                                const double * x,
                                double* grad,
                                void *data) 
    {
        reinterpret_cast<ConstraintFactory*>(data) 
            ->evaluate( constraint_dim, result, n_by_m, x, grad);
    }

  };

}


#endif
