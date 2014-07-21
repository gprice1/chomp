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

#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include "chomputil.h"
#include "../mzcommon/Transform3.h"
#include <vector>

namespace chomp {

  typedef Transform3_t< double > Transform;

  class Constraint {
  public:
  
    virtual ~Constraint();

    virtual size_t numOutputs() =0;

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H) =0;

  };

  class NullConstraint: public Constraint {
  public:

    virtual ~NullConstraint();
  
    virtual size_t numOutputs();

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H);

  };

  class ConstantConstraint: public Constraint {
  public:

    size_t numCons;
    std::vector<size_t> index;
    std::vector<double> value;

    ConstantConstraint(const std::vector<size_t>& consIndex, 
                       const std::vector<double>& consValue);
  
    virtual ~ConstantConstraint();
  
    virtual size_t numOutputs();

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H);

  };

  class TSRConstraint: public Constraint {
  public:
    
    //The main features of the TSR: 
    //    _pose_0_w : the origin to TSR transform
    //    _BW       : the constraint matrix
    //    _pose_w_e : the TSR to end-effector transform
    const Transform _pose_0_w;
    MatX _Bw;
    const Transform _pose_w_e;

    const Transform _pose_w_0, _pose_e_w;
    
    bool _getRotationJacobian, _getTranslationJacobian;

    //dimensionality of different features:
    // Dimensionality of the volume defined by the TSR,
    // Dimensionality of constraint surface. 
    int _dim_volume, _dim_constraint;

    std::vector<int> _dimension_id; // a vector that holds the indices 
                                   //   (in Bw) that correspond to
                                   //   constrained dimensions

    TSRConstraint( Transform & pose_0_w, 
                   MatX & Bw, 
                   Transform & pose_w_e ); 
    

    virtual ~TSRConstraint(){};
  
    virtual size_t numOutputs();
    
    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H);
    
  private:  
    //these are useful helper functions that determine various things.
    // They should not be overwritten.
    inline void calculateDimensionality();
    inline void endeffectorToTSRFrame( const Transform & pose_ee,
                                       double xyzrpy[]);

    //this function takes in a robot state, qt, and returns the position of
    // the relevant end-effector in the world frame. This is equivalent
    //  to the transformation from the end-effector frame to the world
    //  frame.
    virtual void forwardKinematics( const MatX& qt, Transform & pos ) = 0;
    
    //this takes in a state and it gets the jacobian 
    virtual void computeJacobian( const MatX& qt,
                                  const Transform & pose_world_ee,
                                  MatX & jacobian,
                                  std::vector< int > & active_dims
                                  ) = 0;


  };
}

#endif

