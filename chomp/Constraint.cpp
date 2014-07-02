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

#include "Constraint.h"

namespace chomp {

  Constraint::~Constraint() {}

  NullConstraint::~NullConstraint() {}

  size_t NullConstraint::numOutputs() {
    return 0;
  }

  void NullConstraint::evaluateConstraints(const MatX& qt, 
                                           MatX& h, 
                                           MatX& H){
  }

  ConstantConstraint::~ConstantConstraint() {}



  ConstantConstraint::ConstantConstraint(const std::vector<size_t>& consIndex, 
                                         const std::vector<double>& consValue){
    // Check to see if the input is valid
    assert(consIndex.size()==consValue.size()); 
    // Make sure at most 1 constraint per DoF
    index = consIndex;
    std::sort(index.begin(), index.end());
    std::vector<size_t>::iterator it;
    it = std::adjacent_find(index.begin(), index.end());
    assert(it==index.end());
    // Assign values
    index = consIndex;
    value = consValue;
    numCons = consIndex.size();
  }

  size_t ConstantConstraint::numOutputs() {
    return numCons;
  }

  void ConstantConstraint::evaluateConstraints(const MatX& qt, 
                                               MatX& h, 
                                               MatX& H)

  {

    assert(qt.rows() == 1 || qt.cols() == 1);

    size_t DoF = std::max(qt.rows(), qt.cols());

    if (size_t(h.rows()) != numCons || h.cols() != 1) {
      h.resize(numCons, 1);
    }
    
    if (size_t(H.rows()) != numCons || size_t(H.cols()) != DoF) {
      H.resize(numCons, DoF);
    }

    H.setZero();

    for(size_t i=0; i<numCons; i++){
      assert(index[i] < DoF);
      h(i) = qt(index[i])-value[i];
      H(i,index[i]) = 1;
    }

  }
    

  ///////////////////// GENERALIZED TSR CONSTRAINTS ////////////////
    
    size_t TSRConstraint::numOutputs(){
        return _dim_constraint;
    }

    //a basic constructor
    TSRConstraint::TSRConstraint( MatX & pose_0_w,
                                  MatX & Bw,
                                  MatX & pose_w_e ) :
        _pose_0_w( pose_0_w ),
        _Bw( Bw ),
        _pose_w_e( pose_w_e )
    {
        
        // check the dimensions of the matrices to make sure that 
        // they are of the correct shape
        assert( _pose_0_w.cols() == 1 ); 
        assert( _pose_0_w.rows() == 7 );
        
        assert( _pose_w_e.cols() == 1 );
        assert( _pose_w_e.rows() == 7 );
        
        assert( _Bw.cols() == 2 );
        assert( _Bw.rows() == 6 );

        calculateDimensionality();
        calculateInverses();
    }

    void TSRConstraint::calculateDimensionality()
    {
        //get the dimension of the volume by counting all of the
        // dimensions of Bw that have non-zero width. 
        _dim_volume = 0;
        for ( int i = 0; i < 6; i ++ ){

            //if the two values in each row are not the same, then the
            // dimension has width.
            if( _Bw( i, 0 ) != _Bw(i, 1) ){
                
                //for each dimension that has width, increase the 
                // dimension of the volume by 1. 
                _dim_volume ++;
            }
        }
        

        //get the dimension of the constraint by counting up all of the 
        // dimensions that do NOT go from infinity to negative infinity
        _dim_constraint = 0;
        for ( int i = 0; i < 6; i ++ ){

            if ( -HUGE_VAL < _Bw(i, 0 ) || _Bw(i, 1 ) < HUGE_VAL ){

                //store the index of all dimensions that are constrained
                _dimension_id.push_back( i );

                // if a dimension is constrained, count it
                _dim_constraint ++;
            }
        }
    }

    void TSRConstraint::calculateInverses(){
        _pose_w_0 = _pose_0_w.inverse();
        _pose_e_w = _pose_w_e.inverse();
    }


    void TSRConstraint::endeffectorToTSRFrame( const MatX & qt,
                                                     MatX & xyzrpy){
        
        MatX pose_0_endeffector;

        forwardKinematics( qt, pose_0_endeffector );
        // transform the pos into the the correct frame
        
        
        const MatX pose_e_0 = _pose_e_w * _pose_w_0;

        // TODO : finish this function!!

    }
    //Evaluate the constraints for Chompification
    void TSRConstraint::evaluateConstraints(const MatX& qt, 
                                                          MatX& h, 
                                                          MatX& H)
    {
        //the pos is equivalent to the position of the 
        // end-effector in the robot frame. 
        MatX xyzrpy;
        
        endeffectorToTSRFrame( qt, xyzrpy );

        size_t DoF = std::max(qt.rows(), qt.cols());
        
        //format h (constraint value vector) and H (jacobian matrix).
        if (size_t(h.rows()) != _dim_constraint || h.cols() != 1)
        {
          h.resize(_dim_constraint, 1);
        }
        
        if (size_t(H.rows()) != _dim_constraint || size_t(H.cols()) != DoF)         {
          H.resize(_dim_constraint, DoF);
        }

        H.setZero();
        
        for ( int i = 0; i < _dim_constraint; i ++ )
        {
            int dim = _dimension_id[i];
            
            //if the robot's position goes over the TSR's upper bound:
            if ( xyzrpy(dim) > _Bw(dim, 1) ){
                h(i) = xyzrpy( dim ) - _Bw(dim, 1);

                //TODO Add Jacobian.
            }
            //if the robot's position goes below the TSR's lower bound:
            else if ( xyzrpy(dim) < _Bw( dim, 0 ) ){
                h(i) = xyzrpy( dim ) - _Bw(dim, 0);

                //TODO : Add Jacobian information
            }
            //if the robot's position is inside of the TSR's bounds:
            // TODO : Should the Jacobian be set to 0 here?
            //  This could be bad because there are disconintuities in
            //  the derivative of the constraint function
            else {
                h(i) = 0.0;
                H(i, dim) = 0.0;
            }
        }
    }


}

