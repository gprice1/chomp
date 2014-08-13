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

#include "ConstraintFactory.h"
#include "Constraint.h"

namespace chomp {

  ConstraintFactory::~ConstraintFactory(){
    clearConstraints();
  }

  void ConstraintFactory::clearConstraints() {

    for ( std::vector<Constraint*>::iterator it = constraints.begin(); 
          it != constraints.end();
          it ++ )
    {
        if (*it){ delete *it; }
    }
    constraints.clear();
  }

  void ConstraintFactory::getAll(size_t total) {
    clearConstraints();
    
    //fill up the constraints vector.
    constraints.resize( total );
    for (size_t i=0; i<total; i++){
      constraints[i] = getConstraint(i, total);
    }
  }
  

  size_t ConstraintFactory::numOutput()
  {
    //compute the total dimensionality of the constraints
    int constraint_dims = 0;
    for ( std::vector<Constraint*>::iterator it = constraints.begin();
          it != constraints.end();
          it ++ )
    {
        //if there is a constraint, get its dimensions.
        if ( *it ) { constraint_dims += (*it)->numOutputs(); }
    }

    return constraint_dims;
  }



  void ConstraintFactory::evaluate( const MatX& xi, 
                                    MatX& h_tot, 
                                    MatX& H_tot, 
                                    int step)
  {

    size_t DoF = xi.cols();

    assert(size_t(xi.rows()) == constraints.size());

    size_t numCons = 0;

    size_t timesteps = 0;

    for (size_t t=0; t<constraints.size(); t+=step) {
      Constraint* c = constraints[t];
      if (c) {
        numCons += c->numOutputs();
      }
      ++timesteps;
    }

    if (size_t(h_tot.rows()) != numCons || size_t(h_tot.cols()) != 1) {
      h_tot.resize(numCons,1);
    }

    if (size_t(H_tot.rows()) != numCons || size_t(H_tot.cols()) != DoF*timesteps) {
      H_tot.resize(numCons,DoF*timesteps);
    }

    H_tot.setZero(); 

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints

    for (size_t t=0, i=0; t<constraints.size(); t+=step, ++i) {

      Constraint* c = constraints[t];

      //get individual h and H from qt
      if (!c || c->numOutputs()==0){ continue; }
      c->evaluateConstraints(xi.row(t), h, H);

      if (h.rows() == 0) { 
        assert(H.rows() == 0);
        continue;
      }
    
      //stick all the h and H into h_tot and H_tot
      assert((size_t)H.cols()==DoF);
      assert(h.cols() == 1);
      assert((size_t)H.rows()<=constraints[t]->numOutputs());
      assert(H.rows() == h.rows());

      for (size_t r = 0; r<(size_t)h.rows(); r++){
        h_tot(row) = h(r);
        for (size_t j = 0; j<DoF; j++){
          assert( row < size_t(H_tot.rows()) );
          assert( j*timesteps+i < size_t(H_tot.cols()) );
          assert( j < size_t(H.cols()) );
          H_tot(row, j*timesteps+i) = H(r,j);
        }
        row++;
      }

    }  

    assert( row <= (size_t)numCons );
    if ( row < numCons ) {
      h_tot.conservativeResize(row, 1);
      H_tot.conservativeResize(row, DoF*timesteps);
    }

    assert(h_tot.rows() == H_tot.rows());
    assert(h_tot.cols() == 1);
    assert(size_t(H_tot.cols()) == DoF*timesteps);

  }


  void ConstraintFactory::evaluate( ConstMatMap& xi, 
                                    MatMap& h_tot, 
                                    MatMap& H_tot)

  {
    size_t DoF = xi.cols();
    size_t N = xi.rows();
    H_tot.setZero();    
    
    assert(N == constraints.size());

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints
    
    for (size_t t=0, i=0; t<constraints.size(); t++, ++i) {

      Constraint* c = constraints[t];

      //get individual h and H from qt
      if (!c || c->numOutputs()==0){ continue; }

      c->evaluateConstraints(xi.row(t), h, H);

      //stick all the h and H into h_tot and H_tot
      assert(H.cols() == (int)DoF);
      assert(h.cols() == 1);
      assert(H.rows() == (int)constraints[t]->numOutputs());
      assert(H.rows() == h.rows());
      

      for (size_t r = 0; r<(size_t)h.rows(); r++, row++){
        h_tot(row) = h(r);

        for (size_t j = 0; j<DoF; j++){
          assert( H_tot.cols() > (int)row   );
          assert( H_tot.rows() > (int)(j*N + i) );
          
          H_tot(j*N + i, row) = H(r,j);
        }
      }
    }  
  }

  void ConstraintFactory::evaluate( ConstMatMap& xi, 
                                    MatMap& h_tot)
  {
    
    assert(size_t(xi.rows()) == constraints.size());

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints
    
    for (size_t t=0, i=0; t<constraints.size(); t++, ++i) {

      Constraint* c = constraints[t];

      //get individual h and H from qt
      if (!c || c->numOutputs()==0){ continue; }

      c->evaluateConstraints(xi.row(t), h, H);
      
      for (size_t r = 0; r<(size_t)h.rows(); r++, row++){
        h_tot(row) = h(r);
      }
    }  
  }

  void ConstraintFactory::evaluate( unsigned constraint_dim,
                                    double* result,
                                    unsigned n_by_m,
                                    const double * x,
                                    double* grad )
  {
    int N = constraints.size();
    int M = n_by_m / N;

    assert( constraint_dim == numOutput() );

    //put the data into matrix maps.
    ConstMatMap xi( x, N, M);
    MatMap h_total(result, constraint_dim, 1);

    if ( grad == NULL ){
        evaluate( xi, h_total );
    }else{
        MatMap H_total( grad, n_by_m, constraint_dim);
    
        //evaluate the data.
        evaluate( xi, h_total, H_total);
    }
  }

} //namespace
