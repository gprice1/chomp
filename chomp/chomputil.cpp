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

#include "chomputil.h"
#include "Chomp.h"
#include <iomanip>

namespace chomp {

const char* eventTypeString(int eventtype) {
    switch (eventtype) {
    case CHOMP_INIT: return "CHOMP_INIT";
    case CHOMP_GLOBAL_ITER: return "CHOMP_GLOBAL_ITER";
    case CHOMP_LOCAL_ITER: return "CHOMP_LOCAL_ITER";
    case CHOMP_FINISH: return "CHOMP_FINISH";
    case CHOMP_TIMEOUT: return "CHOMP_TIMEOUT";
    case CHOMP_GOALSET_ITER: return "CHOMP_GOALSET_ITER";
    default: return "[INVALID]";
    }
}

ChompObserver::~ChompObserver() {}

int ChompObserver::notify(const Chomp&, 
                            ChompEventType,
                            size_t,
                            double, double, double) { 
    return 0; 
}

DebugChompObserver::~DebugChompObserver() {}

int DebugChompObserver::notify(const Chomp& c, 
                                 ChompEventType e,
                                 size_t iter,
                                 double curObjective, 
                                 double lastObjective,
                                 double constraintViolation) { 
    std::cout << "chomp debug: "
              << "event=" << eventTypeString(e) << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";

    if (std::isnan(curObjective) || std::isinf(curObjective) ||
        std::isnan(lastObjective) || std::isinf(lastObjective)) {
        return 1;
    }
    return 0;
}

ChompGradientHelper::~ChompGradientHelper() {}

ChompCollisionHelper::ChompCollisionHelper(size_t nc,
                                           size_t nw,
                                           size_t nb):
    ncspace(nc), nwkspace(nw), nbodies(nb) {}

ChompCollisionHelper::~ChompCollisionHelper() {}

ChompCollGradHelper::ChompCollGradHelper(ChompCollisionHelper* h,
                                           double g):
    chelper(h), gamma(g)
{
    dx_dq = MatX(h->nwkspace, h->ncspace);
    cgrad = MatX(h->nwkspace, 1);
}

ChompCollGradHelper::~ChompCollGradHelper() {}

double ChompCollGradHelper::addToGradient(const Chomp& c,
                                            MatX& g) {

    q1 = c.getTickBorderRepeat(-1).transpose();
    q2 = c.getTickBorderRepeat(0).transpose();

    double total = 0.0;

    for (int t=0; t<c.N; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = c.getTickBorderRepeat(t+1).transpose();

      cspace_vel = 0.5 * (q2 - q0) * c.inv_dt;        
      cspace_accel = (q0 - 2.0*q1 + q2) * (c.inv_dt * c.inv_dt);

      for (size_t u=0; u < chelper->nbodies; ++u) {

        float cost = chelper->getCost(q1, u, dx_dq, cgrad);
        if (cost > 0.0) {

          wkspace_vel = dx_dq * cspace_vel;

          //this prevents nans from propagating. Several lines below, 
          //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
          if (wkspace_vel.isZero()){ continue; }

          wkspace_accel = dx_dq * cspace_accel;
          
          float wv_norm = wkspace_vel.norm();
          wkspace_vel /= wv_norm;

          // add to total
          double scl = wv_norm * gamma / c.inv_dt;

          total += cost * scl;
          
          P = MatX::Identity(chelper->nwkspace, chelper->nwkspace)
              - (wkspace_vel * wkspace_vel.transpose());

          K = (P * wkspace_accel) / (wv_norm * wv_norm);
         

          //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
          g.row(t) += (scl * (dx_dq.transpose() *
                      (P * cgrad - cost * K)).transpose());
         

        }
      }
    }
    return total;

}

template <class Derived1, class Derived2, class Derived3>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs, // e.g. [1, -4, 6]
             const Eigen::MatrixBase<Derived2>& x,
             const Eigen::MatrixBase<Derived3>& Ax_const ) {

    assert( Ax_const.rows() == x.rows() && Ax_const.cols() == x.cols() );

    Eigen::MatrixBase<Derived3>& Ax = 
      const_cast<Eigen::MatrixBase<Derived3>&>(Ax_const);
    
    //in the case that we are doing goal chomp, n != x.rows(),
    //  so if goal chomp is being done, then, the correct value of n,
    //  should be passed in.
    int n = x.rows();

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    for (int i=0; i<x.rows(); ++i) {

      int j0 = std::max(i-o, int(0));
      int j1 = std::min(i+nc, n );

      Ax.row(i) = x.row(j0) * coeffs(j0-i+o);

      for (int j=j0+1; j<=i; ++j) {
        Ax.row(i) += x.row(j) * coeffs(j-i+o);
      }
    
      for (int j=i+1; j<j1; ++j) {
        Ax.row(i) += x.row(j) * coeffs(i-j+o);
      }

    }
}

template <class Derived1, class Derived2>
void skylineChol(int n,
                 // e.g. [-1,2] or [1,-4,6]
                 const Eigen::MatrixBase<Derived1>& coeffs,
                 Eigen::PlainObjectBase<Derived2>& L) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    L.resize(n, nc);
  
    for (int j=0; j<n; ++j) {
      
      //i1 is the current row, forwarded by the amount of coeffs.
      int i1 = std::min(j+nc, n);
      
      for (int i=j; i<i1; ++i) {

        double sum = 0;

        int k0 = std::max(0,i-o);

        for (int k=k0; k<j; ++k) {
          sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
        }

        if (i == j) {
          L(j,o) = sqrt(coeffs(o) - sum);
        } else {
          L(i,j-i+o) = (coeffs(j-i+o) - sum) / L(j,o);
        }
      }
    }

}
//////////////////////////////////////////////////////////////////////
//This is used by the HMC class to generate random smooth momenta.
template <class Derived1, class Derived2>
void skylineCholMultiplyInverseTranspose(
                               const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const)
{

    const int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    const int nc = L.cols();
    const int o = nc-1;

    for (int i=n-1; i>=0; --i) {
      const int j1 = std::min(i+nc, n);
      for (int j=i+1; j<j1; ++j) {
        x.row(i) -= L(j, i-j+o) * x.row(j); // here j > i so col < row
      }
      x.row(i) /= L(i,o);
    }
}

template <class Derived1, class Derived2>
void skylineCholMultiplyInverse(const Eigen::MatrixBase<Derived1>& L,
                                const Eigen::MatrixBase<Derived2>& x_const)
{
    int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    int nc = L.cols();
    int o = nc-1;

    for (int i=0; i<n; ++i) {
      int j0 = std::max(0, i-o);
      for (int j=j0; j<i; ++j) {
        x.row(i) -= L(i,j-i+o)*x.row(j); // here j < i so col < row
      }
      x.row(i) /= L(i,o);
    }
}


template <class Derived1, class Derived2>
void skylineCholSolve(const Eigen::MatrixBase<Derived1>& L,
                      const Eigen::MatrixBase<Derived2>& x_const)
{

    int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    int nc = L.cols();
    int o = nc-1;

    for (int i=0; i<n; ++i) {
      int j0 = std::max(0, i-o);
      for (int j=j0; j<i; ++j) {
        x.row(i) -= L(i,j-i+o)*x.row(j); // here j < i so col < row
      }
      x.row(i) /= L(i,o);
    }

    for (int i=n-1; i>=0; --i) {
      int j1 = std::min(i+nc, n);
      for (int j=i+1; j<j1; ++j) {
        x.row(i) -= L(j,i-j+o) * x.row(j); // here j > i so col < row
      }
      x.row(i) /= L(i,o);
    }

}


template <class Derived1, class Derived2>
void skylineCholSolveMulti( const Eigen::MatrixBase<Derived1>& L, 
                            const Eigen::MatrixBase<Derived2>& xx_const)
{
    int n = L.rows();
    int m = xx_const.rows() / n;
    assert(xx_const.rows() == m*n);

    Eigen::MatrixBase<Derived2>& xx = 
      const_cast<Eigen::MatrixBase<Derived2>&>(xx_const);

    for (int i=0; i<m; ++i) {
        skylineCholSolve(L, xx.block(n*i, 0, n, xx.cols()));
    }
    
}


template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived2>& x1,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt)
{

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    assert(b_const.rows() == n);
    assert(b_const.cols() == x0.cols());
    assert(b_const.cols() == x1.cols());

    Eigen::MatrixBase<Derived3>& b = 
      const_cast<Eigen::MatrixBase<Derived3>&>(b_const);

    b.setZero();

    double c = 0;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      int i1 = n-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        int t1 = -t0;
        b.row(i0) += coeffs(j)*getPos(x0, t0*dt);
        b.row(i1) += coeffs(j)*getPos(x1, t1*dt);
      }
      c += mydot(b.row(i0), b.row(i0));
      if (i0 != i1) { c += mydot(b.row(i1), b.row(i1)); }
    }

    return 0.5*c;
}
  
  
//this is a diag mul for goal set chomp
template <class Derived1, class Derived2, class Derived3, class Derived4>
void diagMul(const Eigen::MatrixBase<Derived1>& coeffs,
             const Eigen::MatrixBase<Derived2>& gs_coeffs,
             const Eigen::MatrixBase<Derived3>& x,
             const Eigen::MatrixBase<Derived4>& Ax_const ){

    assert( Ax_const.rows() == x.rows() && Ax_const.cols() == x.cols() );

    Eigen::MatrixBase<Derived4>& Ax = 
      const_cast<Eigen::MatrixBase<Derived4>&>(Ax_const);
    
    //in the case that we are doing goal chomp, n != x.rows(),
    //  so if goal chomp is being done, then, the correct value of n,
    //  should be passed in.
    int n = x.rows();

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    const int start_gs = x.rows() - gs_coeffs.rows();

    for (int i=0; i<x.rows(); ++i) {

      int j0 = std::max(i-o, int(0));
      int j1 = std::min(i+nc, n );
      
      Ax.row(i) = x.row(j0) * coeffs(j0-i+o);

      for (int j=j0+1; j<=i; ++j) {
        double coeff;
        if ( i >= start_gs && j >= start_gs ){ 
            coeff = gs_coeffs( i - start_gs, j - start_gs );
        }else { coeff = coeffs(j-i+o); }

        Ax.row(i) += x.row(j) * coeff;
      }
    
      for (int j=i+1; j<j1; ++j) {
        double coeff;
        if ( i >= start_gs && j >= start_gs ){ 
            coeff = gs_coeffs( i - start_gs, j - start_gs );
        }else { coeff = coeffs(i-j+o); }

        Ax.row(i) += x.row(j) * coeff;
      }
    }
}

//this is a skyline chol for goal set chomp
template <class Derived1, class Derived2, class Derived3>
void skylineChol(int n,
                 const Eigen::MatrixBase<Derived1>& coeffs,
                 const Eigen::MatrixBase<Derived2>& gs_coeffs,
                 Eigen::PlainObjectBase<Derived3>& L) {

    int nc = coeffs.size();
    int o = nc-1;

    const int start_gs = n - gs_coeffs.rows();

    L.resize(n, nc);
  
    for (int j=0; j<n; ++j) {
      
      //i1 is the current row, forwarded by the amount of coeffs.
      int i1 = std::min(j+nc, n);
      
      for (int i=j; i<i1; ++i) {

        double sum = 0;

        int k0 = std::max(0,i-o);

        for (int k=k0; k<j; ++k) {
          sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
        }
        
        if (i == j) {
          double coeff;
          if ( i >= start_gs ){
            coeff = gs_coeffs( i - start_gs , i - start_gs);
          } else { coeff = coeffs(o); }

          L(j,o) = sqrt(coeff - sum);

        } else {
          double coeff;
          if ( i >= start_gs && j >= start_gs){
            coeff = gs_coeffs( i - start_gs, j - start_gs);
          } else { coeff = coeffs(j-i+o); }

          L(i,j-i+o) = (coeff - sum) / L(j,o);

        }
      }
    }
}



//This version of createBMatrix is used for goal set chomp.
template <class Derived1, class Derived2, class Derived3>
double createBMatrix(int n, 
                     const Eigen::MatrixBase<Derived1>& coeffs,
                     const Eigen::MatrixBase<Derived2>& x0,
                     const Eigen::MatrixBase<Derived3>& b_const,
                     double dt) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    assert(b_const.rows() == n);
    assert(b_const.cols() == x0.cols());

    Eigen::MatrixBase<Derived3>& b = 
      const_cast<Eigen::MatrixBase<Derived3>&>(b_const);

    b.setZero();

    double c = 0;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        b.row(i0) += coeffs(j)*getPos(x0, t0*dt);
      }
      c += mydot(b.row(i0), b.row(i0));
    }

    return 0.5*c;
}


}//namespace
