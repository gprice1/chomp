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

#include "MotionOptimizer.h"
#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef MZ_HAVE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#endif

using namespace mopt;


class CircleConstraint : public Constraint {
  public:

    virtual size_t numOutputs(){
        return 1;
    }

    virtual void evaluateConstraints( const MatX& qt,
                                      MatX& h, 
                                      MatX& H)
    {

        assert((qt.cols()==1 && qt.rows()==2) || 
               (qt.cols()==2 && qt.rows()==1));

        h.resize(1,1);
        H.resize(1,2);
        
        h(0) = mydot(qt,qt) - 4; 
        
        for(int i=0; i<2; i++){
            H(i) = 2*qt(i);
        }

    }

};


#ifdef MZ_HAVE_CAIRO

class PdfEmitter: public DebugObserver {
  public:

    int frequency;
    int count;

    double sz, scl, m;
    cairo_surface_t* surface;
    cairo_t* cr;

    const char* filename;

    std::ostringstream ostream;

    bool dump_data_to_file;
  

    PdfEmitter(const char* f, bool dump_data): 
        frequency(0), count(0),
        filename(f), dump_data_to_file( dump_data )
    {

        sz = 4*72; 
        scl = sz/10;
        m = 0;

        surface = cairo_pdf_surface_create(filename, sz+2*m, sz+2*m);

        cr = cairo_create(surface);
        ostream << f << "{";
    }
  

    virtual ~PdfEmitter()
    {
        cairo_surface_destroy(surface);
        cairo_destroy(cr);

        std::cout << "wrote " << filename << "\n\n";
    }

    void appendInfoToFile( const std::string & filename )
    {
        if (filename.size() > 0){
            std::ofstream myfile;
            myfile.open (filename.c_str(), std::ios::app );
            
            myfile << ostream.str() << "}\n";

            myfile.close();
        }
    }

    virtual int notify( const OptimizerBase & chomper, 
                        EventType event,
                        size_t iter,
                        double curObjective,
                        double lastObjective,
                        double hmag) 
    {
        if ( dump_data_to_file ){
            ostream << "[" << chomper.problem.getTimesString()
                    << ", iter:" << iter 
                    << ", constraint:" << hmag
                    << ", objective:" << curObjective
                    <<  "], ";
        }else {
            DebugObserver::notify(chomper, event, iter, 
                                   curObjective, lastObjective, hmag);

        }
        
        if ( frequency < 0 ) { return 0; }
        if ( !( (event == INIT)   ||
                (event == FINISH) ||
                (frequency > 0 && iter % frequency == 0 ) ))
        {
            return 0;
        }

        if (count++) { 
            cairo_show_page(cr);
        }

        double x0 = -4;
        double y0 = 6;

#define MAPX(x) (((x)-x0)*scl+m)
#define MAPY(y) ((y0-(y))*scl+m)

        cairo_set_line_width(cr, 1);
        cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
        cairo_arc(cr, MAPX(0), MAPY(0), 2*scl, 0, 2*M_PI);
        cairo_stroke(cr);
        
        const Trajectory & trajectory = chomper.problem.getTrajectory();
        
        for (int i=-1; i<=trajectory.fullN(); ++i) {
            MatX pi;
            if (i < 0) { 
                pi = trajectory.getQ0();
            } else if (i >= trajectory.fullN()) {
                pi = trajectory.getQ1();
            } else {
                pi = trajectory.getFullXi().row(i);
            }
            double u = double(i+1) / (trajectory.fullN()+1);
            cairo_set_source_rgb(cr, 1-u, 0, u);
            cairo_arc(cr, MAPX(pi(0)), MAPY(pi(1)), 2, 0, 2*M_PI);
            cairo_stroke(cr);
        }

        cairo_set_line_width(cr, 2);
        cairo_rectangle(cr, MAPX(-4), MAPY(6), 10*scl, 10*scl);
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_stroke(cr);

          
        return 0;

    }

};

#endif


void usage(int status) {
    std::ostream& ostr = status ? std::cerr : std::cout;
    ostr <<
      "usage: circle_demo OPTIONS\n"
      "\n"
      "OPTIONS:\n"
      "\n"
      "  -l, --algorithm          Algorithm used for optimization\n"
      "  -n, --num-initial        Number of steps for initial trajectory\n"
      "  -t, --num-final          Minimum number of steps for final trajectory\n"
      "  -v, --no-local           Disable local smoothing\n"
      "  -g, --no-global          Disable global smoothing\n"
      "  -m, --no-multigrid       Disable multigrid computation\n"
      "  -e, --error-tol          Relative error tolerance\n"
      "  -a, --alpha              Step size for CHOMP\n"
      "  -d, --dump               Dump recorded data to a file\n"
      "  -o, --objective          Either 'accel' or 'vel' depending on what you want to optimizer for\n"
      "  -p, --pdf                Output PDF's\n"
      "  -k, --covariance         Do covariant optimization\n"
      
      "      --help               See this message.\n";
    exit(status);
}
    

int main(int argc, char** argv) {

    int N = 15;
    int Nmax = 127;
    bool do_covariant = false;
    OptimizationAlgorithm alg = NONE;

    bool doMultigrid = true;
    bool doGlobalSmooth = true;
    bool doLocalSmooth = true;
    int doPDF = -2;
    double alpha = 0.05;
    double errorTol = 1e-7;
    ObjectiveType objective = MINIMIZE_VELOCITY;
    std::string filename = "";
    bool dump_data = false;

    const struct option long_options[] = {
        { "algorithm",         required_argument, 0, 'l' },
        { "num-initial",       required_argument, 0, 'n' },
        { "num-final",         required_argument, 0, 't' },
        { "error-tol",         required_argument, 0, 'e' },
        { "alpha",             required_argument, 0, 'a' },
        { "objective",         required_argument, 0, 'o' },
        { "pdf",               required_argument, 0, 'p' },
        { "dump_to_file",      required_argument, 0, 'd' },
        { "covariance",        no_argument,       0, 'k' },
        { "no-multigrid",      no_argument,       0, 'm' },
        { "no-global",         no_argument,       0, 'g' },
        { "no-local",          no_argument,       0, 'v' },
        { "help",              no_argument,       0, 'h' },
        { 0,                   0,                 0,  0  }
    };

    const char* short_options = "l:n:t:e:a:o:p:d:kmdgvh";
    int opt, option_index;

    while ( (opt = getopt_long(argc, argv, short_options, 
                               long_options, &option_index) ) != -1 ) {

        switch (opt) {
        case 'l': alg = algorithmFromString( optarg ); break;
        case 'k': do_covariant = true; break;
        case 'n': N = atoi(optarg); break;
        case 't': Nmax = atoi(optarg); break;
        case 'e': errorTol = atof(optarg); break; 
        case 'o':
            if ( std::strcmp( optarg, "accel") == 0 ){
                objective = MINIMIZE_ACCELERATION;
            }
            else if ( std::strcmp( optarg, "vel") == 0 ){
                objective = MINIMIZE_VELOCITY;
            }else {
                std::cout << "Input '" << opt << "' is an invalid"
                          << " option for '--objective'\n";
                usage(1);
            }
            break;
        case 'a': alpha = atof( optarg ); break;
        case 'p': doPDF = atoi( optarg ); break;
        case 'd': 
            filename = std::string( optarg );
            dump_data = true;
            break;
        case 'm': doMultigrid = false; break;
        case 'g': doGlobalSmooth = false; break;
        case 'v': doLocalSmooth = false; break;
        case 'h': usage(0); break;
        default:
            std::cout << "opt: " << opt << "\n";
            usage(1);
            break;
        }

    }

    if (!doMultigrid) { 
        N = Nmax; 
        doLocalSmooth = false;
    }
  
    std::cout << "about to optimize with settings:\n";
    std::cout << "  init n:           " << N << "\n";
    std::cout << "  final n:          " << Nmax << "\n";
    std::cout << "  step size:        " << alpha << "\n";
    std::cout << "  error tol:        " << errorTol << "\n";
    std::cout << "  multigrid:        " << (doMultigrid ? "ON":"OFF") << "\n";
    std::cout << "  local smoothing:  " << (doLocalSmooth ? "ON":"OFF") << "\n";
    std::cout << "  global smoothing: " << (doGlobalSmooth ? "ON":"OFF") << "\n\n";
    
    
    DebugObserver obs;
    MotionOptimizer chomper( &obs );

    CircleConstraint c;
    chomper.addConstraint( &c, 0.25, 0.75 );
    chomper.setMaxIterations( 200 );
    
    MatX q0(1,2), q1(1,2);
    q0 << -3, 5;
    q1 << 5, -3;

    chomper.getTrajectory().initialize( q0, q1, N );  
    chomper.setNMax( Nmax );
    
    chomper.getTrajectory().setObjectiveType( MINIMIZE_VELOCITY );
    
    chomper.dontSubsample();

    if (do_covariant ){ chomper.doCovariantOptimization(); }
    
    chomper.setAlpha( alpha );

    chomper.getTrajectory().setObjectiveType( objective );
    
    chomper.setAlgorithm( alg ); 
    
#ifdef MZ_HAVE_CAIRO
    PdfEmitter* pobs = NULL;

    if (doPDF >= -1 ) {

        char pdf_name[1024];

        snprintf(pdf_name, 1024, "circle_%s_%04d_%04d_%f_%s_%s.pdf",
                 algorithmToString( alg ).c_str(),
                 N, Nmax, alpha,
                 objective == MINIMIZE_VELOCITY ? "vel" : "accel",
                 do_covariant ? "covariant" : "non-covariant" );

        pobs = new PdfEmitter(pdf_name, dump_data);

        pobs->frequency = doPDF;
        chomper.setObserver( pobs );
    }
#endif

    chomper.solve();
  
    
#ifdef MZ_HAVE_CAIRO
    if ( pobs ){
        pobs->appendInfoToFile( filename );
        delete pobs;
    }
#endif

  return 0;
}
