
#ifndef _TIMER_H_
#define _TIMER_H_

#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>

#include <unistd.h>
#include <time.h>
#include <sys/time.h>

namespace mopt {

class Timer{

  public:
    bool print;
  
  private:
    class single_timer{
      public:
        double tic, toc;
        double total, elapsed;
        unsigned int count;
        bool isStopped;
        
        single_timer() : 
                total( 0.0 ), count( 0 ),
                isStopped( true ){}
    };
    typedef std::pair< std::string, single_timer> KeyValuePair;
    typedef boost::unordered_map< std::string, single_timer > Map;

    Map timers; 
    

    single_timer * getTimer( std::string & name, bool post_error=true ){
        
        Map::iterator it = timers.find( name );

        //if the object exists, return the thing.
        if ( it != timers.end() ){ return &(it->second); }
        else if ( post_error ){
            std::cout << "No timer exists by the name: " 
                      << name << "\n";
            return NULL;
        }

        //if the timer does not exist, create and return it.
        KeyValuePair new_timer;
        new_timer.first = name;

        std::pair< Map::iterator, bool> inserted_element = 
                       timers.insert( new_timer );

        return &(inserted_element.first->second);
    }


  public:
    
    
    Timer() : print(false){}

    void coutElapsed( std::string name )
    {
        single_timer * t = getTimer(name);
        if ( t ){
            std::cout << "Timer [" << name << "] elapsed time: "
                      << t->elapsed << "s\n";
        }
    }

    static inline double getTime()
    {
        struct timeval time;
        gettimeofday( &time, 0 );
        return (double)time.tv_sec + (double)time.tv_usec / 1e6;
    }

    void coutTotal( std::string name )
    {
        single_timer * t = getTimer(name);
        if ( t ){
             std::cout << "Timer [" << name << "] total time: "
                      << t->total << "s\n";
        }
    }

    
    void start( std::string name )
    {

        single_timer * t = getTimer(name, false);
        if ( print ){
            std::cout << "Starting " << name << std::endl;
        }

        t->isStopped = false;
        t->tic = getTime();
    }

    double stop( std::string name )
    {
        
        single_timer * t = getTimer(name);
        if ( !t) {return 0.0;}

        if ( print ){
            std::cout << "Stopping " << name << std::endl;
        }
        
        t->toc = getTime();

        t->elapsed = t->toc - t->tic;

        if ( t->elapsed > 0 ){ t->total += t->elapsed; }

        t->count ++;
        t->isStopped = true;

        return t->elapsed;
    }

    double reset( std::string name)
    {
        single_timer * t = getTimer(name);
        if ( !t ) {return 0.0;}

        double temp = t->total;
        t->total = 0;
        t->count = 0;
        t->elapsed = 0;
        return temp;
    }

    double getTotal( std::string name )
    {
        single_timer * t = getTimer(name);
        if ( !t ) {return 0.0;}
        return t->total;
    }
    
    double getElapsed( std::string name )
    {
        single_timer * t = getTimer(name);
        if ( !t ) {return 0.0;}
        return t->elapsed;
    }

    //returns true if the timer is already in a start state, elsewise,
    //  it is started.
    bool tryStart( std::string name )
    {
        single_timer * t = getTimer(name, false);

        if ( t->isStopped ){
            t->isStopped = false;
            t->tic = getTime();
            return false;
        }

        return false;
    }
    
    //returns the elapsed time if the timer is already in a
    //  stopped state, elsewise, it is started.
    bool tryStop( std::string name )
    {
        single_timer * t = getTimer(name);

        if ( !t ) {return false;}
        if ( t->isStopped ){
            return true;
        }
        
        t->toc = getTime();
        t->elapsed = t->toc - t->tic;

        if ( t->elapsed > 0 ){ t->total += t->elapsed; }

        t->count ++;
        t->isStopped = true;

        return false;
    }

    unsigned int getCount( std::string name )
    {
        single_timer * t = getTimer(name);
        if ( !t ) {return 0;}

        return t->count;
    }

    inline void wait( double time ) { sleep( time ); }

    void getAllTotal(std::vector< std::pair<std::string, double> > & times)
        const
    {
        times.resize( timers.size() );
        int i = 0;

        for ( Map::const_iterator it = timers.begin();
              it != timers.end();
              it ++ , i ++)
        {
            times[i].first = it->first;
            times[i].second = it->second.total;
        }
    }
};

}//namespace

#endif
