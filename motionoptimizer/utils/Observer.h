#ifndef _OBSERVER_H_
#define _OBSERVER_H_

#include "utils.h"

namespace mopt {

//Forward declaration of optimizer base
class OptimizerBase;

class Observer {
  public:
    virtual ~Observer(){}
    virtual int notify(const OptimizerBase& opt, 
                       EventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
        return 0;
    }
};

class DebugObserver: public Observer {
  public:
    virtual ~DebugObserver(){}
    virtual int notify(const OptimizerBase& c, 
                       EventType e,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
            
         std::string event_string;

         switch (e) {
            case INIT:
                 event_string = "INIT"; break;
            case CHOMP_ITER:
                 event_string = "CHOMP_ITER"; break;
            case CHOMP_LOCAL_ITER:
                 event_string = "CHOMP_LOCAL_ITER"; break;
            case NLOPT_ITER:
                 event_string = "NLOPT_ITER"; break;
            case FINISH:
                 event_string = "FINISH"; break;
            case TIMEOUT:
                 event_string = "TIMEOUT"; break;
            default:
                 event_string = "[INVALID]"; break;
         }

         std::cout << "debug: "
              << "event=" << event_string << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";
        if ( e != INIT && e != FINISH && 
            (std::isnan(curObjective) || std::isinf(curObjective) ||
             std::isnan(lastObjective) || std::isinf(lastObjective)) )
        {
            return 1;
        }
        return 0;
    }
};

} // namespace
#endif
