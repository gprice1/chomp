#ifndef _MOTION_OPTIMIZER_EXCEPTION_H_
#define _MOTION_OPTIMIZER_EXCEPTION_H_

#include <iostream>
#include <exception>

namespace mopt { 

class MotionOptimizerException: public std::exception
{
    virtual const char* what() const throw()
    {
        return "General motion optimization failure.
               " It may or may not be your fault.";
    }
};

class BadParameterException : public MotionOptimizerException
{
    virtual const char* what() const throw()
    {
        return "BadParameterException thrown."
               " You passed a bad parameter to a class.";
    }
};

class ObserverException : public MotionOptimizerException 
{

  private: 
    
  public: 
    
    virtual const char* what() const throw()
    {
        return "ObserverException thrown. The observer threw an exception.";
    }
};

}// namespace

#endif
