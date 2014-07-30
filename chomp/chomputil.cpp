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




}//namespace
