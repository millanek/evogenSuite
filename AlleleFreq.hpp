//
//  AlleleFreq.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef AlleleFreq_hpp
#define AlleleFreq_hpp

#include <stdio.h>
#include "UtilsGeneral.hpp" 

void parseAFoptions(int argc, char** argv);
int AFmain(int argc, char** argv);

void printAFheader(const SetInformation& setInfo, std::ostream* outFileAF);

#endif /* AlleleFreq_hpp */
