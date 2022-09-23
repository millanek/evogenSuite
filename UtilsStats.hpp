//
//  UtilsStats.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef UtilsStats_hpp
#define UtilsStats_hpp

#include "UtilsGeneral.hpp"

std::vector<double> calculatePBSfromAFs(const double p1, const double p2, const double p3, const double p1AlleleCount, const double p2AlleleCount, const double p3AlleleCount);
std::vector<double> calculatePBSnumerators(const GeneralSetCounts* c,const std::vector<string>& PBStrio);
std::vector<double> calculatePBSdenominator(const GeneralSetCounts* c,const std::vector<string>& PBStrio);


double DxyPerSNPfromAFs(double AF1, double AF2);
double calculateDxy(const SetCounts& thisVarCounts, const int n1, const int n2);

template <class T> double calculateFst(const T& fstNumerators, const T& fstDenominators) {
    double numeratorAverage = vector_average(fstNumerators);
    double denominatorAverage = vector_average(fstDenominators);
    double Fst = numeratorAverage/denominatorAverage;
    if (Fst < 0) Fst = 0;
    return Fst;
}


#endif /* UtilsStats_hpp */
