//
//  UtilsStats.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef UtilsStats_hpp
#define UtilsStats_hpp

#include "UtilsGeneral.hpp"
#include "UtilsSetCounts.hpp"

std::vector<double> calculatePBSfromAFs(const double p1, const double p2, const double p3, const double p1AlleleCount, const double p2AlleleCount, const double p3AlleleCount);


double DxyPerSNPfromAFs(double AF1, double AF2);
double DxyPerSNPfromSetAlleles(const GeneralSetCounts* c, const string& set1, const string& set2); 

template <class T> double calculateFst(const T& fstNumerators, const T& fstDenominators, const bool bZeroRounding) {
    double numeratorAverage = vector_average(fstNumerators);
    double denominatorAverage = vector_average(fstDenominators);
    double Fst = numeratorAverage/denominatorAverage;
    if (bZeroRounding && Fst < 0) Fst = 0;
    return Fst;
}


double calculateInbreedingCoefficient(std::vector<int>& individualsWithVariant);
//double calculateChiSqPvalForInbreeding(std::vector<int>& individualsWithVariant);

#endif /* UtilsStats_hpp */
