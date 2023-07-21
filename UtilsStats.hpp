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


inline double calculateFstNumerator(const double p1, const double p2, const int n1, const int n2) {
    double power = pow((p1-p2), 2);
    double fraction1 = (p1*(1-p1))/(n1-1);
    double fraction2 = (p2*(1-p2))/(n2-1);
    double numerator = power - fraction1 - fraction2;
    return numerator;
}

inline double calculateFstDenominator(const double p1, const double p2) {
    double denominator = (p1*(1-p2))+(p2*(1-p1));
    return denominator;
}


template <class T> double calculateFst(const T& fstNumerators, const T& fstDenominators, const bool bZeroRounding) {
    double numeratorAverage = vector_average(fstNumerators);
    double denominatorAverage = vector_average(fstDenominators);
    double Fst = numeratorAverage/denominatorAverage;
    if (bZeroRounding && Fst < 0) Fst = 0;
    return Fst;
}


double calculateInbreedingCoefficient(std::vector<int>& individualsWithVariant);
//double calculateChiSqPvalForInbreeding(std::vector<int>& individualsWithVariant);

// -------------------------------------    BASIC MATH/STATS  ----------------------------------------

// factorial(x): (x! for non-negative integer x) is defined to be gamma(x+1) (as in R)
inline double factorial(double num) {
    if (num < 0) {
        std::cerr << "Can't compute factorial of a negative number " << num << std::endl;
        exit(1);
    }
    return tgamma(num+1);
}

// Calculates the binomial coefficient (n choose k)
inline int choose(int n, int k) {
    double dResult = factorial(n)/(factorial(k)*factorial(n-k));
    int iResult = (int)round(dResult);
    return iResult;
}

// standard deviation of a sample
template <class T> double std_dev(T vector) {
    double mean = vector_average(vector);
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += pow((vector[i] - mean), 2.0);
    }
    double std_dev = sqrt((double)sum/(double)(vector.size()-1));
    return std_dev;
}


inline void copy_except(int i, std::vector<double>& inVec, std::vector<double>& outVec) {
    std::copy(inVec.begin(), inVec.begin() + i, outVec.begin());
    std::copy(inVec.begin() + i + 1, inVec.end(), outVec.begin()+i);
    //std::cerr << "copying:" << i << " "; print_vector_stream(inVec, std::cerr);
    //std::cerr << "copied: " << i << " "; print_vector_stream(outVec, std::cerr);
}

// jackknive standard error
template <class T> double jackknive_std_err(T& vector) {
    std::vector<double> jackkniveAverages;
    std::vector<double> JregionDs; JregionDs.resize(vector.size()-1);
    for (std::vector<double>::size_type i = 0; i != vector.size(); i++) {
        // std::cerr << "copying " << i << std::endl;
        copy_except(i, vector, JregionDs);
        jackkniveAverages.push_back(vector_average(JregionDs));
        JregionDs.clear(); JregionDs.resize(vector.size()-1);
    }
    double jackkniveOverallMean = vector_average(jackkniveAverages);
    double sum = 0;
    for (int i = 0; i < jackkniveAverages.size(); i++) {
        sum += pow((jackkniveAverages[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveAverages.size()-1)/(double)jackkniveAverages.size()) * sum;
    double Dstd_err = sqrt(var);
    return Dstd_err;
}



#endif /* UtilsStats_hpp */
