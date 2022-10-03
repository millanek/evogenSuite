//
//  UtilsStats.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "UtilsStats.hpp" 

// -------------------------------------------------------------------------
// 1. Functions to calculate PBS:

std::vector<double> calculatePBSfromAFs(const double p1, const double p2, const double p3, const double p1AlleleCount, const double p2AlleleCount, const double p3AlleleCount) {
    // std::cerr << "p:\t" << p1 << "\t" << p2 << "\t" << p3 <<  std::endl;
    double Fst12; double Fst13; double Fst23;
    double power12 = pow(p1-p2, 2);
    double power13 = pow(p1-p3, 2);
    double power23 = pow(p2-p3, 2);
    double fraction1 = (p1*(1-p1))/(p1AlleleCount-1);
    double fraction2 = (p2*(1-p2))/(p2AlleleCount-1);
    double fraction3 = (p3*(1-p3))/(p3AlleleCount-1);
    double numerator12 = power12 - fraction1 - fraction2;
    double numerator13 = power13 - fraction1 - fraction3;
    double numerator23 = power23 - fraction2 - fraction3;
    double denominator12 = (p1*(1-p2))+(p2*(1-p1));
    double denominator13 = (p1*(1-p3))+(p3*(1-p1));
    double denominator23 = (p2*(1-p3))+(p3*(1-p2));
    if ((p1 == 0 && p2 == 0) || (p1 == 1 && p2 == 1)) { Fst12 = 0.0; } else { Fst12 = numerator12/denominator12; }
    if ((p1 == 0 && p3 == 0) || (p1 == 1 && p3 == 1)) { Fst13 = 0.0; } else { Fst13 = numerator13/denominator13; }
    if ((p2 == 0 && p3 == 0) || (p2 == 1 && p3 == 1)) { Fst23 = 0.0; } else { Fst23 = numerator23/denominator23; }
    // std::cerr << Fst12 << "\t" << Fst13 << "\t" << Fst23 <<  std::endl;
    if (Fst12 < 0) Fst12 = 0; if (Fst13 < 0) Fst13 = 0; if (Fst23 < 0) Fst23 = 0;
    if (Fst12 == 1) Fst12 = 1 - (Fst12/p1AlleleCount); if (Fst13 == 1) Fst13 = 1 - (Fst13/p1AlleleCount); if (Fst23 == 1) Fst23 = 1 - (Fst23/p2AlleleCount);
    double T12 = -log(1-Fst12); double T13 = -log(1-Fst13); double T23 = -log(1-Fst23);
    double PBS1 = (T12+T13-T23)/2.0;
    double PBS2 = (T12+T23-T13)/2.0;
    double PBS3 = (T13+T23-T12)/2.0;
    if (PBS1 < 0) PBS1 = 0; if (PBS2 < 0) PBS2 = 0; if (PBS3 < 0) PBS3 = 0;
    //std::cerr << PBS1 << "\t" << PBS2 << "\t" << PBS3 <<  std::endl;
    std::vector<double> PBS; PBS.push_back(PBS1); PBS.push_back(PBS2); PBS.push_back(PBS3);
    return PBS;
}



// -------------------------------------------------------------------------
// 2. Functions to calculate Dxy:

double DxyPerSNPfromAFs(double AF1, double AF2) {
    double dxy = AF1*(1-AF2) + AF2*(1-AF1);
    return dxy;
}


double DxyPerSNPfromSetAlleles(const GeneralSetCounts* c, const string& set1, const string& set2) {
    double Dxy = 0;
    const std::vector<int> allelesSet1 = c->setAlleles.at(set1);
    const std::vector<int> allelesSet2 = c->setAlleles.at(set2);
    const int n1 = (int)allelesSet1.size(); const int n2 = (int)allelesSet1.size();
    
    int sumKij = 0;
    for (std::vector<std::string>::size_type i = 0; i != n1; i++) {
        for (std::vector<std::string>::size_type j = 0; j != n2; j++) {
            if (allelesSet1[i] != allelesSet2[j]) sumKij = sumKij + 1;
        }
        Dxy = (double)sumKij/(n1*n2);
    }
    return Dxy;
}


// -------------------------------------------------------------------------
// 3. Functions to calculate Inbreeding coefficient (per SNP):

double calculateInbreedingCoefficient(std::vector<int>& individualsWithVariant) {
    int naa = 0; int nAa = 0; int nAA = 0;
    for (std::vector<int>::size_type i = 0; i != individualsWithVariant.size(); i++) {
        if (individualsWithVariant[i] == 0) naa++;
        if (individualsWithVariant[i] == 1) nAa++;
        if (individualsWithVariant[i] == 2) nAA++;
    }
    
    // Get the proportions of alt-hom and hets
    double pAA = (double)nAA / individualsWithVariant.size();
    double pAa = (double)nAa / individualsWithVariant.size();
    
    // Allele frequencies
    double p = pAA + (0.5 * pAa);
    double q = 1 - p;
    
    // Get the Hardy-Weinberg prediction for expected number of heterozygotes:
    double HWAa = 2*p*q;

    
    // Get the inbreeding coefficient
    double F = (HWAa - pAa) / HWAa;
    return F;
}

/*double calculateChiSqPvalForInbreeding(std::vector<int>& individualsWithVariant) {
    // Get the observed numbers of genotypes
    int naa = 0; int nAa = 0; int nAA = 0;
    for (std::vector<int>::size_type i = 0; i != individualsWithVariant.size(); i++) {
        if (individualsWithVariant[i] == 0) naa++;
        if (individualsWithVariant[i] == 1) nAa++;
        if (individualsWithVariant[i] == 2) nAA++;
    }
    
    // Get the observed proportions of alt-hom and hets
    double pAA = (double)nAA / individualsWithVariant.size();
    double pAa = (double)nAa / individualsWithVariant.size();
    // Allele frequencies
    double p = pAA + (0.5 * pAa);
    double q = 1 - p;
    // Get the Hardy-Weinberg prediction for expected proportion of heterozygotes:
    double HWAa = 2*p*q;
    
    // Get the expected numbers of genotypes:
    double expAA = pow(p, 2) * individualsWithVariant.size();
    double expAa = HWAa * individualsWithVariant.size();
    double expaa = pow(q, 2) * individualsWithVariant.size();
    
    // Perform the chi-Sqared test
    std::vector<double> observed{static_cast<double>(nAA),static_cast<double>(nAa),static_cast<double>(naa)};
    std::vector<double> expected{expAA,expAa,expaa};
    double pVal = pearson_chi_sq_goodness_of_fit(observed, expected, 1);
    return pVal;
} */
