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


std::vector<double> calculatePBSnumerators(const GeneralSetCounts* c, const std::vector<string>& PBStrio) {
    double power12 = pow((c->setAAFs.at(PBStrio[0])-c->setAAFs.at(PBStrio[1])), 2);
    double power13 = pow((c->setAAFs.at(PBStrio[0])-c->setAAFs.at(PBStrio[2])), 2);
    double power23 = pow((c->setAAFs.at(PBStrio[1])-c->setAAFs.at(PBStrio[2])), 2);
    double fraction1 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[0])))/(c->setAlleleCounts.at(PBStrio[0])-1);
    double fraction2 = (c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[1])))/(c->setAlleleCounts.at(PBStrio[1])-1);
    double fraction3 = (c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[2])))/(c->setAlleleCounts.at(PBStrio[2])-1);
    double numerator12 = power12 - fraction1 - fraction2;
    double numerator13 = power13 - fraction1 - fraction3;
    double numerator23 = power23 - fraction2 - fraction3;
    std::vector<double> PBSnumerators; PBSnumerators.push_back(numerator12);
    PBSnumerators.push_back(numerator13); PBSnumerators.push_back(numerator23);
    return PBSnumerators;
}

std::vector<double> calculatePBSdenominator(const GeneralSetCounts* c, const std::vector<string>& PBStrio) {
    double denominator12 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[1])))+(c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[0])));
    double denominator13 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[2])))+(c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[0])));
    double denominator23 = (c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[2])))+(c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[1])));
    std::vector<double> PBSdenominators; PBSdenominators.push_back(denominator12);
    PBSdenominators.push_back(denominator13); PBSdenominators.push_back(denominator23);
    return PBSdenominators;
}


// -------------------------------------------------------------------------
// 2. Functions to calculate XXX:

double DxyPerSNPfromAFs(double AF1, double AF2) {
    double dxy = AF1*(1-AF2) + AF2*(1-AF1);
    return dxy;
}


double calculateDxy(const SetCounts& thisVarCounts, const int n1, const int n2) {
    double Dxy;
    int sumKij = 0;
    for (std::vector<std::string>::size_type i = 0; i != n1/2; i++) {
        for (std::vector<std::string>::size_type j = 0; j != n2/2; j++) {
            if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 0) {
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 0) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 2) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 2) {
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 0) {
                sumKij = sumKij + 4;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 2) {
                sumKij = sumKij + 4;
            }
        }
    }
    Dxy = (double)sumKij/(n1*n2);
    return Dxy;
}
