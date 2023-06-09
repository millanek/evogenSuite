//
//  DistanceMatrix.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 21.09.22.
//

#include "UtilsGeneral.hpp"
#include "UtilsAnnotation.hpp"
#include "UtilsSetInfo.hpp"
#include "UtilsSetCounts.hpp"
#include "UtilsStats.hpp"

#ifndef DistanceMatrix_hpp
#define DistanceMatrix_hpp

int globalStatsMain(int argc, char** argv);
void parseStatsOptions(int argc, char** argv);

void normaliseMatrixBasedOnMissigness(std::vector<std::vector<double>>& dxyMatrix, const std::vector<std::vector<int>>& missignessMatrix, const int totalSites);
void scaleMatrixByAccessibleBP(std::vector<std::vector<double>>& dxyMatrix, const int accessibleBP);
bool bFirstTwoAllelesInformative(const SetInformation& setInfo, const GeneralSetCounts* c, const double maxMissing);

class DistanceMatrices {
public:
    
    DistanceMatrices(const SetInformation& setInfo, const string outFileRoot,const string runName, const bool bBootstrap, const bool bDoubletons) : bootstrapBlockNum(0)  {
        
        int n = (int)setInfo.populations.size();
        initialize_matrix_double(diffMatrix, n);
        initialize_matrix_int(sharedPolymorphisms, n);
        initialize_matrix_int(pairwiseMissingness, n);
        privateVarCounts.resize(n,0);
        
        string dxyFileName = outFileRoot + "_" + runName + "_distanceMatrix.txt";
        pDxyOutFile = new std::ofstream(dxyFileName.c_str());
        
        if (bBootstrap) {
            initialize_matrix_double(thisBootstrapBlock, n);
            initialize_matrix_int(thisBootstrapBlockMissingness, n);
            string bootFileName = outFileRoot + "_" + runName + "_distanceMatrix_boot.txt";
            pBootOutFile = new std::ofstream(bootFileName.c_str());
        }
        
        if (bDoubletons) {
            initialize_matrix_int(doubletonMatrix, n);
            string doubletonFileName = outFileRoot + "_" + runName + "_doubletonMatrix.txt";
            pDoubletonOutFile = new std::ofstream(doubletonFileName.c_str());
        }
    };
    
    void addAllPairwiseDistances(const SetInformation& setInfo, const GeneralSetCounts* c, const double maxMissing, const bool bDoBootstrap) {
        for (int i = 0; i < setInfo.populations.size(); i++) {
            // std::cerr << "setInfo.populations[i]: " <<  setInfo.populations[i] << std::endl;

            for (int j = 0; j <= i; j++) {
                string set1 = setInfo.populations[i];
                string set2 = setInfo.populations[j];
               // std::cerr << "setInfo.populations[j]: " <<  setInfo.populations[j] << std::endl;
                int n1 = c->setRefCounts.at(set1) + vector_sum(c->setAltAlleleCounts.at(set1));
                int n2 = c->setRefCounts.at(set2) + vector_sum(c->setAltAlleleCounts.at(set2));
               // std::cerr << "n1: " <<  n1 << std::endl;
               // std::cerr << "n2: " <<  n2 << std::endl;
                int set1FullSize = 2*(int)setInfo.popToPosMap.at(set1).size();
                int set2FullSize = 2*(int)setInfo.popToPosMap.at(set2).size();
                
                if ( (double)n1/set1FullSize <= (1 - maxMissing) || (double)n2/set2FullSize <= (1 - maxMissing) ) {
                    pairwiseMissingness[i][j]++;
                    if (bDoBootstrap) thisBootstrapBlockMissingness[i][j]++;
                } else {
                    double Dxy = DxyPerSNPfromSetAlleles(c,setInfo.populations[i],setInfo.populations[j]);
                    diffMatrix[i][j] += Dxy;
                    if(bDoBootstrap) thisBootstrapBlock[i][j] += Dxy;
                }
            }
        }
    }
    
    void addBootstrapBlock() {
        bootstrapBlocksDiff[bootstrapBlockNum] = thisBootstrapBlock;
        bootstrapBlocksMissingness[bootstrapBlockNum] = thisBootstrapBlockMissingness;
        reset_matrix_to_zero(thisBootstrapBlock);
        reset_matrix_to_zero(thisBootstrapBlockMissingness);
        bootstrapBlockNum++;
    }
    
    void doBootstrapResampling(const SetInformation& setInfo, const int nReplicates, const int totalSites, const int numAccessibleBP) {
        for (int n = 0; n < nReplicates; n++) {
            std::vector<std::vector<double> > bootReplicateMatrix;
            initialize_matrix_double(bootReplicateMatrix, (int)setInfo.populations.size());
            std::vector<std::vector<int> > bootReplicateMatrixMissingness;
            initialize_matrix_int(bootReplicateMatrixMissingness, (int)setInfo.populations.size());
            for (int b = 0; b < (int)bootstrapBlocksDiff.size(); b++) {
                int block = rand() % bootstrapBlocksDiff.size();
                for (int i = 0; i < (int)diffMatrix.size(); i++) {
                    for (int j = 0; j < (int)diffMatrix.size(); j++) {
                        bootReplicateMatrix[i][j] = bootReplicateMatrix[i][j] + bootstrapBlocksDiff[block][i][j];
                        bootReplicateMatrixMissingness[i][j] = bootReplicateMatrixMissingness[i][j] + bootstrapBlocksMissingness[block][i][j];
                    }
                }
            }
            
            // Normalise the bootstrap diff matrix based on pairwise missingness
            normaliseMatrixBasedOnMissigness(bootReplicateMatrix, bootReplicateMatrixMissingness, totalSites);
            if (numAccessibleBP > -1) scaleMatrixByAccessibleBP(bootReplicateMatrix, numAccessibleBP);
            
            print_vector(setInfo.populations,*pBootOutFile);
            print_matrix<const std::vector<std::vector<double> >&>(bootReplicateMatrix, *pBootOutFile);
            pBootOutFile->close();
        }
    }
    
    void printDxyResults(const SetInformation& setInfo) {
        print_vector(setInfo.populations,*pDxyOutFile);
        pDxyOutFile->precision(10);
        print_matrix<const std::vector<std::vector<double> >&>(diffMatrix, *pDxyOutFile);
    }
    
    void printDoubletonResults(const SetInformation& setInfo) {
        print_vector(setInfo.populations,*pDoubletonOutFile);
        pDoubletonOutFile->precision(10);
        print_matrix<const std::vector<std::vector<int> >&>(doubletonMatrix, *pDoubletonOutFile);
    }
    
    void doubleton_analysis(const SetInformation& setInfo, const GeneralSetCounts* c, double maxMissing) {
        
        
        
        if (c->getOverallCountOfFirstAltAllele() == 2) { // The first alternate allele is a doubleton
            std::vector<int> doubletonPops;
            for (int i = 0; i < (int)setInfo.populations.size(); i++) {
                int countInThisPopulation = c->setAltAlleleCounts.at(setInfo.populations[i])[0];
                if (countInThisPopulation == 2) { doubletonPops.push_back(i); doubletonPops.push_back(i); break; }
                if (countInThisPopulation == 1) { doubletonPops.push_back(i); }
            }
            assert(doubletonPops.size() == 2);
            doubletonMatrix[doubletonPops[0]][doubletonPops[1]]++;
        } else if (c->getOverallCountOfRefAllele() == 2) { // The reference allele is a doubleton
            std::vector<int> doubletonPops;
            for (int i = 0; i < (int)setInfo.populations.size(); i++) {
                int countInThisPopulation = c->setRefCounts.at(setInfo.populations[i]);
                if (countInThisPopulation == 2) { doubletonPops.push_back(i); doubletonPops.push_back(i); break; }
                if (countInThisPopulation == 1) { doubletonPops.push_back(i); }
            }
            assert(doubletonPops.size() == 2);
            doubletonMatrix[doubletonPops[0]][doubletonPops[1]]++;
        }
    }
    
    void privateVars_analysis(const SetInformation& setInfo, const GeneralSetCounts* c, double maxMissing) {
        
        std::vector<int> populationsWithRefPositiveCount;
        // Check if only one population has the ref allele
        for (int i = 0; i < (int)setInfo.populations.size(); i++) {
            int countInThisPopulation = c->setRefCounts.at(setInfo.populations[i]);
            if (countInThisPopulation > 0) populationsWithRefPositiveCount.push_back(i);
        }
        if (populationsWithRefPositiveCount.size() == 1)
            privateVarCounts[populationsWithRefPositiveCount[0]]++;
        
        std::vector<int> populationsWithFirstAltPositiveCount;
        // Check if only one population has the first alt allele
        for (int i = 0; i < (int)setInfo.populations.size(); i++) {
            int countInThisPopulation = c->setAltAlleleCounts.at(setInfo.populations[i])[0];
            if (countInThisPopulation > 0) populationsWithFirstAltPositiveCount.push_back(i);
        }
        if (populationsWithFirstAltPositiveCount.size() == 1)
            privateVarCounts[populationsWithFirstAltPositiveCount[0]]++;
        
        
    }
    
    void sharedPolymorphism_analysis(const SetInformation& setInfo, const GeneralSetCounts* c, double maxMissing) {
        
        for (int i = 0; i < (int)setInfo.populations.size(); i++) {
            int refCountInThisPopulation = c->setRefCounts.at(setInfo.populations[i]);
            int altCountInThisPopulation = c->setAltAlleleCounts.at(setInfo.populations[i])[0];
            if (refCountInThisPopulation > 0 && altCountInThisPopulation > 0) {
                for (int j = 0; j < (int)setInfo.populations.size(); j++) {
                    int refCountInThisPopulation = c->setRefCounts.at(setInfo.populations[i]);
                    int altCountInThisPopulation = c->setAltAlleleCounts.at(setInfo.populations[i])[0];
                    if (refCountInThisPopulation > 0 && altCountInThisPopulation > 0) {
                        sharedPolymorphisms[i][j]++;
                    }
                    
                }
            }
        }
    }

    
    std::vector<std::vector<double>> diffMatrix;
    std::vector<std::vector<int>> doubletonMatrix;
    std::vector<std::vector<int>> sharedPolymorphisms;
    
    std::vector<std::vector<int>> pairwiseMissingness;
    
    std::vector<std::vector<double>> thisBootstrapBlock;
    std::vector<std::vector<int>> thisBootstrapBlockMissingness;
    
    std::unordered_map<int, std::vector<std::vector<double>>>  bootstrapBlocksDiff;
    std::unordered_map<int , std::vector<std::vector<int>>>  bootstrapBlocksMissingness;
    
    std::vector<int> privateVarCounts;
    
    int bootstrapBlockNum;
    
    std::ofstream* pBootOutFile;
    std::ofstream* pDxyOutFile;
    std::ofstream* pDoubletonOutFile;
    std::ofstream* pSharedPolymorphismOutFile;
    
};



#endif /* DistanceMatrix_hpp */
