//
//  UtilsPopulationComparisons.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 06.10.22.
//

#ifndef UtilsPopulationComparisons_hpp
#define UtilsPopulationComparisons_hpp

#include "UtilsGeneral.hpp"

class PopulationComparisons {
public:
    PopulationComparisons(const string& ComparisonsFileName) {
        
        ComparisonsFile = new std::ifstream(ComparisonsFileName.c_str());
        assertFileOpen(*ComparisonsFile, ComparisonsFileName);
    };
    
    std::vector<std::vector<string> > comparisons;
    
    std::ifstream* ComparisonsFile;
    
    std::vector<std::ofstream*> outFiles;
    std::vector<std::ofstream*> outFilesFixedWindow;
    std::vector<std::ofstream*> outFilesGenes;
    std::vector<std::ofstream*> outFilesRegions;
    
    std::vector<int> usedVars;          // Will count the number of used variants for each trio
    
    std::vector<std::vector<std::deque<double>>> resultsSNPwindows;
    std::vector<std::vector<std::vector<double>>> resultsPhysicalWindows;
    std::vector<std::vector<std::vector<double>>> resultsGenes;
    
    protected:
        void initResultsSNPwindows(const int windowSize, const int numSNPwindowDeques) {
            std::deque<double> initDeq(windowSize,0.0); // deque to initialise per-site PBS values
            std::vector<std::deque<double>> initThreeDeques(numSNPwindowDeques,initDeq); // vector of three per-site PBS deques - for each population in the trio, and the fourth is for the coordinates
            std::vector<std::vector<std::deque<double>>> PBSresults(comparisons.size(),initThreeDeques);
            resultsSNPwindows = PBSresults;
        }
        
        void initResultsPhysicalWindows(const int numPhysicalWindowVectors) {
            std::vector<std::vector<double>> initVectorFixed(numPhysicalWindowVectors);
            std::vector<std::vector<std::vector<double>>> PBSfixedWindowResults(comparisons.size(),initVectorFixed);
            resultsPhysicalWindows = PBSfixedWindowResults;
        }
        
        void initResultsGenes(const int numGeneVectors) {
            std::vector<std::vector<double>> initGeneVectors(numGeneVectors); // For the nine PBS columns in the _PBSGenes_ files
            std::vector<std::vector<std::vector<double>>> PBSgeneResults(comparisons.size(),initGeneVectors);
            resultsGenes = PBSgeneResults;
        }
};


#endif /* UtilsPopulationComparisons_hpp */
