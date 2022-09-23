//
//  PBS.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef PBS_hpp
#define PBS_hpp

#include "UtilsGeneral.hpp"
#include "UtilsStats.hpp"

void parsePBSoptions(int argc, char** argv);
int PBSmain(int argc, char** argv);

class PBStrios {
public:
    
    PBStrios(const string& PBStriosFileName, const string& runName, const int windowSize, const int windowStep, const int fixedWindowSize, const bool bAnnotationPresent) {
        std::ifstream* PBStriosFile = new std::ifstream(PBStriosFileName.c_str());
        if (!PBStriosFile) {
            std::cerr << "ERROR: The file " << PBStriosFile << " could not be opened\n";
            exit(1);
        }
        string line;
        while (getline(*PBStriosFile,line)) {
            // std::cerr << line << std::endl;
            std::vector<string> threePops = split(line, '\t'); assert(threePops.size() == 3);
            std::ofstream* outFile = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBS_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
            std::ofstream* outFileFixedWindow = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBS_" + runName + "_FW" + numToString(fixedWindowSize) + ".txt");
            *outFile << "chr\twStart\twEnd\t" << threePops[0] << "\t" << threePops[1] << "\t" << threePops[2] << std::endl;
            *outFileFixedWindow << "chr\twStart\twEnd\t" << threePops[0] << "\t" << threePops[1] << "\t" << threePops[2] << "\t" << "nFwSNPs1" << "\t" << "nFwSNPs2" << "\t" << "nFwSNPs3" << std::endl;
            //outFile->setf(std::ios_base::fixed); // Avoid scientific notation in the coordinates
            outFiles.push_back(outFile); outFilesFixedWindow.push_back(outFileFixedWindow);
            if (bAnnotationPresent) {
                std::ofstream* outFileGenes = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBSGenes_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
                *outFileGenes << "gene\t" << "numSNPsExons\tnumSNPsIntrons\tnumSNPs3kbPromoter\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons\t" << threePops[0] << "_promoter\t" << threePops[1] << "_promoter\t" << threePops[2] << "_promoter" << std::endl;
                outFilesGenes.push_back(outFileGenes);
            }
            trios.push_back(threePops);
        }
        
        // Need to prepare the vectors to hold the PBS values and the coordinates:
        initResultsSNPwindows(windowSize);
        initResultsPhysicalWindows();
        if (bAnnotationPresent) initResultsGenes();
        
        usedVars.resize(trios.size(),0);
    };
    
    std::vector<std::vector<string> > trios;
    
    std::vector<std::ofstream*> outFiles;
    std::vector<std::ofstream*> outFilesFixedWindow;
    std::vector<std::ofstream*> outFilesGenes;
    
    std::vector<int> usedVars;          // Will count the number of used variants for each trio
    
    std::vector<std::vector<std::deque<double>>> resultsSNPwindows;
    std::vector<std::vector<std::vector<double>>> resultsPhysicalWindows;
    std::vector<std::vector<std::vector<double>>> resultsGenes;
    
    void summariseAndOutputPerGene(string geneName);
    
    private:
        void initResultsSNPwindows(const int windowSize) {
            std::deque<double> initDeq(windowSize,0.0); // deque to initialise per-site PBS values
            std::vector<std::deque<double>> initThreeDeques(4,initDeq); // vector of three per-site PBS deques - for each population in the trio, and the fourth is for the coordinates
            std::vector<std::vector<std::deque<double>>> PBSresults(trios.size(),initThreeDeques);
            resultsSNPwindows = PBSresults;
        }
        
        void initResultsPhysicalWindows() {
            std::vector<std::vector<double>> initVectorFixed(3);
            std::vector<std::vector<std::vector<double>>> PBSfixedWindowResults(trios.size(),initVectorFixed);
            resultsPhysicalWindows = PBSfixedWindowResults;
        }
        
        void initResultsGenes() {
            std::vector<std::vector<double>> initGeneVectors(9); // For the nine PBS columns in the _PBSGenes_ files
            std::vector<std::vector<std::vector<double>>> PBSgeneResults(trios.size(),initGeneVectors);
            resultsGenes = PBSgeneResults;
        }
};





#endif /* PBS_hpp */
