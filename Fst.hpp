//
//  Fst.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef Fst_hpp
#define Fst_hpp

#include "UtilsGeneral.hpp"
#include "UtilsAnnotation.hpp"
#include "UtilsStats.hpp"
#include <math.h>

void parseFstOptions(int argc, char** argv);
int fstMain(int argc, char** argv);
bool bPairInformativeThisSNP(const double p1, const double p2);


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


inline double calculateExpectedHeterozygositySimple(const double p) {
    double q = 1 - p;
    double heterozygosity = 1 - (pow(p,2)+pow(q,2));
    return heterozygosity;
}


inline double calculateExpectedHeterozygosityNei78(const double p, const int n) {
    double q = 1 - p;
    double simpleHeterozygosity = 1 - (pow(p,2)+pow(q,2));
    double heterozygosity = (n*simpleHeterozygosity)/(n-1);
    return heterozygosity;
}



// Define column numbers for "scanResults"
#define iFstNum 0
#define iFstDenom 1
#define iDxy 2
#define iPi1 3
#define iPi2 4
#define iCoord 5

class FstPairs {
public:
    
    FstPairs(const string& FstPairsFileName, const string& runName, const int windowSize, const int windowStep, const int fixedWindowSize, const bool bAnnotationPresent) {
        std::ifstream* FstPairsFile = new std::ifstream(FstPairsFileName.c_str());
        if (FstPairsFile->fail()) {
            std::cerr << "ERROR: The file " << FstPairsFileName << " could not be opened\n";
            exit(1);
        }
        std::cout << "\nCalculating statistics for the set pair(s):\n";
        string line;
        while (getline(*FstPairsFile,line)) {
            // std::cerr << line << std::endl;
            std::vector<string> twoPops = split(line, '\t'); assert(twoPops.size() == 2);
            print_vector(twoPops, std::cout);
            std::ofstream* outFile = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_Fst_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
            std::ofstream* outFileFixedWindow = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_" +  "_Fst_" + runName + "_FW" + numToString(fixedWindowSize) + ".txt");
            
            string fstOufileHeader = "chr\twStart\twEnd\tFst\tDxy\t" + twoPops[0] + "_pi\t" + twoPops[1] + "_pi\tAccessible_bp\n";
            *outFile << fstOufileHeader; *outFileFixedWindow << fstOufileHeader;
            outFile->setf(std::ios_base::fixed); // Avoid scientific notation in the coordinates
            outFiles.push_back(outFile); outFilesFixedWindow.push_back(outFileFixedWindow);
      
            /*     if (bAnnotationPresent) {
                std::ofstream* outFileGenes = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBSGenes_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
                *outFileGenes << "gene\t" << "numSNPsExons\tnumSNPsIntrons\tnumSNPs3kbPromoter\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons\t" << threePops[0] << "_promoter\t" << threePops[1] << "_promoter\t" << threePops[2] << "_promoter" << std::endl;
                outFilesGenes.push_back(outFileGenes);
            }*/
            pairs.push_back(twoPops);
        }
        
        // Need to prepare the vectors to hold the PBS values and the coordinates:
        initResultsSNPwindows(windowSize);
        initResultsPhysicalWindows();
        if (bAnnotationPresent) initResultsGenes();
        
        usedVars.resize(pairs.size(),0);
    };
    
    std::vector<std::vector<string> > pairs;
    
    std::vector<std::ofstream*> outFiles;
    std::vector<std::ofstream*> outFilesFixedWindow;
    std::vector<std::ofstream*> outFilesGenes;
    
    std::vector<int> usedVars;          // Will count the number of used variants for each trio
    
    std::vector<std::vector<std::deque<double>>> resultsSNPwindows;
    std::vector<std::vector<std::vector<double>>> resultsPhysicalWindows;
    std::vector<std::vector<std::vector<double>>> resultsGenes;
    
    void summariseAndOutputPerGene(string geneName);
    
    void addSNPresultsToWindows(const int pairNumber, const double thisSNPFstNumerator, const double thisSNPFstDenominator, const double thisSNPDxy, const double thisSNPpi1, const double thisSNPpi2, const double SNPcoordinate) {
        
        resultsSNPwindows[pairNumber][0].push_back(thisSNPFstNumerator); resultsSNPwindows[pairNumber][0].pop_front();
        resultsSNPwindows[pairNumber][1].push_back(thisSNPFstDenominator); resultsSNPwindows[pairNumber][1].pop_front();
        resultsSNPwindows[pairNumber][2].push_back(thisSNPDxy); resultsSNPwindows[pairNumber][2].pop_front();
        resultsSNPwindows[pairNumber][3].push_back(thisSNPpi1); resultsSNPwindows[pairNumber][3].pop_front();
        resultsSNPwindows[pairNumber][4].push_back(thisSNPpi2); resultsSNPwindows[pairNumber][4].pop_front();
        resultsSNPwindows[pairNumber][5].push_back(SNPcoordinate); resultsSNPwindows[pairNumber][5].pop_front();
        
        resultsPhysicalWindows[pairNumber][0].push_back(thisSNPFstNumerator);
        resultsPhysicalWindows[pairNumber][1].push_back(thisSNPFstDenominator);
        resultsPhysicalWindows[pairNumber][2].push_back(thisSNPDxy);
        resultsPhysicalWindows[pairNumber][3].push_back(thisSNPpi1);
        resultsPhysicalWindows[pairNumber][4].push_back(thisSNPpi2);
    }
    
    void finalizeAndOutputSNPwindow(const int pairNumber, const string chr, const double currentSNPcoordinate, const AccessibleGenome* ag) {
        // The starting coordinate of this window
        int startCoord = (int)resultsSNPwindows[pairNumber][iCoord][0];
        if (startCoord < currentSNPcoordinate) {
            int windowLength = ((int)currentSNPcoordinate - startCoord) + 1; // The length of this physical window
            int accessibleInThisWindow;
            if (ag->initialised) {
                accessibleInThisWindow = ag->getAccessibleBPinRegion(chr, startCoord, currentSNPcoordinate);
            } else {
                accessibleInThisWindow = windowLength;
            }
            
            double Fst = calculateFst(resultsSNPwindows[pairNumber][iFstNum], resultsSNPwindows[pairNumber][iFstDenom]);
            double Dxy = vector_average_withRegion(resultsSNPwindows[pairNumber][iDxy], accessibleInThisWindow);
            double Pi1 = vector_average_withRegion(resultsSNPwindows[pairNumber][iPi1], accessibleInThisWindow);
            double Pi2 = vector_average_withRegion(resultsSNPwindows[pairNumber][iPi2], accessibleInThisWindow);
            
            *outFiles[pairNumber] << chr << "\t" << (int)startCoord << "\t" << currentSNPcoordinate << "\t" << Fst << "\t" << Dxy << "\t" << Pi1 << "\t" << Pi2 << "\t" << accessibleInThisWindow << std::endl;
        }
    }
    
    
    void finalizeAndOutputPhysicalWindow(const int pairNumber, const int physicalWindowSize, const string chr, const double currentSNPcoord, const AccessibleGenome* ag, int& thisWindowStart, int& thisWindowEnd) {
        
        if (thisWindowStart < currentSNPcoord) { // As long as the current SNP coord is more than the previous (i.e., we are on the same chromosome; I should check this properly
            int accessibleInThisWindow;
            if (ag->initialised) {
                accessibleInThisWindow = ag->getAccessibleBPinRegion(chr, thisWindowStart, thisWindowEnd);
            } else {
                accessibleInThisWindow = physicalWindowSize;
            }
            
            
            double thisFixedWindowFst = 0; double thisFixedWindowDxy = 0;
            double thisFixedWindowPi1 = 0; double thisFixedWindowPi2 = 0;
            int nFwSNPs = (int)resultsPhysicalWindows[pairNumber][iFstNum].size();
            if (nFwSNPs > 0) {
                thisFixedWindowFst = calculateFst(resultsPhysicalWindows[pairNumber][iFstNum], resultsPhysicalWindows[pairNumber][iFstDenom]);
                thisFixedWindowDxy = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iDxy], accessibleInThisWindow);
                thisFixedWindowPi1 = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iPi1], accessibleInThisWindow);
                thisFixedWindowPi2 = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iPi2], accessibleInThisWindow);
                *outFilesFixedWindow[pairNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\t" << thisFixedWindowFst << "\t" << thisFixedWindowDxy << "\t" << thisFixedWindowPi1 << "\t" << thisFixedWindowPi2 << "\t" << accessibleInThisWindow << std::endl;
            } else { // If this physical window did not have any SNPs, then Fst="NA"
                *outFilesFixedWindow[pairNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\tNA\t" << thisFixedWindowDxy << "\t" << thisFixedWindowPi1 << "\t" << thisFixedWindowPi2 << "\t" << accessibleInThisWindow << std::endl;
            }
            
            for (int i = 0; i < resultsPhysicalWindows[pairNumber].size(); i++) {
                resultsPhysicalWindows[pairNumber][i].clear();
            }
            
            thisWindowStart = thisWindowStart + physicalWindowSize;
            thisWindowEnd = thisWindowEnd + physicalWindowSize;
            
            while (!(thisWindowStart <= currentSNPcoord && thisWindowEnd > currentSNPcoord)) {
                if (ag->initialised) {
                    accessibleInThisWindow = ag->getAccessibleBPinRegion(chr, thisWindowStart, thisWindowEnd);
                } else {
                    accessibleInThisWindow = physicalWindowSize;
                }
                *outFilesFixedWindow[pairNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << accessibleInThisWindow << std::endl;
                thisWindowStart = thisWindowStart + physicalWindowSize;
                thisWindowEnd = thisWindowEnd + physicalWindowSize;
            }
        } else {
            thisWindowStart = 0;
            thisWindowEnd = 0 + physicalWindowSize;
        }
    }
    
    
    
    private:
        void initResultsSNPwindows(const int windowSize) {
            std::deque<double> initDeq(windowSize,0.0); // deque to initialise per-site result values
            // vector of five per-site deques:
            // 1) for Fst numerator; 2) Fst denominator;
            // 3) Dxy; 4) set1 heterozygosity; 5) set2 heterozygosity, and
            // 6) is for the coordinates
            std::vector<std::deque<double>> initSixDeques(6,initDeq);
            std::vector<std::vector<std::deque<double>>> scanResults(pairs.size(),initSixDeques);
            resultsSNPwindows = scanResults;
        }
        
        void initResultsPhysicalWindows() {
            std::vector<std::vector<double>> initVectorFixed(5);
            std::vector<std::vector<std::vector<double>>> PBSfixedWindowResults(pairs.size(),initVectorFixed);
            resultsPhysicalWindows = PBSfixedWindowResults;
        }
        
        void initResultsGenes() {
            std::vector<std::vector<double>> initGeneVectors(9); // For the nine PBS columns in the _PBSGenes_ files
            std::vector<std::vector<std::vector<double>>> PBSgeneResults(pairs.size(),initGeneVectors);
            resultsGenes = PBSgeneResults;
        }
};

#endif /* Fst_hpp */
