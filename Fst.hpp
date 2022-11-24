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
#include "UtilsSetCounts.hpp"
#include "UtilsSetInfo.hpp"
#include "UtilsPopulationComparisons.hpp"
#include <math.h>

void parseFstOptions(int argc, char** argv);
int fstMain(int argc, char** argv);
bool bPairInformativeThisSNP(const double p1, const double p2, const int n1, const int n2,
                             const int set1FullSize, const int set2FullSize, const double maxMissing);


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

// Define sizes of scan result containers:
// 1) for Fst numerator; 2) Fst denominator;
// 3) Dxy; 4) set1 heterozygosity; 5) set2 heterozygosity, and
// 6) is for the coordinates
#define NumFstSNPWindowDeques 6
#define NumFstPhysicalWindowVectors 5
#define NumFstGeneResultVectors 6   // Three Fst values per gene: exons, introns, promoter; keeping Fst numerator and Fst denominator for each

class FstPairs : public PopulationComparisons {
public:
    
    FstPairs(const string& FstPairsFileName, const string& runName, const int windowSize, const int windowStep, const int fixedWindowSize, const bool bAnnotationPresent, const bool bMakeRegionsBed, const double regionBedMin, const bool bRecombMapPresent) : PopulationComparisons(FstPairsFileName) {
        
        std::cout << "\nCalculating statistics for the set pair(s):\n";
        string line;
        while (getline(*ComparisonsFile,line)) {
            // std::cerr << line << std::endl;
            std::vector<string> twoPops = split(line, '\t'); assert(twoPops.size() == 2);
            print_vector(twoPops, std::cout);
            std::ofstream* outFile = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_Fst_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
            std::ofstream* outFileFixedWindow = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_" +  "_Fst_" + runName + "_FW" + numToString(fixedWindowSize) + ".txt");
            
            string fstOufileHeader = "chr\twStart\twEnd\tFst\tDxy\t" + twoPops[0] + "_pi\t" + twoPops[1] + "_pi\tAccessible_bp";
            *outFile << fstOufileHeader;
            if (bRecombMapPresent) *outFile << "\t" << "mean_r";
            *outFile << "\n";
            
            *outFileFixedWindow << fstOufileHeader;
            //outFile->setf(std::ios_base::fixed); // Avoid scientific notation in the coordinates
            outFiles.push_back(outFile); outFilesFixedWindow.push_back(outFileFixedWindow);
      
            if (bAnnotationPresent) {
                std::ofstream* outFileGenes = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_Fst_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
                *outFileGenes << "gene\t" << "numSNPsExons\tnumSNPsIntrons\tnumSNPsPromoter\texons\twIntrons\tpromoter\n" << std::endl;
                outFilesGenes.push_back(outFileGenes);
            }
            
            if (bMakeRegionsBed) {
                std::ofstream* outFileRegions = new std::ofstream(twoPops[0] + "_" + twoPops[1] + "_FstAbove_" +  numToString(regionBedMin) + "_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
                outFilesRegions.push_back(outFileRegions);
            }
            
            comparisons.push_back(twoPops);
        }
        
        // Need to prepare the vectors to hold the PBS values and the coordinates:
        initResultsSNPwindows(windowSize, NumFstSNPWindowDeques);
        initResultsPhysicalWindows(NumFstPhysicalWindowVectors);
        if (bAnnotationPresent) initResultsGenes(NumFstGeneResultVectors);
        
        usedVars.resize(comparisons.size(),0);
    };
    
    std::vector<std::vector<string> > getPairs(){ return comparisons;} // name the method whatever you like.

    
    void addSNPresultsToGene(const int trioNumber, double thisSNPFstNumerator,double thisSNPFstDenominator, const Annotation& a) {
        
        if (a.SNPcategory == "exon") {
            resultsGenes.at(trioNumber)[0].push_back(thisSNPFstNumerator);
            resultsGenes.at(trioNumber)[1].push_back(thisSNPFstDenominator);
        } else if (a.SNPcategory == "intron") {
            resultsGenes.at(trioNumber)[2].push_back(thisSNPFstNumerator);
            resultsGenes.at(trioNumber)[3].push_back(thisSNPFstDenominator);
        } else if (a.SNPcategory == "promoter") {
            resultsGenes.at(trioNumber)[4].push_back(thisSNPFstNumerator);
            resultsGenes.at(trioNumber)[5].push_back(thisSNPFstDenominator);
        }
    }
    
    
    void summariseAndOutputPerGene(const int trioNumber, const string& geneName) {
        int nExonSNPs = (int)resultsGenes.at(trioNumber)[0].size();
        int nIntronSNPs = (int)resultsGenes.at(trioNumber)[1].size();
        int nPromoterSNPs = (int)resultsGenes.at(trioNumber)[2].size();
        
        double ExonFst = 0; if (nExonSNPs > 0)
            ExonFst = (double) vector_sum(resultsGenes.at(trioNumber)[0]) / vector_sum(resultsGenes.at(trioNumber)[1]);
        double IntronFst = 0; if (nIntronSNPs > 0)
            IntronFst = (double) vector_sum(resultsGenes.at(trioNumber)[2]) / vector_sum(resultsGenes.at(trioNumber)[3]);
        double PromoterFst = 0; if (nPromoterSNPs > 0)
            PromoterFst = (double) vector_sum(resultsGenes.at(trioNumber)[4]) / vector_sum(resultsGenes.at(trioNumber)[5]);
        
        *outFilesGenes[trioNumber] << geneName << "\t" << nExonSNPs << "\t" << nIntronSNPs << "\t" << nPromoterSNPs << "\t" << ExonFst << "\t" << IntronFst << "\t" << PromoterFst << std::endl;
        
        for (int j = 0; j < NumFstGeneResultVectors; j++) {
            resultsGenes.at(trioNumber)[j].clear();
            resultsGenes.at(trioNumber)[j].shrink_to_fit();
        }
    }
    
    void addSNPresultsToWindows(const int pairNumber, const double thisSNPFstNumerator, const double thisSNPFstDenominator, const double thisSNPDxy, const double thisSNPpi1, const double thisSNPpi2, const int SNPcoordinate) {
        
        resultsSNPwindows[pairNumber][0].push_back(thisSNPFstNumerator); resultsSNPwindows[pairNumber][0].pop_front();
        resultsSNPwindows[pairNumber][1].push_back(thisSNPFstDenominator); resultsSNPwindows[pairNumber][1].pop_front();
        resultsSNPwindows[pairNumber][2].push_back(thisSNPDxy); resultsSNPwindows[pairNumber][2].pop_front();
        resultsSNPwindows[pairNumber][3].push_back(thisSNPpi1); resultsSNPwindows[pairNumber][3].pop_front();
        resultsSNPwindows[pairNumber][4].push_back(thisSNPpi2); resultsSNPwindows[pairNumber][4].pop_front();
        resultsSNPwindows[pairNumber][5].push_back((double)SNPcoordinate); resultsSNPwindows[pairNumber][5].pop_front();
        
        resultsPhysicalWindows[pairNumber][0].push_back(thisSNPFstNumerator);
        resultsPhysicalWindows[pairNumber][1].push_back(thisSNPFstDenominator);
        resultsPhysicalWindows[pairNumber][2].push_back(thisSNPDxy);
        resultsPhysicalWindows[pairNumber][3].push_back(thisSNPpi1);
        resultsPhysicalWindows[pairNumber][4].push_back(thisSNPpi2);
    }
    
    void finalizeAndOutputSNPwindow(const int pairNumber, const string chr, const int currentSNPcoordinate, const AccessibleGenome* ag, const RecombinationMap* r, const bool bZeroRounding) {
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
            
            double Fst = calculateFst(resultsSNPwindows[pairNumber][iFstNum], resultsSNPwindows[pairNumber][iFstDenom], bZeroRounding);
            double Dxy = vector_average_withRegion(resultsSNPwindows[pairNumber][iDxy], accessibleInThisWindow);
            double Pi1 = vector_average_withRegion(resultsSNPwindows[pairNumber][iPi1], accessibleInThisWindow);
            double Pi2 = vector_average_withRegion(resultsSNPwindows[pairNumber][iPi2], accessibleInThisWindow);
            
            *outFiles[pairNumber] << chr << "\t" << (int)startCoord << "\t" << currentSNPcoordinate << "\t" << Fst << "\t" << Dxy << "\t" << Pi1 << "\t" << Pi2 << "\t" << accessibleInThisWindow;
            
            if (r->initialised) {
                double meanRecomb = r->getMeanRecombinationRate(chr, startCoord, currentSNPcoordinate);
                // std::cerr << "meanRecomb: " << meanRecomb << std::endl;
                *outFiles[pairNumber] << "\t" << meanRecomb;
            }
            
            *outFiles[pairNumber] << std::endl;
        }
    }
    
    
    void finalizeAndOutputPhysicalWindow(const int pairNumber, const int physicalWindowSize, const string chr, const int currentSNPcoord, const AccessibleGenome* ag, int& thisWindowStart, int& thisWindowEnd, const bool bZeroRounding) {
        
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
                thisFixedWindowFst = calculateFst(resultsPhysicalWindows[pairNumber][iFstNum], resultsPhysicalWindows[pairNumber][iFstDenom],bZeroRounding);
                thisFixedWindowDxy = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iDxy], accessibleInThisWindow);
                thisFixedWindowPi1 = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iPi1], accessibleInThisWindow);
                thisFixedWindowPi2 = vector_average_withRegion(resultsPhysicalWindows[pairNumber][iPi2], accessibleInThisWindow);
                *outFilesFixedWindow[pairNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\t" << thisFixedWindowFst << "\t" << thisFixedWindowDxy << "\t" << thisFixedWindowPi1 << "\t" << thisFixedWindowPi2 << "\t" << accessibleInThisWindow << std::endl;
            } else { // If this physical window did not have any SNPs, then Fst="NA"
                *outFilesFixedWindow[pairNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\tNA\t" << thisFixedWindowDxy << "\t" << thisFixedWindowPi1 << "\t" << thisFixedWindowPi2 << "\t" << accessibleInThisWindow << std::endl;
            }
            
            for (int i = 0; i < resultsPhysicalWindows[pairNumber].size(); i++) {
                resultsPhysicalWindows[pairNumber][i].clear();
                resultsPhysicalWindows[pairNumber][i].shrink_to_fit();
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
    
};

#endif /* Fst_hpp */
