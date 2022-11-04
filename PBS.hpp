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
#include "UtilsSetCounts.hpp"
#include "UtilsSetInfo.hpp"
#include "UtilsAnnotation.hpp"
#include "UtilsPopulationComparisons.hpp"

void parsePBSoptions(int argc, char** argv);
int PBSmain(int argc, char** argv);
bool bTrioInformativeThisSNP(const double p1, const double p2, const double p3,
                             const int n1, const int n2, const int n3,
                             const string set1, const string set2, const string set3,
                             const SetInformation& setInfo,
                             const double maxMissing);


#define NumPBS_SNPWindowDeques 4
#define NumPBS_PhysicalWindowVectors 3
#define NumPBS_GeneResultVectors 9

class PBStrios : public PopulationComparisons {
public:
    
    PBStrios(const string& PBStriosFileName, const string& runName, const int windowSize, const int windowStep, const int fixedWindowSize, const bool bAnnotationPresent) : PopulationComparisons(PBStriosFileName) {
        
        string line;
        while (getline(*ComparisonsFile,line)) {
            // std::cerr << line << std::endl;
            std::vector<string> threePops = split(line, '\t'); assert(threePops.size() == 3);
            std::ofstream* outFile = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBS_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
            std::ofstream* outFileFixedWindow = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBS_" + runName + "_FW" + numToString(fixedWindowSize) + ".txt");
            *outFile << "chr\twStart\twEnd\t" << threePops[0] << "\t" << threePops[1] << "\t" << threePops[2] << std::endl;
            *outFileFixedWindow << "chr\twStart\twEnd\t" << threePops[0] << "\t" << threePops[1] << "\t" << threePops[2] << "\t" << "nFwSNPs1" << "\t" << "nFwSNPs2" << "\t" << "nFwSNPs3" << std::endl;
            outFile->setf(std::ios_base::fixed); // Avoid scientific notation in the coordinates
            outFiles.push_back(outFile); outFilesFixedWindow.push_back(outFileFixedWindow);
            if (bAnnotationPresent) {
                std::ofstream* outFileGenes = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBSGenes_" + runName + "_" + numToString(windowSize) + "_" + numToString(windowStep) + ".txt");
                *outFileGenes << "gene\t" << "numSNPsExons\tnumSNPsIntrons\tnumSNPs3kbPromoter\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons\t" << threePops[0] << "_promoter\t" << threePops[1] << "_promoter\t" << threePops[2] << "_promoter" << std::endl;
                outFilesGenes.push_back(outFileGenes);
            }
            comparisons.push_back(threePops);
        }
        
        // Need to prepare the vectors to hold the PBS values and the coordinates:
        initResultsSNPwindows(windowSize, NumPBS_SNPWindowDeques);   // Creating four deques for SNP window results
        initResultsPhysicalWindows(NumPBS_PhysicalWindowVectors);          // Creating three vectors for physical window results
        if (bAnnotationPresent) initResultsGenes(NumPBS_GeneResultVectors);
        
        usedVars.resize(comparisons.size(),0);
    };
    
    
    std::vector<std::vector<string> > getTrios(){ return comparisons;} // name the method whatever you like.
        
    void addSNPresultsToWindows(const int trioNumber, const std::vector<double>& thisSNP_PBS, const int SNPcoordinate) {
        
        for (int i = 0; i < (NumPBS_SNPWindowDeques - 1); i++)
        resultsSNPwindows.at(trioNumber)[i].push_back(thisSNP_PBS[i]);
        resultsSNPwindows[trioNumber][(NumPBS_SNPWindowDeques - 1)].push_back((double)SNPcoordinate);
        
        for (int i = 0; i < NumPBS_SNPWindowDeques; i++) resultsSNPwindows.at(trioNumber)[i].pop_front();
        
        for (int i = 0; i < NumPBS_PhysicalWindowVectors; i++)
            resultsPhysicalWindows.at(trioNumber)[i].push_back(thisSNP_PBS[i]);
    }
    
    void addSNPresultsToGene(const int trioNumber, const std::vector<double>& thisSNP_PBS, const Annotation& a) {
        
        if (a.SNPcategory == "exon") {
            resultsGenes[trioNumber][0].push_back(thisSNP_PBS[0]);
            resultsGenes[trioNumber][1].push_back(thisSNP_PBS[1]);
            resultsGenes[trioNumber][2].push_back(thisSNP_PBS[2]);
        } else if (a.SNPcategory == "intron") {
            resultsGenes[trioNumber][3].push_back(thisSNP_PBS[0]);
            resultsGenes[trioNumber][4].push_back(thisSNP_PBS[1]);
            resultsGenes[trioNumber][5].push_back(thisSNP_PBS[2]);
        } else if (a.SNPcategory == "promoter") {
            resultsGenes[trioNumber][6].push_back(thisSNP_PBS[0]);
            resultsGenes[trioNumber][7].push_back(thisSNP_PBS[1]);
            resultsGenes[trioNumber][8].push_back(thisSNP_PBS[2]);
        }
        
    }
    
    void finalizeAndOutputPhysicalWindow(const int trioNumber, const int physicalWindowSize, const string chr, const int currentSNPcoord, int& thisWindowStart, int& thisWindowEnd) {
        
        if (thisWindowStart < currentSNPcoord) { // As long as the current SNP coord is more than the previous (i.e., we are on the same chromosome; I should check this properly
            int nFwSNPs1 = (int)resultsPhysicalWindows.at(trioNumber)[0].size();
            int nFwSNPs2 = (int)resultsPhysicalWindows.at(trioNumber)[1].size();
            int nFwSNPs3 = (int)resultsPhysicalWindows.at(trioNumber)[2].size();
            
            if ( (nFwSNPs1 > 0) && (nFwSNPs2 > 0) && (nFwSNPs3 > 0) ) {
                double PBSfw1 = vector_average(resultsPhysicalWindows.at(trioNumber)[0]);
                double PBSfw2 = vector_average(resultsPhysicalWindows.at(trioNumber)[1]);
                double PBSfw3 = vector_average(resultsPhysicalWindows.at(trioNumber)[2]);
                *outFilesFixedWindow[trioNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\t" << PBSfw1 << "\t" << PBSfw2 << "\t" << PBSfw3 << "\t" << nFwSNPs1 << "\t" << nFwSNPs2 << "\t" << nFwSNPs3 << std::endl;
            } else { // If this physical window did not have any SNPs, then Fst="NA"
                *outFilesFixedWindow[trioNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\tNA\tNA\tNA\t0\t0\t0\n";
            }
            for (int i = 0; i < resultsPhysicalWindows.at(trioNumber).size(); i++) {
                resultsPhysicalWindows[trioNumber][i].clear();
                resultsPhysicalWindows[trioNumber][i].shrink_to_fit();
            }
        
            thisWindowStart = thisWindowStart + physicalWindowSize;
            thisWindowEnd = thisWindowEnd + physicalWindowSize;
            
            while (!(thisWindowStart <= currentSNPcoord && thisWindowEnd > currentSNPcoord)) {
                *outFilesFixedWindow[trioNumber] << chr << "\t" << thisWindowStart << "\t" << thisWindowEnd << "\tNA\tNA\tNA\t0\t0\t0\n";
                thisWindowStart = thisWindowStart + physicalWindowSize;
                thisWindowEnd = thisWindowEnd + physicalWindowSize;
            }
        } else {
            thisWindowStart = 0;
            thisWindowEnd = 0 + physicalWindowSize;
        }
    }
    
    void summariseAndOutputPerGene(const int trioNumber, const string& geneName) {
        
        int nExonSNPs = (int)resultsGenes.at(trioNumber)[0].size();
        int nIntronSNPs = (int)resultsGenes.at(trioNumber)[3].size();
        int nPromoterSNPs = (int)resultsGenes.at(trioNumber)[6].size();
        
        double ExonPBS1 = 0; double ExonPBS2 = 0; double ExonPBS3 = 0;
        if (nExonSNPs > 0) {
            ExonPBS1 = vector_average(resultsGenes.at(trioNumber)[0]);
            ExonPBS2 = vector_average(resultsGenes.at(trioNumber)[1]);
            ExonPBS3 = vector_average(resultsGenes.at(trioNumber)[2]);
        }
        
        double IntronPBS1 = 0; double IntronPBS2 = 0; double IntronPBS3 = 0;
        if (nIntronSNPs > 0) {
            IntronPBS1 = vector_average(resultsGenes.at(trioNumber)[3]);
            IntronPBS2 = vector_average(resultsGenes.at(trioNumber)[4]);
            IntronPBS3 = vector_average(resultsGenes.at(trioNumber)[5]);
        }
        
        double PromoterPBS1 = 0; double PromoterPBS2 = 0; double PromoterPBS3 = 0;
        if (nPromoterSNPs > 0) {
            PromoterPBS1 = vector_average(resultsGenes.at(trioNumber)[6]);
            PromoterPBS2 = vector_average(resultsGenes.at(trioNumber)[7]);
            PromoterPBS3 = vector_average(resultsGenes.at(trioNumber)[8]);
        }
        *outFilesGenes[trioNumber] << geneName << "\t" << nExonSNPs << "\t" << nIntronSNPs << "\t" << nPromoterSNPs << "\t" << ExonPBS1 << "\t" << ExonPBS2 << "\t" << ExonPBS3 << "\t" << IntronPBS1 << "\t" << IntronPBS2 << "\t" << IntronPBS3 << "\t" << PromoterPBS1 << "\t" << PromoterPBS2 << "\t" << PromoterPBS3 << std::endl;
        
        for (int j = 0; j < NumPBS_GeneResultVectors; j++) {
            resultsGenes.at(trioNumber)[j].clear();
            resultsGenes.at(trioNumber)[j].shrink_to_fit();
        }
    }
    
};





#endif /* PBS_hpp */
