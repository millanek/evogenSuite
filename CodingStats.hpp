//
//  CodingStats.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 28.06.23.
//

#ifndef CodingStats_hpp
#define CodingStats_hpp

#include "UtilsGeneral.hpp"
#include "UtilsSetInfo.hpp"
#include "UtilsCodingStats.hpp"
#include "UtilsStats.hpp"

void parseCodingStatsOptions(int argc, char** argv);
int getCodingStats(int argc, char** argv);

bool isUnrecognisedCodon(string codon);
bool compareFirstInPair(std::pair<int, string> a, std::pair<int, string> b);

class CodonSetCounts {
public:
    CodonSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap) {
        
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setNumValidCodons[it->first]; // Default construct an empty vector there
            setMissingnessProportion[it->first]; // Default construct an empty vector there
            setAAFs[it->first]; // Default construct an empty vector there
        }
    };
    
    // ------------- FIELDS ---------------
      
    // Alternate allele frequencies in each set (for each alternative allele)
    string codon1; string codon2; double pN_jk; double pS_jk;
    std::map<string, double> setAAFs;
    std::map<string, int> setNumValidCodons;
    std::map<string, double> setMissingnessProportion;
    
    
    
    void fillCodonCounts(const SetInformation& setInfo, const std::vector<string>& altCodons, const std::vector<std::pair<int, string>>& uniqueCodonsAndCounts) {
        
        codon1 = uniqueCodonsAndCounts[0].second;
        codon2 = uniqueCodonsAndCounts[1].second;
        
        for (int j = 0; j < setInfo.populations.size(); j++) {
            int c_A1 = 0; int c_A2 = 0;
            for (std::vector<std::string>::size_type i = 0; i != altCodons.size(); i++) {
                try { if (setInfo.posToPopMap.at(i) == setInfo.populations[j]) {
                    if (altCodons[i] == codon1) c_A1++;
                    else if (altCodons[i] == codon2) c_A2++;
                } } catch (const std::out_of_range& e) {
                    std::cerr << "ERROR: Couldn't find sample number from the alignment " << i << " in the SETS file" << std::endl;
                    exit(1);
                }
            }
           // std::cerr << "c_A1: " << c_A1 << std::endl;
           // std::cerr << "c_A2: " << c_A2 << std::endl;
            setNumValidCodons[setInfo.populations[j]] = c_A1 + c_A2;
            setAAFs[setInfo.populations[j]] = (double)c_A2 / (c_A1 + c_A2);
            setMissingnessProportion[setInfo.populations[j]] = 1 - ((c_A1 + c_A2)/(double)setInfo.popToPosMap.at(setInfo.populations[j]).size());
        }
    }
};


class MultipleAlignment {
public:
    MultipleAlignment(const string& fN, int numCombinations) {
        filename = fN;
        readAlignment();
        numSequences = (int)allSeqs.size();
        alignmentLengthNucleotides = (int)allSeqs[0].length();
        validateAlignment();
        altCodons.assign(numSequences, ""); haveStop.assign(numSequences, 0); isUnrecognised.assign(numSequences, 0);
        pairwiseMatrices = new CDSComparisonMatrices(numSequences);
        for (int j = 0; j != numSequences; j++) { altCodons[j].reserve(3);}
        
        N_highFst.resize(numCombinations, 0);
        N_lowFst.resize(numCombinations, 0);
        S_highFst.resize(numCombinations, 0);
        S_lowFst.resize(numCombinations, 0);
        
    }
    
    // To DO:
    // 1) Get some sort of inbreeding coefficient (how often are they homozygous ('fixed' in a species) vs. heterozygous)
    // 2) Would also be good to distinguish derived vs. ancestral alleles?
 // void getStatsAllSequences(const SetInformation* sets, const double tStVratio, bool nonCodingNull) {
    void getStatsAllSequences(const double tStVratio, bool nonCodingNull, const SetInformation& setInfo, const double FstThreshold) {
       // std::cerr << "Collecting gene sequence statistics...." << std::endl;
        clock_t begin = clock();
        
        for (string::size_type i = 0; i != alignmentLengthNucleotides; i=i+3) {
            for (std::vector<std::string>::size_type j = 0; j != allSeqs.size(); j++) {
                altCodons[j] = allSeqs[j].substr(i,3);
                string AA = getAminoAcid(altCodons[j]);
                // If this is an alignment for a coding region:
                // the remainder of the sequence after any stop will be excluded from the calculations for any sequence
                if (AA == "Stop") haveStop[j] = 1;
                if (AA == "Unrecognised codon") isUnrecognised[j] = 1;
            }
            
            // If there is variation at the codon level, we find the unique codons
            if (std::adjacent_find( altCodons.begin(), altCodons.end(), std::not_equal_to<string>() ) != altCodons.end() ) {
                
                addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(altCodons,haveStop,isUnrecognised, pairwiseMatrices);
                std::vector<string> uniqueValidCodons = findUniqueCodons();
                if (uniqueValidCodons.size() <= 1) continue;
                fillAndSortUniqueCodonCounts(uniqueValidCodons);
                
                CodonSetCounts* c = new CodonSetCounts(setInfo.popToPosMap);
                c->fillCodonCounts(setInfo, altCodons, uniqueCodonsAndCounts);
                
                int d = getCodonDistance(c->codon1, c->codon2);
                double N_ijk = calculateN(c->codon1,c->codon2, d, false);
                double n_d_ijk = calculateNd(c->codon1, c->codon2, d);
                double s_d_ijk = d - n_d_ijk;
                
                double N_tS = calculateNtS(c->codon1, c->codon2, d, false);
                double S_tS = 1 - N_tS;
                double tV_N_jk = N_ijk - N_tS;
                double tV_S_jk = 2 - (N_ijk - N_tS);
                
                
                double pN_jk = n_d_ijk/( (2*tStVratio*N_tS) + tV_N_jk);
                double pS_jk = s_d_ijk/( (2*tStVratio*S_tS) + tV_S_jk);
                
                //int d = getCodonDistance(altCodons[j],altCodons[k]);
               // double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
                //double s_d_ijk = d - n_d_ijk;
                //pairwiseMatrices->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_N_jk[j][k])+pairwiseMatrices->tV_N_jk[j][k]);
                // pairwiseMatrices->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_S_jk[j][k])+pairwiseMatrices->tV_S_jk[j][k])
                //p->N_d_jk[j][k] = p->N_d_jk[j][k] + n_d_ijk;
                // p->S_d_jk[j][k] = p->S_d_jk[j][k] + s_d_ijk;

               // N_ijk = calculateN(altCodons[j],altCodons[k], d, false);
               // N_tS = calculateNtS(altCodons[j],altCodons[k], d, false);
            
              //  S_ijk = (3 - N_ijk);
             //   p->N_jk[j][k] = p->N_jk[j][k] + N_ijk; p->S_jk[j][k] = p->S_jk[j][k] + S_ijk;
             //   p->tS_N_jk[j][k] = p->tS_N_jk[j][k] + N_tS;
              //  p->tS_S_jk[j][k] = p->tS_S_jk[j][k] + (1 - N_tS);
              //  p->tV_N_jk[j][k] = p->tV_N_jk[j][k] + (N_ijk - N_tS);
               // p->tV_S_jk[j][k] = p->tV_S_jk[j][k] + (2 - (N_ijk - N_tS));
                
                
                int iCombination = 0;
                for (int ii = 0; ii < setInfo.populations.size(); ii++) {
                    for (int jj = ii+1; jj < setInfo.populations.size(); jj++) {
                        double p1 = c->setAAFs.at(setInfo.populations[ii]);
                        double p2 = c->setAAFs.at(setInfo.populations[jj]);
                     //   std::cerr << "p1: " << p1 << std::endl;
                        
                        double thisFstNumerator = calculateFstNumerator(p1, p2, c->setNumValidCodons.at(setInfo.populations[ii]), c->setNumValidCodons.at(setInfo.populations[jj]));
                        double thisFstDenominator = calculateFstDenominator(p1, p2);
                        double thisFst = thisFstNumerator/thisFstDenominator;
                        if (thisFst < FstThreshold) {
                            N_lowFst[iCombination] += pN_jk; S_lowFst[iCombination] += pS_jk;
                        } else {
                            N_highFst[iCombination] += pN_jk; S_highFst[iCombination] += pS_jk;
                        }
                        iCombination++;
                    }
                }
              //  std::cerr << "thisFst: " << thisFst << std::endl;
                //setInfo.posToPopMap[j];
            }
            
            uniqueCodonsAndCounts.clear();
            altCodons.assign(numSequences, ""); isUnrecognised.assign(numSequences, 0);
        }
        
        double sumPn = 0; double sumPs = 0;
        std::vector<double> pNjackknifeVector;
        std::vector<double> pSjackknifeVector;
        std::vector<double> pN_pSjackknifeVector;
        for (std::vector<std::string>::size_type j = 0; j != numSequences - 1; j++) {
            //std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
            for (std::vector<std::string>::size_type k = j+1; k != numSequences; k++) {
                //double pN_jk = pairwiseMatrices.H1p->N_d_jk[j][k]/pairwiseMatrices.H1p->N_jk[j][k];
                double pN_jk = pairwiseMatrices->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_N_jk[j][k])+pairwiseMatrices->tV_N_jk[j][k]);
             //   validateCalculation(pN_jk,(int)j,(int)k,sumPn);
                sumPn = sumPn + pN_jk;
                double pS_jk = pairwiseMatrices->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_S_jk[j][k])+pairwiseMatrices->tV_S_jk[j][k]);
                sumPs = sumPs + pS_jk;
        
                // make sure each sample is used only once in the jackknife calculation
                // to keep the observations independent
                if ((j % 2 == 0) && (k == j+1)) {
                    pNjackknifeVector.push_back(pN_jk); pSjackknifeVector.push_back(pS_jk);
                    pN_pSjackknifeVector.push_back(pN_jk-pS_jk);
                }
            }
        }
        pN = (2.0/(numSequences*(numSequences-1)))*sumPn;
        pS = (2.0/(numSequences*(numSequences-1)))*sumPs;
        
        pNjackknifeVectorSize = (int)pNjackknifeVector.size();
        
        
        if (pNjackknifeVectorSize > 10) {
            pNstdErr = jackknive_std_err(pNjackknifeVector);
            pSstdErr = jackknive_std_err(pSjackknifeVector);
            pNpSstdErr = jackknive_std_err(pN_pSjackknifeVector);
           
        }
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       // std::cerr << "Time per gene: " << elapsed_secs << std::endl;
    }
    
    void printAlignmentStats(std::ostream& outFile) {
        outFile << filename << "\t";
        outFile << alignmentLengthNucleotides << "\t";
        outFile << pN << "\t";
        outFile << pS << "\t";
        outFile << pNstdErr << "\t";
        outFile << pSstdErr << "\t";
        outFile << pNpSstdErr << std::endl;
    }
    
    void printFstStats(std::ostream& outFile, int iCombination) {
        outFile << filename << "\t";
        outFile << alignmentLengthNucleotides << "\t";
        outFile << N_lowFst[iCombination] << "\t";
        outFile << S_lowFst[iCombination] << "\t";
        outFile << N_highFst[iCombination] << "\t";
        outFile << S_highFst[iCombination] << std::endl;
    }
    
    string filename;
    std::vector<string> allSeqs;
    int numSequences; std::vector<string> sequenceNames;
    int alignmentLengthNucleotides;
    double pN; double pS;
    double pNstdErr; double pSstdErr; double pNpSstdErr;
    int pNjackknifeVectorSize;
    
    std::vector<double> N_highFst; std::vector<double> N_lowFst;
    std::vector<double> S_highFst; std::vector<double> S_lowFst;
    
    
private:
    
    std::vector<string> altCodons; std::vector<std::pair<int, string>> uniqueCodonsAndCounts;
    std::vector<int> haveStop; std::vector<int> isUnrecognised;
    CDSComparisonMatrices* pairwiseMatrices;
    
    
    void readAlignment() {
        std::ifstream* alignmentFile = new std::ifstream(filename.c_str());
        assertFileOpen(*alignmentFile, filename);

        std::string line;
        getline(*alignmentFile, line); assert(line[0] == '>');
        string seqName = line.substr(1);
        sequenceNames.push_back(seqName);
        std::string seqThisIndividual = ""; seqThisIndividual.reserve(10000);
        while (getline(*alignmentFile, line)) {
            if (line[0] != '>') {
                seqThisIndividual.append(line);
            } else {
                seqName = line.substr(1);
                sequenceNames.push_back(seqName);
                allSeqs.push_back(seqThisIndividual);
                seqThisIndividual = "";
            }
        } alignmentFile->close();
        
        allSeqs.push_back(seqThisIndividual);
    }
    
    void validateAlignment() {
        if (allSeqs.size() == 0) {
            std::cerr << "Error: Did not succeed in reading any sequences from the file " << filename << std::endl;
            std::cerr << "Error: Please check - this alignment file is probably empty." << std::endl;
            assert(allSeqs.size() > 0);
        }
        
        if (allSeqs[0].length() % 3 != 0) {
            std::cerr << "Error: The length of the first sequence is not divisible by three " << std::endl;
            std::cerr << "Error: This tool requires coding nucleotide sequence alignment." << std::endl;
            assert(allSeqs[0].length() % 3 == 0); // The gene length must be divisible by three
    
        }
        
        for (std::vector<string>::size_type i = 0; i != numSequences; i++) {
            if (allSeqs[i].length() != alignmentLengthNucleotides) {
                std::cerr << "Error: The length of the first and the " << i + 1 << "th sequence differs." << std::endl;
                std::cerr << "Error: All sequences should be equal length (i.e., indels removed)." << std::endl;
                assert(allSeqs[i].length() == alignmentLengthNucleotides);
            }
        }
    }
    
    void validateCalculation(const double pN_jk, const int j, const int k, const double sumPn) {
        if (std::isnan(pN_jk) || pairwiseMatrices->tS_N_jk[j][k] == 0) {
        std::cerr << "j = " << j << "; k = " << k << std::endl;
        std::cerr << "N_d_jk[j][k] = " << pairwiseMatrices->N_d_jk[j][k] << "; tS_N_jk[j][k] = " << pairwiseMatrices->tS_N_jk[j][k] << std::endl;
        std::cerr << "pN_jk = " << pN_jk << "; sumPn = " << sumPn << std::endl;
            print_matrix(pairwiseMatrices->N_d_jk, std::cerr);
            print_matrix(pairwiseMatrices->S_d_jk, std::cerr);
        }
    }
    
    std::vector<string> findUniqueCodons() {
        std::vector<string> altCodonsCopy = altCodons;
        std::sort(altCodonsCopy.begin(), altCodonsCopy.end());
        vector<string>::iterator it = unique(altCodonsCopy.begin(), altCodonsCopy.end());
        altCodonsCopy.resize(distance(altCodonsCopy.begin(),it));
        it = std::remove_if(altCodonsCopy.begin(), altCodonsCopy.end(), isUnrecognisedCodon);
        altCodonsCopy.resize(distance(altCodonsCopy.begin(),it));
        return altCodonsCopy;
    }
    
    void fillAndSortUniqueCodonCounts(const std::vector<string>& uniqueValidCodons) {
        for (std::vector<std::string>::size_type j = 0; j != uniqueValidCodons.size(); j++) {
            int codonCount = (int)count(altCodons.begin(), altCodons.end(), uniqueValidCodons[j]);
           // std::cerr << "codonCount: " << codonCount << std::endl;
            uniqueCodonsAndCounts.push_back(std::make_pair(codonCount, uniqueValidCodons[j]));
        }
        
        // Sort the vector of pairs
        std::sort(uniqueCodonsAndCounts.begin(), uniqueCodonsAndCounts.end(), compareFirstInPair);
    }
    
};


        /*   //  std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
             if (sets->initialised == true) {
             //    std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
                 if ((sets->set1Loci.count(j) == 1 && sets->set1Loci.count(k) == 1) || (sets->set1Loci.count(k) == 1 && sets->set1Loci.count(j) == 1)) {
                     sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                     sets->withinSet1pN = sets->withinSet1pN + pN_jk + H2pN_jk;
                 } else if ((sets->set2Loci.count(j) == 1 && sets->set2Loci.count(k) == 1) || (sets->set2Loci.count(k) == 1 && sets->set2Loci.count(j) == 1)) {
                     sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                     sets->withinSet2pN = sets->withinSet2pN + pN_jk + H2pN_jk;
                 } else if ((sets->set1Loci.count(j) == 1 && sets->set2Loci.count(k) == 1)) {
                     sets->set1vsSet2pN = sets->set1vsSet2pN + pN_jk + H2pN_jk;
                     sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                 } else if (sets->set1Loci.count(k) == 1 && sets->set2Loci.count(j) == 1) {
                     sets->set1vsSet2pN = sets->set1vsSet2pN + pN_jk + H2pN_jk;
                     sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                 } else if ((sets->set1Loci.count(j) == 1 || sets->set2Loci.count(j) == 1) && sets->set3Loci.count(k) == 1) {
                     sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + pN_jk + H2pN_jk;
                 } else if ((sets->set1Loci.count(k) == 1 || sets->set2Loci.count(k) == 1) && sets->set3Loci.count(j) == 1) {
                     sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + pN_jk + H2pN_jk;
                     //std::cerr << "j = " << j << "; (set 3): " << allSeqs[j] << std::endl;
                     //std::cerr << "k = " << k << "; (set 1 or 2): " << allSeqs[k] << std::endl;
                     //std::cerr << std::endl;
                 } else if ((sets->set3Loci.count(j) == 1 && sets->set3Loci.count(k) == 1) || (sets->set3Loci.count(k) == 1 && sets->set3Loci.count(j) == 1)) {
                     sets->withinSet3pN = sets->withinSet3pN + pN_jk + H2pN_jk;
                 }
             }
         }
     } */
        
        
        
#endif /* CodingStats_hpp */
