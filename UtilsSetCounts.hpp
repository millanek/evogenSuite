//
//  UtilsSetCounts.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 29.09.22.
//

#ifndef UtilsSetCounts_hpp
#define UtilsSetCounts_hpp

#include "UtilsGeneral.hpp"
#include "UtilsSetInfo.hpp"

class VariantInfo {
public:
    VariantInfo(const std::vector<string>& VCFfields) {
        chr = VCFfields[0]; posInt = atoi(VCFfields[1].c_str());
        refAllele = VCFfields[3]; altAlleles = split(VCFfields[4], ',');
        
        if (refAllele.length() > 1) onlyIndel = true;
        else SNPAlleleIndices.push_back(0);
        
        std::vector<std::string>::iterator it = std::find(altAlleles.begin(), altAlleles.end(), "*");
        if (it != altAlleles.end()) starPos = (int) std::distance(altAlleles.begin(), it);
        else starPos = -1; // There is no star among the alternative alleles
        
        for (int i = 0; i < altAlleles.size(); i++) {
            if (altAlleles[i].length() == 1 && i != starPos) SNPAlleleIndices.push_back(i+1);
        }
        if (SNPAlleleIndices.size() == 0) onlyIndel = true;
        
    }
    
    string chr;
    int posInt;
    string refAllele;
    std::vector<string> altAlleles;
    std::vector<int> SNPAlleleIndices;
    
    bool onlyIndel = false;
    
private:
    int starPos;
};


class GeneralSetCounts {
public:
    GeneralSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) {
        
        for (const auto& [setName, positions] : setsToPosMap) {
            setRefCounts[setName] = 0;
            setAlleles[setName]; // Default construct an empty vector there
            setAltAlleleCounts[setName]; // Default construct an empty vector there
            setAAFs[setName]; // Default construct an empty vector there
            setHetCounts[setName] = 0;
            setSNPPairCounts[setName] = 0;
            
            piPerVariantPerSet[setName];
            hetPerVariantPerSet[setName];
            
            //  setDAFs[setName] = -1.0;
        }
    };
    
    // ------------- FIELDS ---------------
    //  Alleles in each set, encoded simply as 0s, 1s, 2s, 3s as in VCF GT fields
    std::map<string, std::vector<int>> setAlleles;
    
    // Counts of the reference allele and each of the alternate alleles
    std::map<string,int> setRefCounts;
    std::map<string,std::vector<int>> setAltAlleleCounts;
    
    // Num
    std::map<string, double> setFirstAlleleMissingness;
    
    // Alternate allele frequencies in each set (for each alternative allele)
    std::map<string,std::vector<double>> setAAFs;
    
    // Pi and heterozygosity per set for this SNP
    std::map<string, double> piPerVariantPerSet;
    std::map<string, double> hetPerVariantPerSet;
    
    // int AAint;
    // std::map<string,double> setDAFs; double averageDAF;// Allele frequencies - derived allele
    
    // ------------- MEMBER FUNCTIONS ---------------
    void getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap, const VariantInfo& v);
    void getAlternateAlleleFreqFromOutgroup(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    
    void calculatePiPerVariantPerSet();
    void calculateHeterozygosityPerVariantPerSet();
    
    int getOverallCountOfFirstAltAllele() const {
        int overallCount = 0;
        for (const auto& [setName, altAlleleCounts] : setAltAlleleCounts) {
            overallCount += altAlleleCounts[0];
        }
        return overallCount;
    }
    
    int getOverallCountOfRefAllele() const {
        int overallCount = 0;
        for (const auto& [setName, refAlleleCount] : setRefCounts) {
            overallCount += refAlleleCount;
        }
        return overallCount;
    }
    
    void fillMissingFistAlleleInfo(const SetInformation& setInfo) {
        for (const auto& [setName, positions] : setInfo.popToPosMap) {
            int setFullSize = ploidy*(int)positions.size();
            int nFirstAllele = setRefCounts.at(setName) + setAltAlleleCounts.at(setName)[0];
            
            setFirstAlleleMissingness[setName] = 1 - ((double)nFirstAllele/setFullSize);
        }
    }
    
private:
    // These are needed to calculate the heterozygosity values
    std::map<string,int> setHetCounts;
    std::map<string,int> setSNPPairCounts;
    
protected:
    // This is needed to use likelihoods (by the GeneralSetCountsWithLikelihoods class)
    double averageFirstAAF;    // Mean allele frequency (across all sets) of the first alternative SNP allele
};



class GeneralSetCountsWithLikelihoods : public GeneralSetCounts {
public:
    GeneralSetCountsWithLikelihoods(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : GeneralSetCounts(setsToPosMap,nSamples), likelihoodsProbabilitiesType(LikelihoodsProbabilitiesAbsent) {
        for (const auto& [setName, positions] : setsToPosMap) {
            setHWEpriorsFromAAFfromGT[setName].assign(3, -1.0);
            setHWEpriorsFromDAFfromGT[setName].assign(3, -1.0);
            setAAFsFromLikelihoods[setName] = -1.0; setDAFsFromLikelihoods[setName] = -1.0;
            setAlleleProbCounts[setName] = 0;
        }
    };
    
    int likelihoodsProbabilitiesType;
    std::map<string,std::vector<double> > setHWEpriorsFromAAFfromGT;
    std::map<string,std::vector<double> > setHWEpriorsFromDAFfromGT;
    std::map<string,double> setAAFsFromLikelihoods; double averageAAFFromLikelihoods; // Allele frequencies - alternative allele
    std::map<string,double> setDAFsFromLikelihoods; double averageDAFFromLikelihoods;// Allele frequencies - derived allele
    std::map<string,int> setAlleleProbCounts; // The number of non-missing alleles for this set in terms of likelihoods/probabilities

    void setHWEpriorsFromAFfromGT();
    std::vector<double> probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species);
    void getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition);
    int checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields);
    int returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag);
    
};


class GeneralSetCountsWithComplements : public GeneralSetCounts {
    public:
    GeneralSetCountsWithComplements(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples, const string thisSNPchr, const string coord) : GeneralSetCounts(setsToPosMap,nSamples) {
        for (const auto& [setName, positions] : setsToPosMap) {
            setAAFsComplement[setName] = -1.0; setDAFsComplement[setName] = -1.0; setAlleleCountsComplement[setName] = 0;
        }
    }
    std::map<string,double> setAAFsComplement; // Allele frequencies - alternative allele, in the complement of the set
    std::map<string,double> setDAFsComplement; // Allele frequencies - derived allele, in the complement of the set
    std::map<string,int> setAlleleCountsComplement; // The number of non-missing alleles for the complement of this set
    
    void getComplementCounts(const std::vector<string>& populationsToUse);
};


class MultiallelicCounts {
public:
    MultiallelicCounts(int nSamples, int starPosInput) : inbreedingCoefficient(0), chiSqPvalForInbreeding(1), bPhased(false) {
        individualsHets.assign(nSamples, -1);
        missingIndividualsDot.assign(nSamples, 0);
        missingIndividualsAny.assign(nSamples, 0);
        missingHaplotypesDot.assign(2*nSamples, 0);
        missingHaplotypesStar.assign(2*nSamples, 0);
        haplotypeVariants.assign(2*nSamples, -1);
        starPos = starPosInput;
    };
    
    int starPos;
    double inbreedingCoefficient;
    double chiSqPvalForInbreeding;
    bool bPhased;
    std::vector<int> individualsHets;
    std::vector<int> haplotypeVariants;
    
    std::vector<int> missingIndividualsDot;  // Individuals which are denoted in the VCF as missing data ./.
    std::vector<int> missingIndividualsAny;  // Individuals with missing data (as above) or which have either haplotype denoted in the VCF as * alleles (missing at the site because of an indel nearby)
    
    std::vector<int> missingHaplotypesDot;  // Haplotypes which are denoted in the VCF as missing data ./.
    std::vector<int> missingHaplotypesStar;   // Haplotypes which are denoted in the VCF as * alleles (missing at the site because of an indel nearby)
    
    void getMultiallelicCounts(const std::vector<string>& genotypes);
    double getPiThisVariant();
    double getHeterozygosityThisVariant();
    
};

#endif /* UtilsSetCounts_hpp */
