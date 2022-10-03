//
//  UtilsSetCounts.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 29.09.22.
//

#include "UtilsSetCounts.hpp"


// Works only on biallelic markers
void GeneralSetCounts::getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap, const VariantInfo& v) {
    
    
    
    // Go through the genotypes and record all SNP alleles
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
        if (genotypes[i][0] == '.' || genotypes[i][2] == '.') continue;
        std::string species = posToSpeciesMap.at(i);
        int firstAllele = genotypes[i][0] - '0'; int secondAllele = genotypes[i][2] - '0';
        bool firstAlleleSNP = false; bool secondAlleleSNP = false;
        if (std::find(v.SNPAlleleIndices.begin(), v.SNPAlleleIndices.end(), firstAllele) != v.SNPAlleleIndices.end()) {
            setAlleles.at(species).push_back(firstAllele);
            firstAlleleSNP = true;
        }
        if (std::find(v.SNPAlleleIndices.begin(), v.SNPAlleleIndices.end(), secondAllele) != v.SNPAlleleIndices.end()) {
            setAlleles.at(species).push_back(secondAllele);
            secondAlleleSNP = true;
        }
        
        // This is for calculating heterozygosity
        if (firstAlleleSNP && secondAlleleSNP) {
            setSNPPairCounts.at(species)++;
            if (firstAllele != secondAllele) setHetCounts.at(species)++;
        }
    }
    //std::cerr << "Went through genotypes" << std::endl;
    
    
    // it->second - A vector with all alleles for this species
    for(std::map<string,std::vector<int>>::iterator it = setAlleles.begin(); it != setAlleles.end(); ++it) {
        std::string species = it->first;
        if(it->second.size() == 0) {        // Missing data for this set
            setAAFs.at(species).push_back(-1);
            setRefCounts.at(species) = 0;
            setAltAlleleCounts.at(species).push_back(0);
        } else {
            setRefCounts.at(species) = (int)std::count(it->second.begin(), it->second.end(), 0);
            for (int i = 1; i < v.SNPAlleleIndices.size(); i++) {
                int thisAAcount = (int)std::count(it->second.begin(), it->second.end(), v.SNPAlleleIndices[i]);
                setAltAlleleCounts.at(species).push_back(thisAAcount);
                double thisAAF = (double)thisAAcount/it->second.size();
                setAAFs.at(species).push_back(thisAAF);
           //     std::cerr << "v.SNPaltAlleleIndices[i]" << v.SNPAlleleIndices[i] << std::endl;
           //     std::cerr << "thisAAcount" << thisAAcount << std::endl;
           //     std::cerr << "it->second.size()" << it->second.size() << std::endl;
            }
        }
    }
}

/*// Works only on biallelic markers
void GeneralSetCounts::getAlternateAlleleFreqFromOutgroup(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    
    
    // If at least one of the outgroup individuals has non-missing data
    // Find out what is the "ancestral allele" - i.e. the one more common in the outgroup
    int AAint;
    try {
        if (setAlleleCounts.at("Outgroup") > 0) {
            if (setRefCounts.at("Outgroup") > setAltCounts.at("Outgroup")) { AAint = 0; }
            else { AAint = 1; }
        }
    } catch (std::out_of_range& e) { AAint = -1; }
    
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (setAlleleCounts.at(it->first) > 0) {
            if (AAint == 0) { // Ancestral allele seems to be the ref, so derived is alt
                setDAFs[it->first] = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
            } else if (AAint == 1) { // Ancestral allele seems to be alt, so derived is ref
                setDAFs[it->first] = (double)setRefCounts.at(it->first)/setAlleleCounts.at(it->first);
            }
        }
    }
} */

void GeneralSetCounts::calculatePiPerVariantPerSet() {
    for(std::map<string,std::vector<int>>::iterator it = setAlleles.begin(); it != setAlleles.end(); ++it) {
        string thisSet = it->first;
        std::vector<int> thisSetAlleles = it->second;
        int n = (int) thisSetAlleles.size();
        double pi = 0;
        if (n > 0) {
            double piSum = 0;
            for (std::vector<std::string>::size_type i = 0; i < n - 1; i++) {
                for (std::vector<std::string>::size_type j = i + 1; j < n; j++) {
                    if (thisSetAlleles[i] != thisSetAlleles[j]) piSum++;
                }
            }
            pi = ( 2.0 / (n*(n-1)) ) * piSum;
        }
        piPerVariantPerSet.at(thisSet) = pi;
    }
}

void GeneralSetCounts::calculateHeterozygosityPerVariantPerSet() {
    for(std::map<string,int>::iterator it = setSNPPairCounts.begin(); it != setSNPPairCounts.end(); ++it) {
        string thisSet = it->first;
        int thisSetSNPpairs = it->second;
        int thisSetHets = setHetCounts.at(thisSet);
        double hetThisVar = (double)thisSetHets/(double)thisSetSNPpairs;
        
        hetPerVariantPerSet.at(thisSet) = hetThisVar;
    }
}

std::vector<double> GeneralSetCountsWithLikelihoods::probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species) {
    std::vector<double> thisProbabilities; thisProbabilities.assign(3, 0.0);
    double multiple0 = thisLikelihoods[0]*setHWEpriorsFromAAFfromGT[species][0];
    double multiple1 = thisLikelihoods[1]*setHWEpriorsFromAAFfromGT[species][1];
    double multiple2 = thisLikelihoods[2]*setHWEpriorsFromAAFfromGT[species][2];
    double sum = multiple0 + multiple1 + multiple2;
    
    thisProbabilities[0] = multiple0/sum;
    thisProbabilities[1] = multiple1/sum;
    thisProbabilities[2] = multiple2/sum;
    
    return thisProbabilities;
}


void GeneralSetCountsWithLikelihoods::setHWEpriorsFromAFfromGT() {
    double AF;
    // Alternative allele frequencies
    for(std::map<string,std::vector<double>>::iterator it = setAAFs.begin(); it != setAAFs.end(); ++it) {
        if (it->second[0] >= 0) AF = it->second[0]; else AF = averageFirstAAF; // This should be average of AFs across populations where it is known
        setHWEpriorsFromAAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromAAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromAAFfromGT[it->first][2] = pow(AF,2);
    }
 /*   // Derived allele frequencies
    for(std::map<string,double>::iterator it = setDAFs.begin(); it != setDAFs.end(); ++it) {
        if (it->second >= 0) AF = it->second; else AF = averageDAF; // This should be average of AFs across populations
        setHWEpriorsFromDAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromDAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromDAFfromGT[it->first][2] = pow(AF,2);
    } */
}

void GeneralSetCountsWithLikelihoods::getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition) {
    if (likelihoodsProbabilitiesType == LikelihoodsProbabilitiesPL || likelihoodsProbabilitiesType == LikelihoodsProbabilitiesGL) {
        setHWEpriorsFromAFfromGT();
    }
    
    for (std::vector<std::string>::size_type i = 0; i < genotypeFields.size(); i++) {
       // std::cerr << "Here posToSpeciesMap.at(i): " << posToSpeciesMap.at(i) << std::endl;
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            // std::cerr << "Here the species is not in the map" << std::endl;
            continue;
        }
       // std::cerr << genotypeFields[i] << std::endl;
        std::string thisLikelihoodsOrProbabilitiesString = split(genotypeFields[i], ':')[likelihoodsOrProbabilitiesTagPosition];
        if (thisLikelihoodsOrProbabilitiesString == "." || thisLikelihoodsOrProbabilitiesString == "0,0,0") continue;
        
        else {
            setAlleleProbCounts.at(species) += 2;
            std::vector<double> thisLikelihoodsOrProbabilities = splitToDouble(thisLikelihoodsOrProbabilitiesString,',');
            std::vector<double> thisProbabilities;
            switch (likelihoodsProbabilitiesType)
            {
                case LikelihoodsProbabilitiesPL:
                    transformFromPhred(thisLikelihoodsOrProbabilities);
                   // print_vector(thisLikelihoodsOrProbabilities, std::cerr);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGL: transformFromGL(thisLikelihoodsOrProbabilities);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGP:
                    thisProbabilities = thisLikelihoodsOrProbabilities;
                    break;
            }
            if (setAAFsFromLikelihoods.at(species) == -1) setAAFsFromLikelihoods.at(species) = 0;
            setAAFsFromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
        }
    }
    
  /*  for(std::map<string,double>::iterator it = setAAFsFromLikelihoods.begin(); it != setAAFsFromLikelihoods.end(); ++it) {
        if (setAAFsFromLikelihoods.at(it->first) != -1) {
            double AF = it->second/setAlleleProbCounts.at(it->first);
            it->second = AF;
            if (AAint == AncestralAlleleRef) {
                setDAFsFromLikelihoods.at(it->first) = AF;
            } else if (AAint == AncestralAlleleAlt) {
                setDAFsFromLikelihoods.at(it->first) = (1 - AF);
            }
        }
    }
     */
}

int GeneralSetCountsWithLikelihoods::returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag) {
    // Find the position of GQ (genotype quality) in the genotypeData vector below
    std::vector<std::string>::iterator TAGit; int TAGi = std::numeric_limits<int>::min();
    TAGit = find (format.begin(), format.end(), tag);
    if (TAGit == format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample GQ info" << std::endl;
    } else {
        TAGi = (int)std::distance( format.begin(), TAGit );
        //hasGQ = true;
    }
    return TAGi;
}

int GeneralSetCountsWithLikelihoods::checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields) {
    std::vector<std::string> format = split(vcfLineFields[8], ':');
    if (format.size() == 1) return LikelihoodsProbabilitiesAbsent; // The GT tag must be present in the first place
    
    int likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GP");
    if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGP; }
    else {
        likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GL");
        if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGL; }
        else {
            likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "PL");
            if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesPL; }
        }
    }
    return likelihoodsOrProbabilitiesTagPosition;
}

void MultiallelicCounts::getMultiallelicCounts(const std::vector<std::string>& genotypes) {
    
  /*  std::cerr << "genotypes.size(): " << genotypes.size() << std::endl;
    std::cerr << "missingIndividualsDot.size(): " << missingIndividualsDot.size() << std::endl;
    std::cerr << "missingIndividualsAny.size(): " << missingIndividualsAny.size() << std::endl;
    std::cerr << "missingHaplotypesDot.size(): " << missingHaplotypesDot.size() << std::endl;
    std::cerr << "starpos: " << starPos << std::endl;
   */
    
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
     //   std::cerr << "genotypes[i]: " << genotypes[i] << std::endl;
        
        
        if (genotypes[i][0] == '.' || genotypes[i][2] == '.') { missingIndividualsDot[i] = 1; missingIndividualsAny[i] = 1; }
     //   std::cerr << "Missing individuals assigned: " << std::endl;
        if (genotypes[i][0] == starPos || genotypes[i][2] == starPos) missingIndividualsAny[i] = 1;
          
        if (genotypes[i][0] == '.') { missingHaplotypesDot[2*i] = 1;}
        else if (genotypes[i][0] == starPos) { missingHaplotypesStar[2*i] = 1; }
        else haplotypeVariants[2*i] = genotypes[i][0];
        
        if (genotypes[i][2] == '.')  { missingHaplotypesDot[(2*i)+1] = 1;}
        else if (genotypes[i][2] == starPos) { missingHaplotypesStar[(2*i)+1] = 1; }
        else haplotypeVariants[(2*i)+1] = genotypes[i][2];
        
        if (missingIndividualsAny[i] == 0 && genotypes[i][0] != genotypes[i][2]) {
            individualsHets[i] = 1;
        }
        
        if (missingIndividualsAny[i] == 0 && genotypes[i][0] == genotypes[i][2]) {
            individualsHets[i] = 0;
        }
        
    }
}

void GeneralSetCountsWithComplements::getComplementCounts(const std::vector<string>& populationsToUse) {
    // Now fill in the allele frequencies
    for(std::map<string,std::vector<int>>::iterator it = setAltAlleleCounts.begin(); it != setAltAlleleCounts.end(); ++it) {
        int complementAAFCount = 0; int complementTotalAlleleCount = 0;
        for (int i = 0; i < populationsToUse.size(); i++) {
            if (it->first != populationsToUse[i]) {
                complementAAFCount += setAltAlleleCounts.at(populationsToUse[i])[0];
                complementTotalAlleleCount += setRefCounts.at(populationsToUse[i]) + setAltAlleleCounts.at(populationsToUse[i])[0];
            }
        }
        setAAFsComplement[it->first] = (double)complementAAFCount/complementTotalAlleleCount;
        setAlleleCountsComplement[it->first] = complementTotalAlleleCount;
    }
}
