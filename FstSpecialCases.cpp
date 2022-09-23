//
//  FstSpecialCases.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 23.09.22.
//

#include "FstSpecialCases.hpp"



std::cerr << "Calculating Fst using variants from: " << opt::vcfFile << std::endl;
std::cerr << "using a sliding window of size: " << opt::windowSize << " variants and sliding in steps of: " << opt::windowStep << std::endl;

// Load up the annotation file if provided
Annotation wgAnnotation(opt::annotFile, false);
if (!opt::annotFile.empty()) {
    string snpCategoryFstFileName = opt::runName + "SNPcategory_fst.txt";
    snpCategoryFstFile = new std::ofstream(snpCategoryFstFileName.c_str());
    *snpCategoryFstFile << "SNPcategory" << "\t" << "thisSNPFst" << "\t" << "thisSNPDxy" << "\t" << "scaffold" << "\t" << "position" << std::endl;
}

if (!opt::ancSets.empty()) {
    ancSetsFile = new std::ifstream(opt::ancSets);
    string ancOutFileName = fileRoot + opt::runName + "ancestralSNPs_fst.txt";
    ancSetsOutFile = new std::ofstream(ancOutFileName);
    *ancSetsOutFile << "scaffold" << "\t" << "position" << "\t" << "AncAllelePopulation" << "\t" << "Fst" << "\t" << "ancSet1_segregating" << "\t" << "ancSet2_segregating" << std::endl;
    string ancSet1String; string ancSet2String;
    getline(*ancSetsFile, ancSet1String);
    getline(*ancSetsFile, ancSet2String);
    ancSet1 = split(ancSet1String, ','); ancSet2 = split(ancSet2String, ',');
    std::sort(ancSet1.begin(),ancSet1.end()); std::sort(ancSet2.begin(),ancSet2.end());
}

if (opt::regAbove > 0) {
    string regionsAboveFstFileName = fileRoot + opt::runName + "_w_" + numToString(opt::windowSize) + "_fst_above" + numToString(opt::regAbove) + ".txt";
    regionsAboveFstFile = new std::ofstream(regionsAboveFstFileName.c_str());
}

string fileRoot = "";
std::ofstream* snpCategoryFstFile;
std::ofstream* regionsAboveFstFile; bool inRegAbove = false;
std::ifstream* ancSetsFile; std::ofstream* ancSetsOutFile;
std::vector<string> ancSet1; std::vector<string> ancSet2;

if (opt::regAbove > 0) *regionsAboveFstFile << "scaffold" << "\t" << "Start" << "\t" << "End" << std::endl;

if (!opt::ancSets.empty()) {
    ancSet1Loci = locateSet(sampleNames, ancSet1);
    ancSet2Loci = locateSet(sampleNames, ancSet2);
    std::cerr << "Ancestral Set1 loci: " << std::endl;
    print_vector(ancSet1Loci, std::cerr);
    std::cerr << "Ancestral Set2 loci: " << std::endl;
    print_vector(ancSet2Loci, std::cerr);
    n1anc = ancSet1Loci.size() * 2; n2anc = ancSet2Loci.size() * 2;
}


if (!opt::annotFile.empty()) {
/*    string SNPcategory = wgAnnotation.getCategoryOfSNP(scaffold, loc);
    double thisSNPFst = FstNumerator/FstDenominator;
    *snpCategoryFstFile << SNPcategory << "\t" << thisSNPFst << "\t" << thisSNPDxy << "\t" << scaffold << "\t" << loc << std::endl;
 */
}


//    countedVariantNumber++;
    
   /* if (!opt::ancSets.empty()) {
        double thisSNPFst = FstNumerator/FstDenominator;
        if (thisSNPFst < 0) { thisSNPFst = 0; }
        string AA = split(info[info.size()-1],'=')[1];
        //std::cerr << "AA=" << " " << AA << std::endl;
        FourSetCounts c;
        if (AA == fields[3]) {
            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"ref");
            *ancSetsOutFile << scaffold << "\t" << fields[1] << "\t" << c.set1daAF-c.set2daAF << "\t" << thisSNPFst << "\t";
            if (c.set3daAF > 0 & c.set3daAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
            if (c.set4daAF > 0 & c.set4daAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
        } else if (AA == fields[4]) {
            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"alt");
            *ancSetsOutFile << scaffold << "\t" << fields[1] << "\t" << c.set1daAF-c.set2daAF << "\t" << thisSNPFst << "\t";
            if (c.set3daAF > 0 & c.set3daAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
            if (c.set4daAF > 0 & c.set4daAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
            // std::cerr << "AA=alt" << " " << c.set1daAF << " " << c.set2daAF << std::endl;
        } else {
            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"N");
            *ancSetsOutFile << scaffold << "\t" << fields[1] << "\t" << "-888" << "\t" << thisSNPFst << "\t";
            if (c.set3AltAF > 0 & c.set3AltAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
            if (c.set4AltAF > 0 & c.set4AltAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
        }
    }
    
    if (opt::windowSize == 1) {
        double Fst = FstNumerator/FstDenominator;
        if (Fst < 0) Fst = 0;
        *pFst << countedVariantNumber << "\t" << scaffold + "\t" + fields[1] << "\t" << Fst << "\t" << thisSNPDxy << std::endl;
        
    } else if ((opt::windowSize > 0) && (countedVariantNumber % opt::windowStep == 0) && countedVariantNumber >= opt::windowSize) {
        std::vector<double> windowFstNumerators(fstNumerators.end()-opt::windowSize, fstNumerators.end());
        std::vector<double> windowFstDenominators(fstDenominators.end()-opt::windowSize, fstDenominators.end());
        double windowFst = calculateFst(windowFstNumerators, windowFstDenominators); if (windowFst < 0) windowFst = 0;
        std::vector<double> windowDxyVec(DxyVector.end()-opt::windowSize, DxyVector.end());
        double windowDxy = vector_average(windowDxyVec);
        if (opt::windowSize == opt::windowStep) {
            std::vector<string> s = split(windowStartEnd, '\t');
            if (s[0] == scaffold) {
                windowStartEnd = windowStartEnd + "\t" + fields[1];
                windowEnd = atoi(fields[1].c_str());
                double windowDxyIncNonSeg = vector_average_withRegion(windowDxyVec, windowEnd-windowStart);
                *pFst << countedVariantNumber-opt::windowSize+1 << "\t" << windowStartEnd << "\t" << windowFst << "\t" << windowDxy << "\t" << windowDxyIncNonSeg << "\t" << windowFstDenominators.size() << std::endl;
                if (opt::regAbove > 0) {
                    if (windowFst >= opt::regAbove && !inRegAbove) {
                        inRegAbove = true;
                        *regionsAboveFstFile << s[0] << "\t" << s[1] << "\t";
                    } else if (windowFst < opt::regAbove && inRegAbove) {
                        inRegAbove = false;
                        *regionsAboveFstFile << s[1] << std::endl;
                    }
                }
            }
        } else {
            *pFst << countedVariantNumber-opt::windowSize+1 << "\t" << windowMiddleVariant << "\t" << windowFst << "\t" << windowDxy << "\t" << windowFstDenominators.size() << std::endl;
        }
       
        
       
    }
    */


void getVariantCountsForFst(const std::vector<std::string>& fields, SetCounts* thisVariantCounts, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci) {

    int numSamples = (int)fields.size()-NUM_NON_GENOTYPE_COLUMNS;
    thisVariantCounts->individualsWithVariant.assign(numSamples,0);
    thisVariantCounts->missingGenotypesPerIndividual.assign(numSamples,false);
    thisVariantCounts->haplotypesWithVariant.assign(numSamples*2,0);
    
    thisVariantCounts->set1individualsWithVariant.assign(set1_loci.size(),0);
    thisVariantCounts->set2individualsWithVariant.assign(set2_loci.size(),0);
    int n1 = (int)(set1_loci.size()*2); int n2 = (int)(set2_loci.size()*2);
    thisVariantCounts->set1HaplotypeVariant.assign(n1,0);
    thisVariantCounts->set2HaplotypeVariant.assign(n2,0);
    thisVariantCounts->set1_n_withoutMissing = n1;
    thisVariantCounts->set2_n_withoutMissing = n2;
    
    std::vector<std::string> altAlleles = split(fields[4], ',');
    thisVariantCounts->n_alt_alleles = (int)altAlleles.size();
    int alleleAsMissing = -1;
    for (int i = 0; i < (int)altAlleles.size(); i++) {
        if (altAlleles[i] == "*") {
            thisVariantCounts->n_alt_alleles = thisVariantCounts->n_alt_alleles - 1;
            alleleAsMissing = i + 1;
        }
        if (altAlleles[i].length() > 1) {
            thisVariantCounts->bIndel = true;
        }
    }
    
    if (thisVariantCounts->n_alt_alleles == 1) {
        int alt = 1;

        if (alleleAsMissing == 1) alt = 2;

        int set1i = 0; int set2i = 0; int set1hapI = 0; int set2hapI = 0;
        std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
        //std::cerr << fields[0] << "\t" << fields[1] << "\tgenotypes.size()" << genotypes.size() << std::endl;
        for (const size_t i : set1_loci) {
            int v1int = genotypes[i][0] - '0'; int v2int = genotypes[i][2] - '0';
            if (v1int == alt) {
                thisVariantCounts->set1Count++; thisVariantCounts->set1individualsWithVariant[set1i]++;
                thisVariantCounts->set1HaplotypeVariant[set1hapI]++;
            } else if (genotypes[i][0] == '.' || v1int == alleleAsMissing) {
               //std::cerr << fields[0] << "\t" << fields[1] << "\t" << genotypes[i][0] << "\t" << i << std::endl;
                thisVariantCounts->missingGenotypesPerIndividual[i] = true;
                thisVariantCounts->set1_n_withoutMissing = thisVariantCounts->set1_n_withoutMissing - 1;
            }
            if (v2int == alt) {
                thisVariantCounts->set1Count++; thisVariantCounts->set1individualsWithVariant[set1i]++;
                thisVariantCounts->set1HaplotypeVariant[set1hapI+1]++;
            } else if (genotypes[i][2] == '.' || v2int == alleleAsMissing) {
                thisVariantCounts->missingGenotypesPerIndividual[i] = true;
                //std::cerr << fields[0] << "\t" << fields[1] << "\t" << genotypes[i][2] << std::endl;
                thisVariantCounts->set1_n_withoutMissing = thisVariantCounts->set1_n_withoutMissing - 1;
            }
            set1i++; set1hapI = set1hapI+2;
        }
        for (const size_t i : set2_loci) {
            int v1int = genotypes[i][0] - '0'; int v2int = genotypes[i][2] - '0';
            if (genotypes[i][0] - '0' == alt) {
                thisVariantCounts->set2Count++; thisVariantCounts->set2individualsWithVariant[set2i]++;
                thisVariantCounts->set2HaplotypeVariant[set2hapI]++;
            } else if (genotypes[i][0] == '.' || v1int == alleleAsMissing) {
                thisVariantCounts->missingGenotypesPerIndividual[i] = true;
                thisVariantCounts->set2_n_withoutMissing = thisVariantCounts->set2_n_withoutMissing - 1;
            }
            if (genotypes[i][2] - '0' == alt) {
                thisVariantCounts->set2Count++; thisVariantCounts->set2individualsWithVariant[set2i]++;
                thisVariantCounts->set2HaplotypeVariant[set2hapI+1]++;
            } else if (genotypes[i][2] == '.' || v2int == alleleAsMissing) {
                thisVariantCounts->missingGenotypesPerIndividual[i] = true;
                thisVariantCounts->set2_n_withoutMissing = thisVariantCounts->set2_n_withoutMissing - 1;
            }
            set2i++; set2hapI = set2hapI+2;
        }
        for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
            if (genotypes[i][0] - '0' == alt)
                thisVariantCounts->overall++;
            if (genotypes[i][2] - '0' == alt)
                thisVariantCounts->overall++;
        }
    }
}
