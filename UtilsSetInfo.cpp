//
//  UtilsSetInfo.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 29.09.22.
//

#include "UtilsSetInfo.hpp"


// ---------------- Functions for dealing with SETS/POPULATIONS files

std::vector<size_t> locateSet(const std::vector<std::string>& sample_names, const std::vector<std::string>& set, bool printWarnings) {
    std::vector<size_t> setLocs;
    for (std::vector<std::string>::size_type i = 0; i != set.size(); i++) {
        std::vector<std::string>::const_iterator it = std::find(sample_names.begin(), sample_names.end(), set[i]);
        if (it == sample_names.end() && printWarnings) {
            std::cerr << "WARNING: Did not find the sample: \"" << set[i] << "\"" << std::endl;
            print_vector(sample_names, std::cerr,',');
        } else {
            size_t loc = std::distance(sample_names.begin(), it);
            setLocs.push_back(loc);
        }
    }
    return setLocs;
} 

std::vector<size_t> complementIndices(const size_t fullVectorSize, const std::vector<size_t>& originalIndices) {
    std::vector<size_t> indexComplement;
    std::vector<bool> isInOriginal(fullVectorSize, false);
    for (const size_t j : originalIndices) {
        isInOriginal[j] = true;
    }
    indexComplement.reserve(fullVectorSize - originalIndices.size());
    for (size_t i = 0; i < fullVectorSize; i++) {
        if (!isInOriginal[i])
        indexComplement.push_back(i);
    }
    return indexComplement;
}

size_t locateOneSample(std::vector<std::string>& sample_names, const std::string toFind) {
    size_t pos = 0;
    std::vector<std::string>::iterator it = std::find(sample_names.begin(), sample_names.end(), toFind);
    if (it == sample_names.end()) {
        std::cerr << "Did not find the sample: " << toFind << std::endl;
    } else {
        pos = std::distance(sample_names.begin(), it);
    }
    return pos;
}

void SetInformation::linkSetsAndVCFpositions(const std::vector<std::string>& sampleNames, bool printWarnings) {
    // print_vector_stream(sampleNames, std::cerr);
    for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
        try { posToPopMap[i] = IDsToPopMap.at(sampleNames[i]); } catch (const std::out_of_range& oor) {
            if (printWarnings) std::cerr << "WARNING: The sample " << sampleNames[i] << " is in the VCF but not assigned in the SETS.txt file" << std::endl;
        }
    }
    // Iterate over all the keys in the map to find the samples in the VCF:
    // Give an error if no sample is found for a species:
    for(std::map<string, std::vector<string>>::const_iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
        string sp =  it->first;
        //std::cerr << "sp " << sp << std::endl;
        std::vector<string> IDs = it->second;
        std::vector<size_t> spPos = locateSet(sampleNames, IDs,printWarnings);
        if (spPos.empty() && printWarnings) {
            std::cerr << "WARNING: Did not find any samples in the VCF for \"" << sp << "\"" << std::endl;
        }
        popToPosMap[sp] = spPos;
    }
    
}
