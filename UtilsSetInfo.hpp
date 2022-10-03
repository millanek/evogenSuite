//
//  UtilsSetInfo.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 29.09.22.
//

#ifndef UtilsSetInfo_hpp
#define UtilsSetInfo_hpp

#include "UtilsGeneral.hpp"

// ---------------- Functions and a class for dealing with SETS/POPULATIONS files

// a) Find which fields in the VCF file corresponds to samples listed in a given set
std::vector<size_t> locateSet(const std::vector<std::string>& sample_names, const std::vector<std::string>& set);
// b) locating one sample
size_t locateOneSample(std::vector<std::string>& sample_names, const std::string toFind);

//
class SetInformation {
public:
    
    SetInformation(string setsFileName) {
        std::ifstream* setsFile = new std::ifstream(setsFileName.c_str());
        if (setsFile->fail()) {
            std::cerr << "ERROR: The file " << setsFileName << " could not be opened\n";
            exit(1);
        }
        string line;
        while (getline(*setsFile, line)) {
            // std::cerr << line << std::endl;
            std::vector<string> ID_Pop = split(line, '\t');
            popToIDsMap[ID_Pop[1]].push_back(ID_Pop[0]);
            IDsToPopMap[ID_Pop[0]] = ID_Pop[1];
            //std::cerr << ID_Species[1] << "\t" << ID_Species[0] << std::endl;
        }
        
        for(std::map<string,std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
            populations.push_back(it->first);
            // std::cerr << it->first << std::endl;
        } std::cerr << "There are " << populations.size() << " sets (populations/species) " << std::endl;
    };
    
    
    
    std::vector<string> populations;
    std::map<string, string> IDsToPopMap;
    std::map<string, std::vector<string>> popToIDsMap;
    std::map<string, std::vector<size_t>> popToPosMap;
    std::map<size_t, string> posToPopMap;

    void linkSetsAndVCFpositions(const std::vector<std::string>& sampleNames);

};

#endif /* UtilsSetInfo_hpp */
