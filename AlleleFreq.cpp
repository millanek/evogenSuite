//
//  AlleleFreq.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "AlleleFreq.hpp"
//#include "process_vcf_annotation_tools.h"
#include <deque>

#define SUBPROGRAM "alleleFreq"

#define DEBUG 1

static const char *AF_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt\n"
"Calculate the Allele Frequencies per population/species from a VCF \n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -g, --use-genotype-probabilities        (optional) use genotype probabilities (GP tag) if present\n"
"                                               if GP not present calculate genotype probabilities from likelihoods (GL or PL tags) using a Hardy-Weinberg prior\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:";


static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "use-genotype-probabilities", no_argument, NULL, 'g'},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string runName = "out";
    static bool useGenotypeProbabilities = false;
}


int AFmain(int argc, char** argv) {
    parseAFoptions(argc, argv);
    
    SetInformation setInfo(opt::setsFile);
    
    int totalVariantNumber = 0;
    std::vector<std::string> fields;
    int reportProgressEvery = 10000; string chr; string coord; double coordDouble;
    std::clock_t startTime = std::clock();
    std::ostream* outFileAF = createWriter(stripExtension(opt::setsFile) + "_" + opt::runName + "_AF" + ".txt");
    
    string line; // for reading the input files
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
            printAFheader(setInfo, outFileAF);
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1]; coordDouble = stringToDouble(coord);
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size(), chr, coord);
            c->getSetVariantCountsSimple(genotypes, setInfo.posToPopMap);
            // std::cerr << "Here:" << totalVariantNumber << std::endl;

            
            *outFileAF << chr << "\t" << coord << "\t" << refAllele << "\t" << altAllele;
            if (opt::useGenotypeProbabilities) {
                int likelihoodsOrProbabilitiesTagPosition = c->checkForGenotypeLikelihoodsOrProbabilities(fields);
                if (likelihoodsOrProbabilitiesTagPosition == LikelihoodsProbabilitiesAbsent) {
                    printMissingLikelihoodsWarning(fields[0], fields[1]);
                    opt::useGenotypeProbabilities = false;
                } else c->getAFsFromGenotypeLikelihoodsOrProbabilities(genotypes,setInfo.posToPopMap,likelihoodsOrProbabilitiesTagPosition);
                for(std::map<string,double>::iterator iter =  c->setAAFsFromLikelihoods.begin(); iter != c->setAAFsFromLikelihoods.end(); ++iter) {
                    *outFileAF << "\t" << iter->second;
                }
            
            } else {
                for(std::map<string,double>::iterator iter =  c->setAAFs.begin(); iter != c->setAAFs.end(); ++iter) {
                    *outFileAF << "\t" << iter->second;
                }
            }
            *outFileAF << "\n";

            genotypes.clear(); genotypes.shrink_to_fit();
            
            delete c;
        }
    }
    
    outFileAF->flush();
    return 0;
}



void parseAFoptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'g': opt::useGenotypeProbabilities = true; break;
            case 'h':
                std::cout << AF_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << AF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}


void printAFheader(const SetInformation& setInfo, std::ostream* outFileAF) {
    *outFileAF << "chr" << "\t" << "coord" << "\t" << "ref" << "\t" << "alt";
    for(std::map<string,std::vector<size_t>>::const_iterator iter =  setInfo.popToPosMap.begin(); iter != setInfo.popToPosMap.end(); ++iter) {
        *outFileAF << "\t" << iter->first;
    }
    *outFileAF << "\n";
}
