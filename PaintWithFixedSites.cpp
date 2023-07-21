//
//  PainWithFixedSites.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 14.06.23.
//

#include "PaintWithFixedSites.hpp"

#define SUBPROGRAM "PaintFixed"
#define MIN_SETS 2

#define DEBUG 1

static const char *DISTOUT_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt Hybrids.txt SOURCE_POP1,SOURCE_POP2\n"
"Find the genotype of the individuals listed in Hybrids.txt at sites fixed between SOURCE_POP1 and SOURCE_POP2:\n"
"The POPULATIONS.txt file; This file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The Hybrids.txt should contain IDs of the individuals whose genomes are to be painted - one per line\n"
"SOURCE_POP1,SOURCE_POP2 should specify two populations used for defining the fixed sites\n"
"All output will be in one file: .txt\n"
"\n"
HelpOption RunNameOption MaxMissOption
//"       -i, --allow-indels-and-multiallelics    (optional) also calculate the PBS score for indels, and multiallelics\n"
//"                                               for multiallelics, the first alternate allele is considered\n"
//"                                               sites where the ALT allele is simply '*' are ignored\n"
// "       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
//"       --no-scaling                            (optional) do not scale distances by the proportion of missing data in the ingroup (default: scale)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:m:";

enum { OPT_ACC_GEN_BED, OPT_AF  };

static const struct option longopts[] = {
    { "fixedW",   required_argument, NULL, 'f' },
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "maxMissing", required_argument, NULL, 'm' },
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string hybridsFile;
    static std::vector<string> sourcePops;
    static string runName = "";
    static double maxMissing = 0.2;
}

int PaintFixedMain(int argc, char** argv) {
    parsePaintFixedOptions(argc, argv);
    
    SetInformation setInfo(opt::setsFile, MIN_SETS);
    
    string line; // for reading the input files
    
    // 1) Dividing up and assigning the populations
    std::vector<string> hybrids; // Will hold the names of the hybrid individuals
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* hybridsFile = new std::ifstream(opt::hybridsFile.c_str());
    std::ofstream* outFile = new std::ofstream("PaintFixed" + opt::runName + ".txt");
    *outFile << "chr\tpos";
    
    while (getline(*hybridsFile,line)) {
        hybrids.push_back(line);
        *outFile << "\t" << line;
    }  *outFile << std::endl;
    
    /*
    // 3) Praparing the containers to hold the results
    // 3a) for per-SNP results
    std::vector<std::vector<double>> initVectorFixed(ingroups.size());
    std::vector<std::vector<std::vector<double>>> DxyFixedWindowPerSNP(outgroups.size(),initVectorFixed);
    // 3b) for calculated Dxy values
    std::vector<int> ingroupsMissingInit(ingroups.size(),0);
    std::vector<std::vector<int>> missingDist(outgroups.size(),ingroupsMissingInit);
    
    // 3b) for calculated Dxy values
    std::vector<double> ingroupsResultInit(ingroups.size(),0.0);
    std::vector<std::vector<double>> DxyFixedWindowAveraged(outgroups.size(),ingroupsResultInit);
    */
    
    // 4) Going through the VCF
    // 4a) Preparing useful variables
    int totalVariantNumber = 0;
    
    std::vector<std::string> fields; std::vector<std::string> sampleNames;
    int reportProgressEvery = 1000; std::clock_t startTime = std::clock();
    // 4b) Looping through the file
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            sampleNames.assign(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
    
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            c->getSetVariantCounts(genotypes, setInfo.posToPopMap, v);
            c->fillMissingFistAlleleInfo(setInfo);
            genotypes.clear(); genotypes.shrink_to_fit();
            
            double source_1_AAF = c->setAAFs.at(opt::sourcePops[0])[0];
            double source_2_AAF = c->setAAFs.at(opt::sourcePops[1])[0];
            double source_1_missingness = c->setFirstAlleleMissingness.at(opt::sourcePops[0]);
            double source_2_missingness = c->setFirstAlleleMissingness.at(opt::sourcePops[1]);

            if (source_1_missingness < opt::maxMissing && source_2_missingness < opt::maxMissing) {
                if ((source_1_AAF == 0 && source_2_AAF == 1)
                    || (source_1_AAF == 1 && source_2_AAF == 0)) {
                    *outFile << v.chr << "\t" << v.posInt;
                    for (int i = 0; i != hybrids.size(); i++) {
                        double hybridAAF = c->setAAFs.at(hybrids[i])[0];
                        
                        if (hybridAAF == -1) { *outFile << "\tNA"; }
                        else {
                            if (source_1_AAF == 0) *outFile << "\t" << hybridAAF;
                            else *outFile << "\t" << 1 - hybridAAF;
                        }
                    }
                    *outFile << std::endl;
                }
            }
            delete c;
        }
    }
    
    return 0;
    
}



void parsePaintFixedOptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'm': arg >> opt::maxMissing; break;
            case 'h':
                std::cout << DISTOUT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 4) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 4)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DISTOUT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::hybridsFile = argv[optind++];
    string sourcePopsString = argv[optind++];
    opt::sourcePops = split(sourcePopsString,',');
    assert(opt::sourcePops.size() == 2);
}

