//
//  DistanceMatrix.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 21.09.22.
//

#include "DistanceMatrix.hpp"

#define SUBPROGRAM "GlobalPairStats"
#define MIN_SETS 2

static const char *STATS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE POPULATIONS.txt\n"
"Calculating various global (e.g. genome-wide) statistics from a VCF file\n"
"By default calculates a distance matrix that can be used e.g. for a Neighbour-Joining tree\n"
"\n"
HelpOption RunNameOption MaxMissOption
"\n"
"       You can add an extra option specifying the number of accessible bp; the diff_matrix will be defined in d_XY\n"
"       --numAccessibleBP=NUM                   (optional) Number of accessible bp in the region\n"
"                                               the diff_matrix will be defined in d_XY\n"
"       --bootstrap=NUM_REP,BLOCKSIZE           (default 100,1000) Generate NUM_REP distance matrices resampling with replacement in blocks of BLOCKSIZE variants\n"
"\n"
"       Additional statistics that can optionally be calculated are:\n"
"       --doubletons                            Calculates the distribution of doubletons between/within sets\n"
"       --private-variants                      Calculates the numbers of private variants in each of the populations\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_NUM_ACCESSIBLE, OPT_DOUBLETON, OPT_SHARED_POL, OPT_BOOTSTRAP, OPT_PRIVATE_VARS };

static const char* shortopts = "hn:";



static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "maxMissing", required_argument, NULL, 'm' },
    { "numAccessibleBP", required_argument, NULL, OPT_NUM_ACCESSIBLE },
    { "bootstrap", required_argument,    NULL, OPT_BOOTSTRAP },
    { "doubletons", no_argument,    NULL, OPT_DOUBLETON },
    { "private-variants", no_argument,    NULL, OPT_PRIVATE_VARS },
    { "shared-polymophisms", no_argument, NULL, OPT_SHARED_POL },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    // Individual or population files required for outputtng statistics
    static string setsFile = "";

    // Boolean flags indicating which statistics should be calculated
    static bool bSharedPolymorph = false;
    static bool countPrivateVars = false;
    static bool bDoubleton = false;
    static int numAccessibleBP = -1;
    
    static string runName = "";
    static double maxMissing = 0.2;

    static int bootstrapBlockSize;
    static int bootstrapReps;
    static bool bDoBootstrap = false;
}

int globalStatsMain(int argc, char** argv) {
    parseStatsOptions(argc, argv);
    string fileRoot = stripExtension(opt::vcfFile);
    
    SetInformation setInfo(opt::setsFile, MIN_SETS);
    
    DistanceMatrices d(setInfo, fileRoot, opt::runName, opt::bDoBootstrap, opt::bDoubleton);
    
    std::cerr << "Calculating statistics from: " << opt::vcfFile << std::endl;
    if (opt::bSharedPolymorph) {
        std::cerr << "Counts of hets (per individual) shared with others are going to be output to: " << fileRoot + ".sharedHets.txt" << std::endl;
    }
    if (opt::bDoubleton) {
        std::cerr << "The distribution of doubletons is going to be output to: " <<  fileRoot + ".doubletons.txt" << std::endl;
    }
    if (opt::countPrivateVars) {
        std::cerr << "The numbers of private variants fixed in groups defined in " << opt::setsFile << " are going to be output to: " <<  fileRoot + "_" + stripExtension(opt::setsFile) + ".privateFixedVars.txt" << std::endl;
    }
    
    // Start reading from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    string line; std::vector<string> fields;
    std::clock_t startTime = std::clock();
    int totalVariantNumber = 0; int reportProgressEvery = 1000;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            // Read the sample names, initialise output variables, and (if supplied) locate populations
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
            
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            fields = split(line, '\t');
            
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            c->getSetVariantCounts(genotypes, setInfo.posToPopMap, v);
            //c->calculatePiPerVariantPerSet();
            genotypes.clear(); genotypes.shrink_to_fit();
            
            d.addAllPairwiseDistances(setInfo, c, opt::maxMissing);
            
            if (opt::bDoBootstrap && totalVariantNumber % opt::bootstrapBlockSize == 0) d.addBootstrapBlock();
            
            if (bFirstTwoAllelesInformative(setInfo,c,opt::maxMissing)) {
                if (opt::bDoubleton) d.doubleton_analysis(setInfo, c, opt::maxMissing);
                if (opt::countPrivateVars) d.privateVars_analysis(setInfo, c, opt::maxMissing);
                if (opt::bSharedPolymorph) d.sharedPolymorphism_analysis(setInfo, c, opt::maxMissing);
            }
            
            delete c;
        }
    }
    
    // Finished processing the VCF. Now summarising + normalising the results:
    normaliseMatrixBasedOnMissigness(d.diffMatrix, d.pairwiseMissingness, totalVariantNumber);
    if (opt::numAccessibleBP > -1) scaleMatrixByAccessibleBP(d.diffMatrix, opt::numAccessibleBP);
    // Printing the results
    d.printDxyResults(setInfo);
    // Bootstrap if requested
    if (opt::bDoBootstrap) d.doBootstrapResampling(setInfo, opt::bootstrapReps, totalVariantNumber, opt::numAccessibleBP);
    
        
    // Printing doubletons
   // if (opt::bDoubleton)  print_doubleton_distribution(fileRoot, pop_unique, doubletons);
    
    /*
    // Printing het counts
    if (opt::countHets)
        print_het_counts(fileRoot, sampleNames, hetCounts, hetsSharedWithOthers);
    
    // Printing statistics for private variants
    if (opt::countPrivateVars)
        print_privateFixedVarsSummary(fileRoot, populationsToUse, opt::populationsFile, privateVarCounts);
    
     */
    
    return 0;
}

void parseStatsOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> bootstrapRepsBlocksize;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'm': arg >> opt::maxMissing; break;
            case OPT_DOUBLETON: opt::bDoubleton = true; break;
            case OPT_SHARED_POL: opt::bSharedPolymorph = true; break;
            case OPT_BOOTSTRAP:
                bootstrapRepsBlocksize = split(arg.str(), ',');
                opt::bootstrapBlockSize = atoi(bootstrapRepsBlocksize[0].c_str());
                opt::bootstrapReps = atoi(bootstrapRepsBlocksize[1].c_str());
                opt::bDoBootstrap = true;
                break;
            case OPT_PRIVATE_VARS: opt::countPrivateVars = true; break;
            case OPT_NUM_ACCESSIBLE: arg >> opt::numAccessibleBP; break;
            case 'h':
                std::cout << STATS_USAGE_MESSAGE;
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
    
    if (bootstrapRepsBlocksize.size() != 2) {
        std::cerr << "Error in the --bootstrap option: It needs two numbers separated by a comma; e.g. '--bootstrap 100,1000\n";
        die = true;
    }
        
    if (die) {
        std::cout << "\n" << STATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}


void normaliseMatrixBasedOnMissigness(std::vector<std::vector<double>>& dxyMatrix, const std::vector<std::vector<int>>& missignessMatrix, const int totalSites) {
    for (int i = 0; i < (int)dxyMatrix.size(); i++) {
        for (int j = 0; j < (int)dxyMatrix.size(); j++) {
            double proportionUsed = 1 - ((double)missignessMatrix[i][j]/totalSites);
            dxyMatrix[i][j] = dxyMatrix[i][j]/proportionUsed;
        }
    }
}

void scaleMatrixByAccessibleBP(std::vector<std::vector<double>>& dxyMatrix, const int accessibleBP) {
    for (int i = 0; i < dxyMatrix.size(); i++) {
        for (int j = 0; j < dxyMatrix.size(); j++) {
            dxyMatrix[i][j] =  (double)dxyMatrix[i][j]/accessibleBP;
        }
    }
}

bool bFirstTwoAllelesInformative(const SetInformation& setInfo, const GeneralSetCounts* c, const double maxMissing) {
    // Make sure that ref and first alt allele make at least (1 - maxMissing) proportion of the dataset
    if (c->getOverallCountOfRefAllele() + c->getOverallCountOfFirstAltAllele() > (ploidy * setInfo.IDsToPopMap.size() * (1 - maxMissing) ))
        return true;
    else
        return false;
}

