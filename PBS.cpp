//
//  PBS.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "PBS.hpp"

#define SUBPROGRAM "PBS"
#define MIN_SETS 3

#define DEBUG 1

static const char *PBS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt PBS_trios.txt\n"
"Calculate the PBS statistic from:\n"
"Sequencing of 50 human exomes reveals adaptation to high altitude. Science 329, 75–78 (2010).\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The PBS_trios.txt should contain names of three populations for which the PBS will be calculated:\n"
"POP1   POP2    POP3\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_PBS_SIZE_STEP.txt\n"
"\n"
HelpOption RunNameOption MaxMissOption
"       -f, --fixedW sizeKb                     fixed window size (default: 10kb)\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --annot=ANNOTATION.gffExtract           (optional) gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs PBS per gene (only exons, with introns, and with 3kb upstream)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hw:n:f:m:";

enum { OPT_ANNOT  };

static const struct option longopts[] = {
    { "fixedW",   required_argument, NULL, 'f' },
    { "window",   required_argument, NULL, 'w' },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "maxMissing", required_argument, NULL, 'm' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile; static string setsFile; static string PBStriosFile;
    //
    static string annotFile;
    static string runName = "";
    static int physicalWindowSize = 10000;
    static int windowSize = 20;
    static int windowStep = 10;
    static double maxMissing = 0.2;
}


int PBSmain(int argc, char** argv) {
    parsePBSoptions(argc, argv);
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile, MIN_SETS);
    
    // Initialise the PBS trios
    PBStrios t(opt::PBStriosFile, opt::runName, opt::windowSize, opt::windowStep, opt::physicalWindowSize, !opt::annotFile.empty());
    
    // Load up the annotation file if provided
    Annotation wgAnnotation(opt::annotFile, false);
    
    int currentWindowStart = 0; int currentWindowEnd = currentWindowStart + opt::physicalWindowSize;
    
    std::vector<std::string> fields;
    std::clock_t startTime = std::clock();
    int totalVariantNumber = 0; int reportProgressEvery = 1000;
    
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
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            c->getSetVariantCounts(genotypes, setInfo.posToPopMap,v);
            genotypes.clear(); genotypes.shrink_to_fit();
            // std::cerr << "Here:" << totalVariantNumber << std::endl;
            
            // find if we are in a gene:
            if (wgAnnotation.initialised) wgAnnotation.getSNPgeneDetails(v.chr, v.posInt);
            
            // std::cerr << coord << "\t";
            // print_vector_stream(SNPgeneDetails, std::cerr);
            // Now calculate the PBS stats:
            for (int i = 0; i != t.getTrios().size(); i++) {
                string set1 = t.getTrios()[i][0]; string set2 = t.getTrios()[i][1]; string set3 = t.getTrios()[i][2];
                
                double p1 = c->setAAFs.at(set1)[0];
                double p2 = c->setAAFs.at(set2)[0];
                double p3 = c->setAAFs.at(set3)[0];
                
                int n1 = c->setRefCounts.at(set1) + c->setAltAlleleCounts.at(set1)[0];
                int n2 = c->setRefCounts.at(set2) + c->setAltAlleleCounts.at(set2)[0];
                int n3 = c->setRefCounts.at(set3) + c->setAltAlleleCounts.at(set3)[0];
                
                if (bTrioInformativeThisSNP(p1,p2,p3,n1,n2,n3,set1,set2,set3,setInfo,opt::maxMissing)) t.usedVars[i]++;
                else continue;
                
                std::vector<double> thisSNP_PBS = calculatePBSfromAFs(p1,p2,p3,n1,n2,n3); 
                
                t.addSNPresultsToWindows(i,thisSNP_PBS,v.posInt);
                
                if (wgAnnotation.currentGene != "") t.addSNPresultsToGene(i, thisSNP_PBS, wgAnnotation);
                
                // Check if we are at the step if the SNP window:
                if (t.usedVars[i] > opt::windowSize && (t.usedVars[i] % opt::windowStep == 0)) {
                    // std::cerr << PBSresults[i][0][0] << std::endl;
                    *t.outFiles[i] << v.chr << "\t" << (int)t.resultsSNPwindows[i][3][0] << "\t" << v.posInt << "\t" << vector_average(t.resultsSNPwindows[i][0]) << "\t" << vector_average(t.resultsSNPwindows[i][1]) << "\t" << vector_average(t.resultsSNPwindows[i][2]) << std::endl;
                }
                
                // Check if we are still in the same physical window:
                if (v.posInt > currentWindowEnd || v.posInt < currentWindowStart) {
                    t.finalizeAndOutputPhysicalWindow(i, opt::physicalWindowSize, v.chr, v.posInt, currentWindowStart, currentWindowEnd);
                }
                
                // Check if we are still in the same gene:
                if (wgAnnotation.bUpdateGene() == true) {
                    t.summariseAndOutputPerGene(i, wgAnnotation.previousGene);
                    wgAnnotation.previousGene = wgAnnotation.currentGene;
                }
            }
            
            delete c;
        }
    }
    
    return 0;
    
}

void parsePBSoptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case 'n': arg >> opt::runName; break;
            case 'm': arg >> opt::maxMissing; break;
            case OPT_ANNOT: arg >> opt::annotFile; break;
            case 'h':
                std::cout << PBS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << PBS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::PBStriosFile = argv[optind++];
    
}


bool bTrioInformativeThisSNP(const double p1, const double p2, const double p3,
                             const int n1, const int n2, const int n3,
                             const string set1, const string set2, const string set3,
                             const SetInformation& setInfo,
                             const double maxMissing) {
    
    bool bDataInformative = true;
    
    // Missing allele frequencies
    if (p1 == -1) bDataInformative = false;
    if (p2 == -1) bDataInformative = false;
    if (p3 == -1) bDataInformative = false;
    
    // Uniformative allele frequencies
    if (p1 == 0 && p2 == 0 && p3 == 0) bDataInformative = false;
    if (p1 == 1 && p2 == 1 && p3 == 1) bDataInformative = false;
    
    
    // Too much missingness based on user threshold
    if (n1 <= 1 || n2 <= 1 || n3 <= 1) bDataInformative = false;
    int set1FullSize = ploidy*(int)setInfo.popToPosMap.at(set1).size();
    int set2FullSize = ploidy*(int)setInfo.popToPosMap.at(set2).size();
    int set3FullSize = ploidy*(int)setInfo.popToPosMap.at(set2).size();
    if ( (double)n1/set1FullSize <= (1 - maxMissing) || (double)n2/set2FullSize <= (1 - maxMissing) || (double)n3/set3FullSize <= (1 - maxMissing) ) bDataInformative = false;
    
    return bDataInformative;
}
