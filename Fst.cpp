//
//  Fst.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

// TO DO:
// - Make sure that chromosome transitions are handled correctly
// - Need to deal with a window size of 1 SNP
// - When using annotation, make possible that Dxy and Pi are also calculated and that accessible BP per gene/feature are output
// - Deal with ancestral sets
// - Get a bed file with regions above a certain level
// - If a population has entirely missing data for a SNP, pi will be 0

#include "Fst.hpp"

#define SUBPROGRAM "Fst"
#define MIN_SETS 2

static const char *FST_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt Fst_pairs.txt\n"
"Calculate Fst statistic from a vcf file\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The Fst_pairs.txt should contain names of two populations for which the Fst will be calculated:\n"
"POP1   POP2\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_Fst_SIZE_STEP.txt\n"

"\n"
HelpOption RunNameOption MaxMissOption
"       -f, --fixedW SIZE                       (default: 0) fixed window size - will do a fixed window calculation if SIZE > 0bp\n"
"       -w SIZE,STEP --window=SIZE,STEP         (default: 20,10) the parameters of the sliding window: contains SIZE SNPs and move by STEP\n"
//"       --ancSets=ANCESTRAL_SAMPLE_SETS.txt     (optional) two sets of samples that form outgroup populations\n"
//"                                               for particular Fst levels asks whether the SNPs are segregating in the outgroups\n"
"       --annot=ANNOTATION.gffExtract           (optional) gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"       -g, --genMap=GENETIC_MAP.optimise           (optional) A genetic map in the pyRho format\n"
//"       --regions-above=minFst                  (optional, requires -w) outputs the boundaries of regions whose Fst in windows of size set in -w is at least minFst\n"
//"                                               the output file has the suffix '_fst_above_minFst.txt'\n"
"       --accessibleGenomeBED=BEDfile.bed       (optional) a bed file specifying the regions of the genome where we could call SNPs\n"
"                                               this is used when calculating nucleotide diversity (pi) and absolute sequence divergence (d_XY)\n"
"       -z, --noRoundingToZero                  (optional) Do not round negative Fst values to zero\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1, OPT_ANNOT, OPT_ANC_SETS, OPT_REG_ABOVE, OPT_ACC_GEN_BED  };

static const char* shortopts = "hn:w:f:m:g:z";

static const struct option longopts[] = {
//    { "ancSets",   required_argument, NULL, OPT_ANC_SETS },
    { "window",   required_argument, NULL, 'w' },
    { "regions-above", required_argument, NULL, OPT_REG_ABOVE },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { "fixedW", required_argument, NULL, 'f' },
    { "run-name",   required_argument, NULL, 'n' },
    { "maxMissing", required_argument, NULL, 'm' },
    { "genMap",   required_argument, NULL, 'g' },
    { "help",   no_argument, NULL, 'h' },
    { "noRoundingToZero",   no_argument, NULL, 'z' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string FstPairsFile;

    static string accesibleGenBedFile; static string annotFile; static string genMapFile;
    static int windowSize = 20; static int windowStep = 10;
    static int physicalWindowSize = 0;
    static string runName = "";
    static double maxMissing = 0.2;
    static bool bZeroRounding = true;
    static bool bMakeRegionsBed = false; static double regAbove = 1;
    //static string ancSets;
}

int fstMain(int argc, char** argv) {
    parseFstOptions(argc, argv);
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile, MIN_SETS);
    
    // Get the population pairs
    FstPairs p(opt::FstPairsFile, opt::runName, opt::windowSize, opt::windowStep, opt::physicalWindowSize, !opt::annotFile.empty(), opt::bMakeRegionsBed, opt::regAbove, !opt::genMapFile.empty());
    
    // Read the accessible genome bed file if provided
    AccessibleGenome* ag = new AccessibleGenome(opt::accesibleGenBedFile);
    
    // Read the genetic map if provided
    RecombinationMap* r = new RecombinationMap(opt::genMapFile);
    
    // Load up the annotation file if provided
    Annotation wgAnnotation(opt::annotFile, false);
    
    // Open connection to read from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    string line; std::vector<string> fields;
    int currentWindowStart = 0; int currentWindowEnd = currentWindowStart + opt::physicalWindowSize;
    string chr; string coord;
    int totalVariantNumber = 0; int reportProgressEvery = 10000; std::clock_t startTime = std::clock();
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
            
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            
           // std::cerr << "Variant N:" << totalVariantNumber << std::endl;
            
            fields = split(line, '\t');
            
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            
         //  if (totalVariantNumber == 6015) std::cerr <<  line << std::endl;
            
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            //std::cerr << "Counts created: " << std::endl;
            c->getSetVariantCounts(genotypes, setInfo.posToPopMap, v);
            //std::cerr << "Summarised all counts: " << std::endl;
            c->calculatePiPerVariantPerSet();
            //if (totalVariantNumber >= 43) std::cerr << "Got pi for all populations: " << std::endl;
            genotypes.clear(); genotypes.shrink_to_fit();
            
            for (int i = 0; i != p.getPairs().size(); i++) {
                string set1 = p.getPairs()[i][0]; string set2 = p.getPairs()[i][1];
                
                double p1 = c->setAAFs.at(set1)[0]; double p2 = c->setAAFs.at(set2)[0];
                int n1 = c->setRefCounts.at(set1) + c->setAltAlleleCounts.at(set1)[0];
                int n2 = c->setRefCounts.at(set2) + c->setAltAlleleCounts.at(set2)[0];
                int set1FullSize = ploidy*(int)setInfo.popToPosMap.at(set1).size();
                int set2FullSize = ploidy*(int)setInfo.popToPosMap.at(set2).size();
                
                if (bPairInformativeThisSNP(p1,p2,n1,n2,set1FullSize,set2FullSize, opt::maxMissing)) p.usedVars[i]++;
                else continue;
                
           /*     if (totalVariantNumber == 6015) {
                    std::cerr << "p1: " << p1 << std::endl;
                    std::cerr << "p2: " << p2 << std::endl;
                    std::cerr << "n1: " << n1 << std::endl;
                    std::cerr << "n2: " << n2 << std::endl;
                    std::cerr << "totalVariantNumber: " << totalVariantNumber << std::endl;
                } */
                
                double thisSNPFstNumerator = calculateFstNumerator(p1, p2, n1, n2);
                double thisSNPFstDenominator = calculateFstDenominator(p1, p2);
                double thisSNPDxy = DxyPerSNPfromSetAlleles(c, set1, set2);
                double thisSNPpi1 = c->piPerVariantPerSet.at(set1);
                double thisSNPpi2 = c->piPerVariantPerSet.at(set2);
                
                
                /*  if (isnan(thisSNPFstNumerator) || isnan(thisSNPFstDenominator)) {
                     std::cerr << "thisSNPFstNumerator: " << thisSNPFstNumerator << std::endl;
                     std::cerr << "thisSNPFstDenominator: " << thisSNPFstDenominator << std::endl;
                     std::cerr << "totalVariantNumber: " << totalVariantNumber << std::endl;
                 } */
                
                p.addSNPresultsToWindows(i,thisSNPFstNumerator,thisSNPFstDenominator, thisSNPDxy, thisSNPpi1,thisSNPpi2,v.posInt);
                
                if (opt::physicalWindowSize > 0) p.addSNPresultsToPhysicalWindows(i,thisSNPFstNumerator,thisSNPFstDenominator, thisSNPDxy, thisSNPpi1,thisSNPpi2,v.posInt);
                
                if (wgAnnotation.currentGene != "") p.addSNPresultsToGene(i, thisSNPFstNumerator,thisSNPFstDenominator, wgAnnotation);
                
                if (p.usedVars[i] > opt::windowSize && (p.usedVars[i] % opt::windowStep == 0)) {
                    p.finalizeAndOutputSNPwindow(i, v.chr, v.posInt, ag, r, opt::bZeroRounding);
                }
                
                // Check if we are still in the same physical window...
                if (opt::physicalWindowSize > 0) {
                    if (v.posInt > currentWindowEnd || v.posInt < currentWindowStart) {
                        p.finalizeAndOutputPhysicalWindow(i, opt::physicalWindowSize, v.chr, v.posInt, ag, currentWindowStart, currentWindowEnd, r, opt::bZeroRounding);
                    }
                }
                
                // Check if we are still in the same gene:
                if (wgAnnotation.bUpdateGene() == true) {
                    p.summariseAndOutputPerGene(i, wgAnnotation.previousGene);
                    wgAnnotation.previousGene = wgAnnotation.currentGene;
                }
                
            }
            
            delete c;
        }
    }
    
    /*double Fst = calculateFst(fstNumerators, fstDenominators);
    std::cerr << "Fst: " << Fst << std::endl; */
    
    clock_t end = clock();
    double elapsed_secs = double(end - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << elapsed_secs << std::endl;
    
    return 0;
    
}

void parseFstOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
          //  case OPT_ANC_SETS: arg >> opt::ancSets; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                if (windowSizeStep.size() != 2) {
                    std::cerr << "ERROR: The -w option requires two numbers separated by a comma; e.g. '-w 20,10'\n";
                    die = true;
                }
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case OPT_REG_ABOVE: opt::bMakeRegionsBed = true; arg >> opt::regAbove; break;
            case OPT_ANNOT: arg >> opt::annotFile; break;
            case OPT_ACC_GEN_BED: arg >> opt::accesibleGenBedFile; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'g': arg >> opt::genMapFile; break;
            case 'm': arg >> opt::maxMissing; break;
            case 'n': arg >> opt::runName; break;
            case 'z': opt::bZeroRounding = false; break;
            case '?': die = true; break;
            case 'h': std::cout << FST_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    if (argc - optind < 3)
    {
        std::cerr << "too few arguments\n";
        std::cerr << "You need to provide a VCF file, a POPULATIONS.txt file, and a Fst_pairs.txt\n";
        die = true;
    }
    
    if (opt::windowStep > opt::windowSize) {
        std::cerr << "ERROR: In the -w option the window step cannot be higher than window size (you can set -w 1,1 for per variant Fst)\n";
        die = true;
    }
    
    validateMissOption(opt::maxMissing);
    
    if (opt::regAbove <= 0 ||  opt::regAbove > 1) {
        std::cerr << "Error in the --regions-above option: The value has to be in the interval (0,1) \n";
        die = true;
    }
    
    
    if (die) {
        std::cout << "\n" << FST_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::FstPairsFile = argv[optind++];
}

bool bPairInformativeThisSNP(const double p1, const double p2, const int n1, const int n2,
                             const int set1FullSize, const int set2FullSize, const double maxMissing) {
    bool bDataInformative = true;
    if (p1 == -1) bDataInformative = false;  // If any member of the pair has entirely missing data, just move on to the next pair
    if (p2 == -1) bDataInformative = false;
    if (p1 == 0 && p2 == 0) bDataInformative = false;
    if (p1 == 1 && p2 == 1) bDataInformative = false;
    if (n1 <= 1 || n2 <= 1) bDataInformative = false;
    if ( (double)n1/set1FullSize <= (1 - maxMissing) || (double)n2/set2FullSize <= (1 - maxMissing) ) bDataInformative = false;
    return bDataInformative;
}
