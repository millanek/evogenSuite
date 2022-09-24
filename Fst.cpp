//
//  Fst.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

// TO DO:
// - Make sure that chromosome transitions are handled correctly
// - Make multiallelic SNPs work well and as default: * as missing data; indel with SNPs also as missing data
// - Need to deal with a window size of 1 SNP
// - Make a possibility of using annotation (make consistent with PBS?)
// - Deal with ancestral sets
// - Get a bed file with regions above a certain level

#include "Fst.hpp"

#define SUBPROGRAM "fst"

static const char *FST_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt Fst_pairs.txt\n"
"Calculate Fst statistic from a vcf file\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The Fst_pairs.txt should contain names of two populations for which the Fst will be calculated:\n"
"POP1   POP2\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_Fst_SIZE_STEP.txt\n"

"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name(s)\n"
"       -f, --fixedW sizeKb                     fixed window size (default: 10kb)\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --ancSets=ANCESTRAL_SAMPLE_SETS.txt     (optional) two sets of samples that form outgroup populations\n"
"                                               for particular Fst levels asks whether the SNPs are segregating in the outgroups\n"
"       --annot=ANNOTATION.gffExtract           (optional) gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs the location of SNPs with particular Fst levels with respect to exons, introns, UTRs, non-coding regions\n\n"
"       --regions-above=minFst                  (optional, requires -w) outputs the boundaries of regions whose Fst in windows of size set in -w is at least minFst\n"
"                                               the output file has the suffix '_fst_above_minFst.txt'\n"

"       --accessibleGenomeBED=BEDfile.bed       (optional) a bed file specifying the regions of the genome where we could call SNPs\n:"
"                                               this is used when calculating nucleotide diversity (pi) and absolute sequence divergence (d_XY) in fixed windows\n"
"       -i, --allow-indels-and-multiallelics    (optional) also calculate the PBS score for indels, and multiallelics\n"
"                                               for multiallelics, the first alternate allele is considered\n"
"                                               sites where the ALT allele is simply '*' are ignored\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1, OPT_ANNOT, OPT_ANC_SETS, OPT_REG_ABOVE, OPT_ACC_GEN_BED  };

static const char* shortopts = "hn:w:if:";

static const struct option longopts[] = {
    { "ancSets",   required_argument, NULL, OPT_ANC_SETS },
    { "window",   required_argument, NULL, 'w' },
    { "regions-above", required_argument, NULL, OPT_REG_ABOVE },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { "fixedW", required_argument, NULL, 'f' },
    { "samples",   required_argument, NULL, 's' },
    { "run-name",   required_argument, NULL, 'n' },
    { "help",   no_argument, NULL, 'h' },
    { "allow-indels-and-multiallelics",   no_argument, NULL, 'i' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string FstPairsFile;

    static string ancSets;
    static string accesibleGenBedFile;
    static int windowSize = 20;
    static int windowStep = 10;
    static int physicalWindowSize = 10000;
    static double regAbove = 0;
    static string annotFile;
    static string runName = "";
    static bool allowIndels = false;
}

int fstMain(int argc, char** argv) {
    parseFstOptions(argc, argv);
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile);
    
    // Get the population pairs
    FstPairs p(opt::FstPairsFile, opt::runName, opt::windowSize, opt::windowStep, opt::physicalWindowSize, !opt::annotFile.empty());
    
    // Read the accessible genome bed file if provided
    AccessibleGenome* ag = new AccessibleGenome(opt::accesibleGenBedFile);
    
    // Open connection to read from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    string line; std::vector<string> fields;
    int currentWindowStart = 0; int currentWindowEnd = currentWindowStart + opt::physicalWindowSize;
    string chr; string coord; int coordInt;
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
            
            //std::cerr << "Variant N:" << totalVariantNumber << std::endl;
            
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1]; coordInt = atoi(coord.c_str());
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4]; bool ignoreSite = false;
            if (altAllele == "*") ignoreSite = true;
            if (!opt::allowIndels) {
                if (refAllele.length() > 1 || altAllele.length() > 1) ignoreSite = true;
            }
            if (ignoreSite) {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            if (coordInt >= 3447835) {
                std::cerr << "coordDouble: " << coordInt << std::endl;
            }
            
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            if (coordInt >= 3447835) {
                std::cerr << "Counts created: " << std::endl;
            }
            c->getSetVariantCountsSimple(genotypes, setInfo.posToPopMap);
            if (coordInt >= 3447835) {
                std::cerr << "Summarised all counts: " << std::endl;
            }
            c->calculatePiPerVariantPerSet();
            if (coordInt >= 3447835) {
                std::cerr << "Got pi for all populations: " << std::endl;
            }
            genotypes.clear(); genotypes.shrink_to_fit();
            
            for (int i = 0; i != p.pairs.size(); i++) {
                string set1 = p.pairs[i][0]; string set2 = p.pairs[i][1];
                
                
                double p1 = c->setAAFs.at(set1); double p2 = c->setAAFs.at(set2);
                if (bPairInformativeThisSNP(p1,p2)) p.usedVars[i]++;
                else continue;
                int n1 = c->setAlleleCounts.at(set1); int n2 = c->setAlleleCounts.at(set2);
                
                if (coordInt >= 3447835) {
                    std::cerr << "p1: " << p1 << std::endl;
                    std::cerr << "p2: " << p2 << std::endl;
                }
                
                double thisSNPFstNumerator = calculateFstNumerator(p1, p2, n1, n2);
                double thisSNPFstDenominator = calculateFstDenominator(p1, p2);
                double thisSNPDxy = DxyPerSNPfromAFs(p1, p2);
                double thisSNPpi1 = c->piPerVariantPerSet.at(set1);
                double thisSNPpi2 = c->piPerVariantPerSet.at(set2);
                
                p.addSNPresultsToWindows(i,thisSNPFstNumerator,thisSNPFstDenominator, thisSNPDxy, thisSNPpi1,thisSNPpi2,coordInt);
                
                
                if (p.usedVars[i] > opt::windowSize && (p.usedVars[i] % opt::windowStep == 0)) {
                    p.finalizeAndOutputSNPwindow(i, chr, coordInt, ag);
                }
                
                // Check if we are still in the same physical window...
                if (coordInt > currentWindowEnd || coordInt < currentWindowStart) {
                    p.finalizeAndOutputPhysicalWindow(i, opt::physicalWindowSize, chr, coordInt, ag, currentWindowStart, currentWindowEnd);
                }
                
            }
        }
    }
    
    /*double Fst = calculateFst(fstNumerators, fstDenominators);
    std::cerr << "Fst: " << Fst << std::endl; */
    
    clock_t end = clock();
    double elapsed_secs = double(end - startTime) / CLOCKS_PER_SEC;
    std::cerr << "Time taken: " << elapsed_secs << std::endl;
    
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
            case OPT_ANC_SETS: arg >> opt::ancSets; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case OPT_REG_ABOVE: arg >> opt::regAbove; break;
            case OPT_ANNOT: arg >> opt::annotFile; break;
            case OPT_ACC_GEN_BED: arg >> opt::accesibleGenBedFile; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'n': arg >> opt::runName; break;
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
        std::cerr << "Error in the -w option: The window step cannot be higher than window size (you can set -w 1,1 for per variant Fst)\n";
        die = true;
    }
    
    
    if (opt::regAbove > 0 && opt::windowSize == 0) {
        std::cerr << "The window size (-w option) needs to be set for --regions-above option to work (you can set -w 1,1 for per variant Fst)\n";
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

bool bPairInformativeThisSNP(const double p1, const double p2) {
    bool bDataInformative = true;
    if (p1 == -1) bDataInformative = false;  // If any member of the pair has entirely missing data, just move on to the next pair
    if (p2 == -1) bDataInformative = false;
    if (p1 == 0 && p2 == 0) bDataInformative = false;
    if (p1 == 1 && p2 == 1) bDataInformative = false;
    return bDataInformative;
}
