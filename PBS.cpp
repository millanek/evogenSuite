//
//  PBS.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "PBS.hpp"
#include "UtilsAnnotation.hpp"

#define SUBPROGRAM "PBS"

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
"       -h, --help                              display this help and exit\n"
"       -f, --fixedW sizeKb                     fixed window size (default: 10kb)\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --annot=ANNOTATION.gffExtract           (optional) gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs PBS per gene (only exons, with introns, and with 3kb upstream)\n"
"       -i, --allow-indels-and-multiallelics    (optional) also calculate the PBS score for indels, and multiallelics\n"
"                                               for multiallelics, the first alternate allele is considered\n"
"                                               sites where the ALT allele is simply '*' are ignored\n"
"       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hw:n:f:i";

enum { OPT_ANNOT  };

static const struct option longopts[] = {
    { "fixedW",   required_argument, NULL, 'f' },
    { "window",   required_argument, NULL, 'w' },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "allow-indels-and-multiallelics",   no_argument, NULL, 'i' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile; static string setsFile; static string PBStriosFile;
    //
    static string annotFile;
    static string runName = "";
    static int fixedWindowSize = 10000;
    static int windowSize = 20;
    static int windowStep = 10;
    static bool allowIndels = false;
}


int PBSmain(int argc, char** argv) {
    parsePBSoptions(argc, argv);
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile);
    
    // Initialise the PBS trios
    PBStrios t(opt::PBStriosFile, opt::runName, opt::windowSize, opt::windowStep, opt::fixedWindowSize, !opt::annotFile.empty());
    
    // Load up the annotation file if provided
    Annotation wgAnnotation(opt::annotFile, false);
    
    int currentWindowStart = 0; int currentWindowEnd = currentWindowStart + opt::fixedWindowSize;
    std::string currentGene = ""; std::string previousGene = "";
    
    std::vector<std::string> fields;
    std::clock_t startTime = std::clock();
    int totalVariantNumber = 0; int reportProgressEvery = 1000; string chr; string coord; double coordDouble;
    
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
            
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1]; coordDouble = stringToDouble(coord);
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
            
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size(), chr, coord);
            c->getSetVariantCountsSimple(genotypes, setInfo.posToPopMap);
            genotypes.clear(); genotypes.shrink_to_fit();
            // std::cerr << "Here:" << totalVariantNumber << std::endl;
            
            // find if we are in a gene:
            std::vector<string> SNPgeneDetails = wgAnnotation.getSNPgeneDetails(chr, atoi(coord.c_str()));
            if (SNPgeneDetails[0] != "") {
                currentGene = SNPgeneDetails[0];
                if (previousGene == "") previousGene = currentGene;
            }
            
            // Check if we are still in the same physical window...
            if (coordDouble > currentWindowEnd || coordDouble < currentWindowStart) {
                for (int i = 0; i != t.trios.size(); i++) {
                    int nFwSNPs1 = (int)t.resultsPhysicalWindows[i][0].size(); int nFwSNPs2 = (int)t.resultsPhysicalWindows[i][1].size(); int nFwSNPs3 = (int)t.resultsPhysicalWindows[i][2].size();
                    double PBSfw1 = 0; if (nFwSNPs1 > 0) { PBSfw1 = vector_average(t.resultsPhysicalWindows[i][0]); }
                    double PBSfw2 = 0; if (nFwSNPs2 > 0) { PBSfw2 = vector_average(t.resultsPhysicalWindows[i][1]); }
                    double PBSfw3 = 0; if (nFwSNPs3 > 0) { PBSfw3 = vector_average(t.resultsPhysicalWindows[i][2]); }
                    *t.outFilesFixedWindow[i] << chr << "\t" << currentWindowStart << "\t" << currentWindowEnd << "\t" << PBSfw1 << "\t" << PBSfw2 << "\t" << PBSfw3 << "\t" << nFwSNPs1 << "\t" << nFwSNPs2 << "\t" << nFwSNPs3 << std::endl;
                    t.resultsPhysicalWindows[i][0].clear(); t.resultsPhysicalWindows[i][1].clear(); t.resultsPhysicalWindows[i][2].clear();
                }
                if (coordDouble > currentWindowEnd) {
                    currentWindowStart = currentWindowStart + opt::fixedWindowSize; currentWindowEnd = currentWindowEnd + opt::fixedWindowSize;
                } else if (coordDouble < currentWindowStart) {
                    currentWindowStart = 0; currentWindowEnd = 0 + opt::fixedWindowSize;
                }
            }
            
            // std::cerr << coord << "\t";
            // print_vector_stream(SNPgeneDetails, std::cerr);
            // Now calculate the PBS stats:
            double p1; double p2; double p3;
            for (int i = 0; i != t.trios.size(); i++) {
                p1 = c->setAAFs.at(t.trios[i][0]); //assert(p_S1 == pS1test);
                if (p1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                p2 = c->setAAFs.at(t.trios[i][1]); //assert(p_S2 == pS2test);
                if (p2 == -1) continue;
                p3 = c->setAAFs.at(t.trios[i][2]); // assert(p_S3 == pS3test);
                if (p3 == -1) continue;
                
                if (p1 == 0 && p2 == 0 && p3 == 0) { continue; }
                if (p1 == 1 && p2 == 1 && p3 == 1) { continue; }
                t.usedVars[i]++;
                
                std::vector<double> thisSNP_PBS = calculatePBSfromAFs(p1,p2,p3,
                                                                      c->setAlleleCounts.at(t.trios[i][0]),
                                                                      c->setAlleleCounts.at(t.trios[i][1]),
                                                                      c->setAlleleCounts.at(t.trios[i][2]));
                
                t.resultsSNPwindows[i][0].push_back(thisSNP_PBS[0]); t.resultsSNPwindows[i][1].push_back(thisSNP_PBS[1]); t.resultsSNPwindows[i][2].push_back(thisSNP_PBS[2]);
                t.resultsSNPwindows[i][3].push_back(coordDouble);
                t.resultsSNPwindows[i][0].pop_front(); t.resultsSNPwindows[i][1].pop_front(); t.resultsSNPwindows[i][2].pop_front(); t.resultsSNPwindows[i][3].pop_front();
                
                t.resultsPhysicalWindows[i][0].push_back(thisSNP_PBS[0]); t.resultsPhysicalWindows[i][1].push_back(thisSNP_PBS[1]); t.resultsPhysicalWindows[i][2].push_back(thisSNP_PBS[2]);
                
                if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] != "") {
                    if (SNPgeneDetails[1] == "exon") {
                        t.resultsGenes[i][0].push_back(thisSNP_PBS[0]); t.resultsGenes[i][1].push_back(thisSNP_PBS[1]); t.resultsGenes[i][2].push_back(thisSNP_PBS[2]);
                    } else if (SNPgeneDetails[1] == "intron") {
                        t.resultsGenes[i][3].push_back(thisSNP_PBS[0]); t.resultsGenes[i][4].push_back(thisSNP_PBS[1]); t.resultsGenes[i][5].push_back(thisSNP_PBS[2]);
                    } else if (SNPgeneDetails[1] == "promoter") {
                        t.resultsGenes[i][6].push_back(thisSNP_PBS[0]); t.resultsGenes[i][7].push_back(thisSNP_PBS[1]); t.resultsGenes[i][8].push_back(thisSNP_PBS[2]);
                    }
                }}
                
                if (t.usedVars[i] > opt::windowSize && (t.usedVars[i] % opt::windowStep == 0)) {
                    // std::cerr << PBSresults[i][0][0] << std::endl;
                    *t.outFiles[i] << chr << "\t" << (int)t.resultsSNPwindows[i][3][0] << "\t" << coord << "\t" << vector_average(t.resultsSNPwindows[i][0]) << "\t" << vector_average(t.resultsSNPwindows[i][1]) << "\t" << vector_average(t.resultsSNPwindows[i][2]) << std::endl;
                }
                // }
            }
            if (!opt::annotFile.empty()) { if (previousGene != "" && currentGene != previousGene) {
                t.summariseAndOutputPerGene(previousGene);
                previousGene = currentGene;
            }}
            
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
            case 'f': arg >> opt::fixedWindowSize; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case 'n': arg >> opt::runName; break;
            case 'i': opt::allowIndels = true; break;
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

void PBStrios::summariseAndOutputPerGene(string geneName) {
    for (int i = 0; i != trios.size(); i++) {
        int nExonSNPs = (int)resultsGenes[i][0].size(); int nIntronSNPs = (int)resultsGenes[i][3].size(); int nPromoterSNPs = (int)resultsGenes[i][6].size();
        double ExonPBS1 = 0; double ExonPBS2 = 0; double ExonPBS3 = 0;
        if (nExonSNPs > 0) { ExonPBS1 = vector_average(resultsGenes[i][0]); ExonPBS2 = vector_average(resultsGenes[i][1]); ExonPBS3 = vector_average(resultsGenes[i][2]);}
        double IntronPBS1 = 0; double IntronPBS2 = 0; double IntronPBS3 = 0;
        if (nIntronSNPs > 0) { IntronPBS1 = vector_average(resultsGenes[i][3]); IntronPBS2 = vector_average(resultsGenes[i][4]); IntronPBS3 = vector_average(resultsGenes[i][5]);}
        double PromoterPBS1 = 0; double PromoterPBS2 = 0; double PromoterPBS3 = 0;
        if (nPromoterSNPs > 0) { PromoterPBS1 = vector_average(resultsGenes[i][6]); PromoterPBS2 = vector_average(resultsGenes[i][7]); PromoterPBS3 = vector_average(resultsGenes[i][8]);}
        *outFilesGenes[i] << geneName << "\t" << nExonSNPs << "\t" << nIntronSNPs << "\t" << nPromoterSNPs << "\t" << ExonPBS1 << "\t" << ExonPBS2 << "\t" << ExonPBS3 << "\t" << IntronPBS1 << "\t" << IntronPBS2 << "\t" << IntronPBS3 << "\t" << PromoterPBS1 << "\t" << PromoterPBS2 << "\t" << PromoterPBS3 << std::endl;
        for (int j = 0; j <= 8; j++) { resultsGenes[i][j].clear(); }
    }
}
