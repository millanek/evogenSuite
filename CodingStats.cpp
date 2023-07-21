//
//  CodingStats.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 28.06.23.
//

#include "CodingStats.hpp"

#define SUBPROGRAM "CodingStats"
#define MIN_SETS 1

static const char *CODINGSTATS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-a multiple_alignment.fa | -l list_of_multiple_aligment_files.txt> SETS.txt\n"
"Calculate statistics out of multiple aligments of gene sequences (e.g. from the output of 'evo getCodingSeq')\n"
"No checks are performed on the sequences\n"
"\n"
HelpOption RunNameOption
"       -t,   --tStV RATIO                      observed genome-wide tS/tV ratio in the dataset\n"
"       -a,   --alignment FILE.fa               a multiple alignment file (either -a or -l required)\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"       -u,   --nonCodingNull                   the alignment(s) contain non-coding sequences - used to derive null distributions of the statistics\n"
"       --pNofGroups=GROUPS-SAMPLES.txt         outputs the average pN betwen the first two subgroups defined in GROUPS-SAMPLES.txt and\n"
"       -f,   --Fst THRESHOLD                   (default: 0.3) Threshold for between-population comparisons\n"
"                                               between the two subgroups joined and the third group\n"
"       --genomeWide_dXY=MATRIX.txt             a matrix that can be used to normalise the pairwise scores to be\n"
"                                               per-unit of sequence divergence\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hf:a:l:t:n:u";

enum { OPT_DXY_MATRIX, OPT_PN_GROUPS };

static const struct option longopts[] = {
    { "Fst",   required_argument, NULL, 'f' },
    { "alignment",   required_argument, NULL, 'a' },
    { "tStV",   required_argument, NULL, 't' },
    { "listOfFiles",   required_argument, NULL, 'l' },
    { "nonCodingNull",   no_argument, NULL, 'u' },
    { "run-name",   required_argument, NULL, 'n' },
    { "pNofGroups", required_argument, NULL, OPT_PN_GROUPS },
    { "genomeWide_dXY",   required_argument, NULL, OPT_DXY_MATRIX },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
        static string setsFile = "";
        static string alignmentFile = "";
        static string alignmentListFile = "";
        static string pNgroupsFile = "";
        static string runName = "RN";
        string genomeWide_DxyMatrixFile = "";
        static double FstThreshold = 0.3;
        static bool nonCodingNull = false;
        static double tStVratio = 0.5;      // Under equal mutation probability,
                                            //half transitions will be observed compared with transversions
}

std::vector<string> getAlignmentFilenames(const string& alignmentListFile) {
    std::vector<string> allAligmentFiles;
    std::ifstream* alignmentList = new std::ifstream(opt::alignmentListFile.c_str());
    std::string line;
    while (getline(*alignmentList, line)) {
        allAligmentFiles.push_back(line);
    }
    return allAligmentFiles;
}


int getCodingStats(int argc, char** argv) {
    parseCodingStatsOptions(argc, argv);
    // string pcaVectorsFileName;
    
    SetInformation setInfo(opt::setsFile, MIN_SETS);
    int nCombinations = choose((int)setInfo.populations.size(),2);
    
    std::cerr << "Calculating gene coding statistics" << std::endl;
    std::cerr << "for the alignments in: " << opt::alignmentListFile << std::endl;
    std::vector<string> allAligmentFiles = getAlignmentFilenames(opt::alignmentListFile);
    
    string statsFileName = stripExtension(opt::alignmentListFile) + "_" + opt::runName + "_NewStats.txt";
    std::ofstream* statsFile = new std::ofstream(statsFileName.c_str());
    *statsFile << "alignment" << "\t" << "ntLengh" << "\t" << "pN" << "\t" << "pS" << "\t" << "pNstdErr" << "\t" << "pSstdErr" << "\t" << "pNpSstdErr" << std::endl;
    
    std::vector<std::ofstream*> fstOutFiles;
    for (int i = 0; i < setInfo.populations.size(); i++) {
        for (int j = i+1; j < setInfo.populations.size(); j++) {
            string fstFileName = stripExtension(opt::alignmentListFile) + "_" + opt::runName + "_FstStats_" + setInfo.populations[i] + "_" + setInfo.populations[j] + "_" + numToString(opt::FstThreshold) + ".txt";
            std::ofstream* fstFile = new std::ofstream(fstFileName.c_str());
            *fstFile << "alignment" << "\t" << "ntLengh" << "\t" << "N_lowFst" << "\t" << "S_lowFst" << "\t" << "N_highFst" << "\t" << "S_highFst" << std::endl;
            fstOutFiles.push_back(fstFile);
        }
    }
    
    for (std::vector<std::string>::size_type i = 0; i != allAligmentFiles.size(); i++) {
        MultipleAlignment* ma = new MultipleAlignment(allAligmentFiles[i], nCombinations);
        setInfo.linkSetsAndVCFpositions(ma->sequenceNames, false);
        // std::cerr << "setInfo.posToPopMap[0]: " << setInfo.posToPopMap[0] << std::endl;
        ma->getStatsAllSequences(opt::tStVratio, false, setInfo, opt::FstThreshold);
        ma->printAlignmentStats(std::cout);
        ma->printAlignmentStats(*statsFile);
        for (int j = 0; j < nCombinations; j++) {
            ma->printFstStats(*fstOutFiles[j], j);
        }
        delete ma;
        
      //  getStatsBothPhasedHaps(ma->allSeqs, ma->allSeqsH2, statsThisGene, combinedVectorForPCA, sets, opt::tStVratio, opt::nonCodingNull);
       // std::cerr << "got stats for: " << allAligmentFiles[i] << std::endl;
        /*
        if (opt::pNgroupsFile != "") {
            int ns1 = (int)sets->set1Loci.size(); int ns2 = (int)sets->set2Loci.size(); int ns3 = (int)sets->set3Loci.size();
            statsThisGene.push_back(numToString(sets->withinSet1pN/(2*ns1*(ns1-1))));
            statsThisGene.push_back(numToString(sets->withinSet2pN/(2*ns2*(ns2-1))));
            statsThisGene.push_back(numToString(sets->withinSet3pN/(2*ns3*(ns3-1))));
            statsThisGene.push_back(numToString(sets->set1vsSet2pN/(2*ns1*ns2)));
            statsThisGene.push_back(numToString(sets->sets1and2vsSet3pN/(2*(ns1+ns2)*ns3)));
            statsThisGene.push_back(numToString(sets->withinSet1andSet2pN/(2*(ns1+ns2)*(ns1+ns2-1))));
            sets->withinSet1pN = 0;
            sets->withinSet2pN = 0;
            sets->withinSet3pN = 0;
            sets->set1vsSet2pN = 0;
            sets->sets1and2vsSet3pN = 0;
            sets->withinSet1andSet2pN = 0;
        }
         */
    }
    
    return 0;
}

void parseCodingStatsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'f': arg >> opt::FstThreshold; break;
            case 'a': arg >> opt::alignmentFile; break;
            case 'l': arg >> opt::alignmentListFile; break;
            case 't': arg >> opt::tStVratio; break;
            case 'u': opt::nonCodingNull = true; break;
            case 'n': arg >> opt::runName; break;
            case OPT_DXY_MATRIX: arg >> opt::genomeWide_DxyMatrixFile; break;
            case OPT_PN_GROUPS: arg >> opt::pNgroupsFile; break;
            case 'h':
                std::cout << CODINGSTATS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::alignmentFile == "" && opt::alignmentListFile == "") {
        std::cerr << "Either -a or -l options must be specified\n";
        die = true;
    }
    
    if (opt::alignmentFile != "" && opt::alignmentListFile != "") {
        std::cerr << "The -a and -l options can't both be specified\n";
        die = true;
    }
    
    if (opt::FstThreshold < 0 || opt::FstThreshold > 1) {
        std::cerr << "The -f (--Fst) option can only have the values between 0 and 1\n";
        std::cerr << "At the moment I don't support other than haploid and diploid species\n";
        die = true;
    }
    
    if (opt::tStVratio <= 0) {
        std::cerr << "The tS/tV ratio must be positive\n";
        die = true;
    }
    
    
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << CODINGSTATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    opt::setsFile = argv[optind++];
    
}


bool isUnrecognisedCodon(string codon) {
    if (getAminoAcid(codon) == "Unrecognised codon") {
        return true;
    } else {
        return false;
    }
}

bool compareFirstInPair(std::pair<int, string> a, std::pair<int, string> b) {
    return a.first > b.first;
}
