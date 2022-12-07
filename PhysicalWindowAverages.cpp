//
//  PhysicalWindowAverages.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 05.12.22.
//

#include "PhysicalWindowAverages.hpp"

#define SUBPROGRAM "PhysicalWindowAverages"

static const char *PW_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE(.txt|.bed) [INPUT_FILE2(.txt|.bed)]\n"
"Calculate average statistics for fixed size physical windows\n"
"The INPUT_FILE should have at least four columns:\n"
"CHR    START   END VALUE\n"
"However, additional value columns are allowed:\n"
"CHR    START   END VALUE_1  VALUE_2  ... VALUE_N\n"

"Also, multiple input files can be provided and the averages will be combined into one output file\n"
"\n"
"If the input file suffix is .bed, then the left coordinate is assumed to be zero-based (starting at 0)\n"
"With any other suffix both coordinates are assumed to be one-based (starting at 1)\n"

//"The output is POP1_POP2_Fst_SIZE_STEP.txt\n"

"\n"
HelpOption RunNameOption
"       -p, --physicalWindow SIZE               (default: 1000) fixed window size in bp\n"
"\n"
"ACCESSIBLE_GENOME options:\n"
"       --accessibleGenomeBED=BEDfile.bed       (optional) a bed file specifying the regions of the genome where we could call SNPs\n"
"                                               this is used when calculating nucleotide diversity (pi) and absolute sequence divergence (d_XY)\n"
"       -m, --minimumAccessibleBp MIN_BP        (optional; default=100) The minimum number of accessible basepairs in a window\n"
"                                               if this minimum number is not reached for a window, then the average does not get calculated\n"
"                                               instead the program outputs 'NA' for the values\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1, OPT_ACC_GEN_BED  };

static const char* shortopts = "hn:p:m:";

static const struct option longopts[] = {
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { "physicalWindow", required_argument, NULL, 'p' },
    { "run-name",   required_argument, NULL, 'n' },
    { "minimumAccessibleBp", required_argument, NULL, 'm' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static std::vector<string> inputFiles;

    static string accesibleGenBedFile;
    static int physicalWindowSize = 1000;
    static string runName = "";
    static double minAccessible = 0;
}

template <typename T> int roundToNearestValue(T num, int roundingValue)
{
    int d_i = num;
    int halfRoundingValue = roundingValue/2;
    return ((d_i % roundingValue) < halfRoundingValue) ? d_i - (d_i % roundingValue) : d_i + (roundingValue - (d_i % roundingValue));
}

int pwMain(int argc, char** argv) {
    parsePwOptions(argc, argv); 
    
    
    // Read the accessible genome bed file if provided
    AccessibleGenome* ag = new AccessibleGenome(opt::accesibleGenBedFile);
    
    // Read the input file provided
    vector<IntervalFile*> fs;
    for (int i = 0; i < opt::inputFiles.size(); i++) {
        IntervalFile* thisFile = new IntervalFile(opt::inputFiles[i]);
        fs.push_back(thisFile);
    }
    
    //std::cerr << "f->BedFeatureMap.size\t" << f->BedFeatureMap.size() << "\n";
    
    for (map<string, vector<vector<int>> >::const_iterator it = fs.front()->BedFeatureMap.begin(); it != fs.front()->BedFeatureMap.end(); it++) {
        string chromosome = it->first;
        vector< vector<int> > allCoords = it->second;
        //std::cerr << "allCoords.size()" << allCoords.size() << "\n";
        int lastCoord = allCoords[1].back();
        //std::cerr << "lastCoord\t" << lastCoord << "\t";
        int maxCoordToOutput = roundToNearestValue(lastCoord, opt::physicalWindowSize);
        for (int i = 1; i < maxCoordToOutput; i = i + opt::physicalWindowSize) {
            int physicalWindowEnd = i + opt::physicalWindowSize - 1;
            std::cout << chromosome << "\t" << i << "\t" << physicalWindowEnd << "\t";
            for (int j = 0; j < fs.size() - 1; j++) {
                print_vector(fs[j]->getMeanValuesForRegion(chromosome, i, physicalWindowEnd),std::cout,'\t',false);
                std::cout << "\t";
            }
            print_vector(fs.back()->getMeanValuesForRegion(chromosome, i, physicalWindowEnd),std::cout,'\t',true);
        }
    }
    
    return 0;
}

void parsePwOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
          //  case OPT_ANC_SETS: arg >> opt::ancSets; break;
            case OPT_ACC_GEN_BED: arg >> opt::accesibleGenBedFile; break;
            case 'p': arg >> opt::physicalWindowSize; break;
            case 'm': arg >> opt::minAccessible; break;
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << PW_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::minAccessible < 0) {
        std::cerr << "ERROR: The minimum number of accessible bp (-m option) cannot be negative \n";
        die = true;
    }
    
    if (opt::minAccessible > 0 && opt::accesibleGenBedFile.empty()) {
        std::cerr << "ERROR: The -m option only makes sense together with the --accessibleGenomeBED option\n";
        die = true;
    }
    int nInputFiles = argc - optind;
    if (nInputFiles < 1)
    {
        std::cerr << "ERROR: Too few arguments\n";
        std::cerr << "You need to provide at least one input file\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << PW_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < nInputFiles; i++) {
        // Parse the input filenames
        opt::inputFiles.push_back(argv[optind++]);
    }
}
