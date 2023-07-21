//
//  evogenSuite.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include <iostream>
#include <stdlib.h>
#include "UtilsGeneral.hpp"

#include "AlleleFreq.hpp"
#include "CodingStats.hpp"
#include "DistanceToOutgroups.hpp"
#include "DistanceMatrix.hpp"
#include "Fst.hpp"
#include "PaintWithFixedSites.hpp"
#include "PBS.hpp"
#include "PhysicalWindowAverages.hpp"

/*
#include "process_vcf_stats.h"
#include "process_vcf_filter.h"

#include "process_vcf_variant_sharing.h"
#include "process_vcf_testing.h"
#include "process_vcf_massoko.h"
#include "process_vcf_get_sequences.h"
#include "process_vcf_sequenom.h"
#include "process_vcf_coding_sequences.h"
#include "evo_shared_variation.h"
#include "process_vcf_search_sex.h"
#include "process_vcf_mt_sequences.h"
#include "process_vcf_stats_testing.h"
#include "process_vcf_reorder.h"
#include "process_vcf_vcf_from_sequenom.h"
#include "process_vcf_fst.h"
#include "process_vcf_merge.h"
#include "process_vcf_cbs.h"
#include "process_vcf_get_aa_seq.h"
#include "process_vcf_fill_aa.h"
#include "process_vcf_join_multiFasta.h"
#include "process_vcf_shortRNA.h"
#include "process_vcf_linkGeneNames.h"
#include "evo_codingStats_from_alignment.h"
#include "evo_DNA_to_Protein.h"
#include "evo_protein_SegregatingSites.h"
#include "evo_codingSeqs_fromGenomes.h"
#include "evo_fullAnnotationExtract.h"
#include "evo_diversity_subsampling.h"
#include "remove_lowercase.h"
#include "evo_permute_codons.h"
#include "evo_ABS.h"
#include "evo_FstAgainstAll.h"
#include "evo_combineVCFs.h"
#include "evo_AlleleFeq.h"
#include "evo_agpToNewFasta.h"
#include "evo_distanceToOutgroups.h"
#include "evo_diversityFromHaps.h"
#include "evo_getInformativePairs.h"
#include "evo_findDiscordantPairs.h"
#include "evo_getInformativeReadsFromSam.h"
#include "evo_findDiscordantPairsFromSAM.h"

 */

//#define TESTING 1


#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.1 r1"


static const char *VERSION_MESSAGE =
PROGRAM_BIN " software Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR "\n"
"\n";

static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           AlleleFreq              Calculate allele frequencies per population\n"
"           CodingStats             Calculate statistics from a multiple sequence alignment of coding sequences\n"
"           DistOutgroups           Calculate the distance between ingroups/outgroups in windows\n"
"           Fst                     Calculating Fst values in windows and per-gene\n"
"           GlobalPairs             Calculating global (e.g. genome-wide) statistics from a VCF file. Main output is a distance matrix.\n"
"           PaintFixed              Find the genotypes of hybrid individuals at sites fixed between source populations\n"
"           PBS                     Calculating the Population Branch Statistics in windows and per-gene\n"
"           PhysicalWindowAverages  Getting (weighted) averages of statistics in fixed physical windows (e.g. 1kb windows)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char **argv) {
    
    if(argc <= 1)
    {
        std::cout << USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help" || command == "-h")
        {
            std::cout << USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << VERSION_MESSAGE;
            return 0;
        }
        
        if(command == "AlleleFreq")
            AFmain(argc - 1, argv + 1);
        else if (command == "CodingStats")
            getCodingStats(argc - 1, argv + 1);
        else if (command == "DistOutgroups")
            DistOutMain(argc - 1, argv + 1);
        else if(command == "Fst")
            fstMain(argc - 1, argv + 1);
        else if(command == "GlobalPairs")
            globalStatsMain(argc - 1, argv + 1);
        else if (command == "PaintFixed")
            PaintFixedMain(argc - 1, argv + 1);
        else if(command == "PBS")
            PBSmain(argc - 1, argv + 1);
        else if(command == "PhysicalWindowAverages")
            pwMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        if (iFstNum == 0) 
        return 0;
    }
}
