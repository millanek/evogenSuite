//
//  UtilsCodingStats.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 28.06.23.
//

#ifndef UtilsCodingStats_hpp
#define UtilsCodingStats_hpp

#include "UtilsGeneral.hpp"

class CDSComparisonMatrices {
public:
    CDSComparisonMatrices(int numSamples) {
        initialize_matrix_double(N_d_jk, numSamples); initialize_matrix_double(N_jk, numSamples);
        initialize_matrix_double(S_d_jk, numSamples); initialize_matrix_double(S_jk, numSamples);
        initialize_matrix_double(tS_N_jk, numSamples); initialize_matrix_double(tS_N_jk, numSamples);
        initialize_matrix_double(tS_S_jk, numSamples); initialize_matrix_double(tS_S_jk, numSamples);
        initialize_matrix_double(tV_N_jk, numSamples); initialize_matrix_double(tV_N_jk, numSamples);
        initialize_matrix_double(tV_S_jk, numSamples); initialize_matrix_double(tV_S_jk, numSamples);
    }
    std::vector<std::vector<double> > N_d_jk; std::vector<std::vector<double> > N_jk;
    std::vector<std::vector<double> > S_d_jk; std::vector<std::vector<double> > S_jk;
    // Separate matrices to incorporate unequal tS/tV mutation probabilities:
    std::vector<std::vector<double> > tS_N_jk; std::vector<std::vector<double> > tS_S_jk;
    std::vector<std::vector<double> > tV_N_jk; std::vector<std::vector<double> > tV_S_jk;
};

static std::unordered_map<std::string, bool> s_mapCodonPairToSynonymous =
{
    // Phe
    { "TTTTTC", true }, { "TTCTTT", true },
    // Leu
    { "TTATTG", true }, { "TTACTA", true }, { "TTGTTA", true }, { "TTGCTG", true },
    { "CTATTA", true }, { "CTACTC", true }, { "CTACTG", true }, { "CTACTT", true },
    { "CTCCTA", true }, { "CTCCTG", true }, { "CTCCTT", true },
    { "CTGTTG", true }, { "CTGCTC", true }, { "CTGCTA", true }, { "CTGCTT", true },
    { "CTTCTA", true }, { "CTTCTC", true }, { "CTTCTG", true },
    // Ile
    { "ATAATC", true }, { "ATAATT", true },
    { "ATCATA", true }, { "ATCATT", true },
    { "ATTATA", true }, { "ATTATC", true },
    // Val
    { "GTAGTC", true }, { "GTAGTG", true }, { "GTAGTT", true },
    { "GTCGTA", true }, { "GTCGTG", true }, { "GTCGTT", true },
    { "GTGGTA", true }, { "GTGGTC", true }, { "GTGGTT", true },
    { "GTTGTA", true }, { "GTTGTC", true }, { "GTTGTG", true },
    // Ser
    { "TCATCC", true }, { "TCATCG", true }, { "TCATCT", true },
    { "TCCTCA", true }, { "TCCTCG", true }, { "TCCTCT", true },
    { "TCGTCA", true }, { "TCGTCC", true }, { "TCGTCT", true },
    { "TCTTCA", true }, { "TCTTCC", true }, { "TCTTCG", true },
    { "AGCAGT", true }, { "AGTAGC", true },
    // Pro
    { "CCACCC", true }, { "CCACCG", true }, { "CCACCT", true },
    { "CCCCCA", true }, { "CCCCCG", true }, { "CCCCCT", true },
    { "CCGCCA", true }, { "CCGCCC", true }, { "CCGCCT", true },
    { "CCTCCA", true }, { "CCTCCC", true }, { "CCTCCG", true },
    // Thr
    { "ACAACC", true }, { "ACAACG", true }, { "ACAACT", true },
    { "ACCACA", true }, { "ACCACG", true }, { "ACCACT", true },
    { "ACGACA", true }, { "ACGACC", true }, { "ACGACT", true },
    { "ACTACA", true }, { "ACTACC", true }, { "ACTACG", true },
    // Ala
    { "GCAGCC", true }, { "GCAGCG", true }, { "GCAGCT", true },
    { "GCCGCA", true }, { "GCCGCG", true }, { "GCCGCT", true },
    { "GCGGCA", true }, { "GCGGCC", true }, { "GCGGCT", true },
    { "GCTGCA", true }, { "GCTGCC", true }, { "GCTGCG", true },
    // Tyr
    { "TACTAT", true }, { "TATTAC", true },
    // Stop
    { "TAATAG", true }, { "TAATGA", true },
    { "TAGTAA", true }, { "TGATAA", true },
    // His
    { "CACCAT", true }, { "CATCAC", true },
    //
    { "CAACAG", true }, { "CAGCAA", true },
    //
    { "AACAAT", true }, { "AATAAC", true },
    //
    { "AAAAAG", true }, { "AAGAAA", true },
    //
    { "GACGAT", true }, { "GATGAC", true },
    //
    { "GAAGAG", true }, { "GAGGAA", true },
    // Cys
    { "TGCTGT", true }, { "TGTTGC", true },
    // Arg
    { "CGACGC", true }, { "CGACGG", true }, { "CGACGT", true }, { "CGAAGA", true },
    { "CGCCGA", true }, { "CGCCGG", true }, { "CGCCGT", true },
    { "CGGCGC", true }, { "CGGCGA", true }, { "CGGCGT", true }, { "CGGAGG", true },
    { "CGTCGA", true }, { "CGTCGC", true }, { "CGTCGG", true },
    { "AGACGA", true }, { "AGAAGG", true }, { "AGGAGA", true }, { "AGGCGG", true },
    // Gly
    { "GGAGGC", true }, { "GGAGGG", true }, { "GGAGGT", true },
    { "GGCGGA", true }, { "GGCGGG", true }, { "GGCGGT", true },
    { "GGGGGA", true }, { "GGGGGC", true }, { "GGGGGT", true },
    { "GGTGGA", true }, { "GGTGGC", true }, { "GGTGGG", true },
};

static std::unordered_map<std::string, double> s_mapCodonToExpDist =
{
    { "TTT", 8.0/3.0 }, { "TTC", 8.0/3.0 },
    { "TTA", 7.0/3.0 }, { "TTG", 7.0/3.0 },                                             // Leu
    { "CTA", 5.0/3.0 }, { "CTG", 5.0/3.0 }, { "CTC", 2.0 }, { "CTT", 2.0 },             // Leu
    { "ATA", 7.0/3.0 }, { "ATC", 7.0/3.0 }, { "ATT", 7.0/3.0 },
    { "ATG", 3.0 },
    { "GTA", 2.0 }, { "GTC", 2.0 }, { "GTG", 2.0 }, { "GTT", 2.0 },
    { "TCA", 2.0 }, { "TCC", 2.0 }, { "TCG", 2.0 }, { "TCT", 2.0 },                     // Ser
    { "AGC", 8.0/3.0 }, { "AGT", 8.0/3.0 },                                             // Ser
    { "CCA", 2.0 }, { "CCC", 2.0 }, { "CCG", 2.0 }, { "CCT", 2.0 },
    { "ACA", 2.0 }, { "ACC", 2.0 }, { "ACG", 2.0 }, { "ACT", 2.0 },
    { "GCA", 2.0 }, { "GCC", 2.0 }, { "GCG", 2.0 }, { "GCT", 2.0 },
    { "TAC", 8.0/3.0 }, { "TAT", 8.0/3.0 },
    { "TAA", 7.0/3.0 }, { "TAG", 8.0/3.0 }, { "TGA", 8.0/3.0 },                         // Stop
    { "CAC", 8.0/3.0 }, { "CAT", 8.0/3.0 },
    { "CAA", 8.0/3.0 }, { "CAG", 8.0/3.0 },
    { "AAT", 8.0/3.0 }, { "AAC", 8.0/3.0 },
    { "AAA", 8.0/3.0 }, { "AAG", 8.0/3.0 },
    { "GAC", 8.0/3.0 }, { "GAT", 8.0/3.0 },
    { "GAA", 8.0/3.0 }, { "GAG", 8.0/3.0 },
    { "TGC", 8.0/3.0 }, { "TGT", 8.0/3.0 },
    { "TGG", 3.0 },
    { "CGA", 5.0/3.0 }, { "CGG", 5.0/3.0 }, { "CGC", 2.0 }, { "CGT", 2.0 },             // Arg
    { "AGA", 7.0/3.0 }, { "AGG", 7.0/3.0 },                                             // Arg
    { "GGA", 2.0 }, { "GGC", 2.0 }, { "GGG", 2.0 }, { "GGT", 2.0 }
};

static std::unordered_map<std::string, double> s_mapCodonToExpTs =
{
    { "TTT", 2.0/3.0 }, { "TTC", 2.0/3.0 },
    { "TTA", 1.0/3.0 }, { "TTG", 1.0/3.0 },                                             // Leu
    { "CTA", 1.0/3.0 }, { "CTG", 1.0/3.0 }, { "CTC", 2.0/3.0 }, { "CTT", 2.0/3.0 },     // Leu
    { "ATA", 1.0 }, { "ATC", 2.0/3.0 }, { "ATT", 2.0/3.0 },                             // Ile
    { "ATG", 1.0 },                                                                     // Met
    { "GTA", 2.0/3.0 }, { "GTC", 2.0/3.0 }, { "GTG", 2.0/3.0 }, { "GTT", 2.0/3.0 },     // Val
    { "TCA", 2.0/3.0 }, { "TCC", 2.0/3.0 }, { "TCG", 2.0/3.0 }, { "TCT", 2.0/3.0 },     // Ser
    { "AGC", 2.0/3.0 }, { "AGT", 2.0/3.0 },                                             // Ser
    { "CCA", 2.0/3.0 }, { "CCC", 2.0/3.0 }, { "CCG", 2.0/3.0 }, { "CCT", 2.0/3.0 },     // Pro
    { "ACA", 2.0/3.0 }, { "ACC", 2.0/3.0 }, { "ACG", 2.0/3.0 }, { "ACT", 2.0/3.0 },     // Thr
    { "GCA", 2.0/3.0 }, { "GCC", 2.0/3.0 }, { "GCG", 2.0/3.0 }, { "GCT", 2.0/3.0 },     // Ala
    { "TAC", 2.0/3.0 }, { "TAT", 2.0/3.0 },                                             // Tyr
    { "TAA", 1.0/3.0 }, { "TAG", 2.0/3.0 }, { "TGA", 2.0/3.0 },                         // Stop
    { "CAC", 2.0/3.0 }, { "CAT", 2.0/3.0 },                                             // His
    { "CAA", 2.0/3.0 }, { "CAG", 2.0/3.0 },                                             // Gln
    { "AAT", 2.0/3.0 }, { "AAC", 2.0/3.0 },                                             // Asn
    { "AAA", 2.0/3.0 }, { "AAG", 2.0/3.0 },                                             // Lys
    { "GAC", 2.0/3.0 }, { "GAT", 2.0/3.0 },                                             // Asp
    { "GAA", 2.0/3.0 }, { "GAG", 2.0/3.0 },                                             // Glu
    { "TGC", 2.0/3.0 }, { "TGT", 2.0/3.0 },                                             // Cys
    { "TGG", 1.0 },                                                                     // Trp
    { "CGA", 2.0/3.0 }, { "CGG", 2.0/3.0 }, { "CGC", 2.0/3.0 }, { "CGT", 2.0/3.0 },     // Arg
    { "AGA", 2.0/3.0 }, { "AGG", 2.0/3.0 },                                             // Arg
    { "GGA", 2.0/3.0 }, { "GGC", 2.0/3.0 }, { "GGG", 2.0/3.0 }, { "GGT", 2.0/3.0 }      // Gly
};

inline std::string getAminoAcid(const std::string& nuc)
{
    if (nuc == "TTT" || nuc == "TTC")
        return "Phe";
    if (nuc == "TTA" || nuc == "TTG" || nuc == "CTA" || nuc == "CTC" || nuc == "CTG" || nuc == "CTT")
        return "Leu";
    if (nuc == "ATA" || nuc == "ATC" || nuc == "ATT")
        return "Ile";
    if (nuc == "ATG")
        return "Met";
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return "Val";
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT")
        return "Ser";
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return "Pro";
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return "Thr";
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return "Ala";
    if (nuc == "TAC" || nuc == "TAT")
        return "Tyr";
    if (nuc == "TAA" || nuc == "TAG" || nuc == "TGA")
        return "Stop";
    if (nuc == "CAC" || nuc == "CAT")
        return "His";
    if (nuc == "CAA" || nuc == "CAG")
        return "Gln";
    if (nuc == "AAC" || nuc == "AAT")
        return "Asn";
    if (nuc == "AAA" || nuc == "AAG")
        return "Lys";
    if (nuc == "GAC" || nuc == "GAT")
        return "Asp";
    if (nuc == "GAA" || nuc == "GAG")
        return "Glu";
    if (nuc == "TGC" || nuc == "TGT")
        return "Cys";
    if (nuc == "TGG")
        return "Trp";
    if (nuc == "CGA" || nuc == "CGC" || nuc == "CGG" || nuc == "CGT")
        return "Arg";
    if (nuc == "AGC" || nuc == "AGT")
        return "Ser";
    if (nuc == "AGA" || nuc == "AGG")
        return "Arg";
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return "Gly";
    else
        return "Unrecognised codon";
}

inline int getCodonDistance(const std::string& refCdn, const std::string& altCdn) {
    int numDiffs = 0;
    assert(refCdn.length() == altCdn.length());
    for (std::string::size_type i = 0; i != refCdn.length(); i++) {
        if (refCdn[i] != altCdn[i]) {
            numDiffs = numDiffs + 1;
        }
    }
    return numDiffs;
}

double calculateNd(const std::string& refCdn, const std::string& altCdn, int diffNum);
double calculateN(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral);
double calculateNtS(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral);

inline void addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(const std::vector<string>& altCodons, const std::vector<int>& haveStop, const std::vector<int>& isUnrecognised, CDSComparisonMatrices* p) {

    for (std::vector<std::string>::size_type j = 0; j != altCodons.size() - 1; j++) {
        if (haveStop[j] == 1 || isUnrecognised[j] == 1)
            continue; // only consider this individual if it did not have a premature stop codon and there are no unrecognised letters
        for (std::vector<std::string>::size_type k = j+1; k != altCodons.size(); k++) {
            if (haveStop[k] == 1 || isUnrecognised[k] == 1)
                continue;
            
            double N_ijk; double S_ijk;
            double N_tS;
            if (altCodons[j] == altCodons[k]) {
                N_ijk = s_mapCodonToExpDist[altCodons[j]];
                N_tS = s_mapCodonToExpTs[altCodons[j]];
            } else {
                int d = getCodonDistance(altCodons[j],altCodons[k]);
                double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
                double s_d_ijk = d - n_d_ijk;
                p->N_d_jk[j][k] = p->N_d_jk[j][k] + n_d_ijk;
                p->S_d_jk[j][k] = p->S_d_jk[j][k] + s_d_ijk;

                N_ijk = calculateN(altCodons[j],altCodons[k], d, false);
                N_tS = calculateNtS(altCodons[j],altCodons[k], d, false);
            }
            S_ijk = (3 - N_ijk);
            p->N_jk[j][k] = p->N_jk[j][k] + N_ijk; p->S_jk[j][k] = p->S_jk[j][k] + S_ijk;
            p->tS_N_jk[j][k] = p->tS_N_jk[j][k] + N_tS;
            p->tS_S_jk[j][k] = p->tS_S_jk[j][k] + (1 - N_tS);
            p->tV_N_jk[j][k] = p->tV_N_jk[j][k] + (N_ijk - N_tS);
            p->tV_S_jk[j][k] = p->tV_S_jk[j][k] + (2 - (N_ijk - N_tS));
        }
    }
}


#endif /* UtilsCodingStats_hpp */
