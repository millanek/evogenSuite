//
//  UtilsIUPAC.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef UtilsIUPAC_hpp
#define UtilsIUPAC_hpp

#include <iostream>
#include "UtilsGeneral.hpp"

// IUPAC alphabet tools

// Complement a base using the full IUPAC alphabet (allows for lowercase bases)
char complementIUPAC(char c);
// Reverse complement a sequence using the full iupac alphabet (allows for lowercase bases)
std::string reverseComplementIUPAC(const std::string& seq);

// Return IUPAC nucleotide codes for heterozygous calls
std::string getAmbiguityCode(char b1, char b2);

// Given reference DNA base and IUPAC encoded het - returns DNA base for the alteranative allele
inline char disambiguateIUPAC(char refBase, char IUPACbase) {
    char altBase;
    switch (std::toupper(IUPACbase)) // Allows for lowercase bases
    {
        case 'K': altBase = (refBase == 'G') ? 'T' : 'G'; break;
        case 'M': altBase = (refBase == 'A') ? 'C' : 'A'; break;
        case 'R': altBase = (refBase == 'A') ? 'G' : 'A'; break;
        case 'S': altBase = (refBase == 'C') ? 'G' : 'C'; break;
        case 'W': altBase = (refBase == 'A') ? 'T' : 'A'; break;
        case 'Y': altBase = (refBase == 'C') ? 'T' : 'C'; break;
        default: assert(false);
    }
    return altBase;
}

// Returns both DNA bases for het ambiguous codes
inline std::string returnHetIUPAC(char IUPACbase) {
    char firstBase; char secondBase;
    switch (std::toupper(IUPACbase)) // Allows for lowercase bases
    {
        case 'K': firstBase = 'T'; secondBase = 'G'; break;
        case 'M': firstBase = 'C'; secondBase = 'A'; break;
        case 'R': firstBase = 'G'; secondBase = 'A'; break;
        case 'S': firstBase = 'G'; secondBase = 'C'; break;
        case 'W': firstBase = 'T'; secondBase = 'A'; break;
        case 'Y': firstBase = 'T'; secondBase = 'C'; break;
        default: assert(false);
    }
    std::string bases = ""; bases = firstBase + secondBase;
    return bases;
}

inline bool isDNAonly(char base) {
    bool DNAonly;
    switch (std::toupper(base))  // Allows for lowercase bases
    {
        case 'A': DNAonly = true; break;
        case 'C': DNAonly = true; break;
        case 'G': DNAonly = true; break;
        case 'T': DNAonly = true; break;
        default: DNAonly = false;
    }
    return DNAonly;
}

inline bool isDNAonlySeq(const std::string& seq) {
    bool DNAonly;
    for (std::string::size_type i = 0; i != seq.length(); i++) {
        DNAonly = isDNAonly(seq[i]);
        if (!DNAonly)
            return false;
    }
    return true;
}

#endif /* UtilsIUPAC_hpp */
