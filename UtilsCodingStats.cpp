//
//  UtilsCodingStats.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 28.06.23.
//

#include "UtilsCodingStats.hpp"

double calculateNd(const std::string& refCdn, const std::string& altCdn, int diffNum) {
    assert(refCdn.length() == altCdn.length());
    std::string refStep = refCdn;
    int countNS = 0;
    double Nd = 0;
    
    if (diffNum == 1) {
        if (s_mapCodonPairToSynonymous.count(refCdn+altCdn) == 0)
            Nd = 1.0;
    }
    
    if (diffNum == 2) {
        std::vector<std::string::size_type> diffPos;
        for (std::string::size_type i = 0; i != refCdn.length(); i++) {
            if (refCdn[i] != altCdn[i]) {
                diffPos.push_back(i);
            }
        }
        refStep[diffPos[0]] = altCdn[diffPos[0]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0)
            countNS++;
        refStep = refCdn; // the reverse order of mutations
        refStep[diffPos[1]] = altCdn[diffPos[1]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0)
            countNS++;
        Nd = countNS/2.0;
    }
    
    if (diffNum == 3) { // six different mutation pathways (orders of mutations)
        refStep[0] = altCdn[0]; // 1
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        string refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 2
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 3
        refStep[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 4
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 5
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 6
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        Nd = countNS/6.0;
    }
    assert(Nd <= diffNum);
    return Nd;
}


double calculateN(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral) {
    assert(refCdn.length() == altCdn.length()); 
    assert(diffNum >= 0 && diffNum <= 3);
    double N = 0;
    
    if (diffNum == 0) {
        return s_mapCodonToExpDist[refCdn];
    }
    
    if (diffNum == 1) {
        if (refAncestral) {
            return s_mapCodonToExpDist[refCdn];
        } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
            return (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[altCdn])/2;
        }
    } else {
        if (diffNum == 2) {
            // Find where the diffs are
            std::vector<std::string::size_type> diffPos;
            for (std::string::size_type i = 0; i != refCdn.length(); i++) {
                if (refCdn[i] != altCdn[i]) {
                    diffPos.push_back(i);
                }
            }
            
            // Then calculate N for the possible mutation paths:
            double Nsum = 0;
            // e.g. TAA -> TGA -> TGG
            std::string stepCdn = refCdn;
            stepCdn[diffPos[0]] = altCdn[diffPos[0]];
            Nsum = (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn])/2;
            
            // e.g. TAA -> TAG -> TGG
            stepCdn = refCdn;
            stepCdn[diffPos[1]] = altCdn[diffPos[1]];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn])/2);
            
            if (refAncestral) {
                return Nsum / 2;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. TAA <- TGA <- TGG
                stepCdn = altCdn;
                stepCdn[diffPos[0]] = refCdn[diffPos[0]];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn])/2);
                
                // e.g. TAA <- TAG <- TGG
                stepCdn = altCdn; // the reverse order of mutations
                stepCdn[diffPos[1]] = refCdn[diffPos[1]];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn])/2);
                // Get the average N for the four mutation paths:
                return Nsum / 4;
            }
        }
        
        // This could surely be simplified but I write out all the possible mutation paths explicitly
        // one could even have a lookup table for all the three letter pairs and what N scores they give
        if (diffNum == 3) {
            double Nsum = 0;
            
            // e.g. AAA -> TAA -> TGA -> TGG
            std::string stepCdn = refCdn; stepCdn[0] = altCdn[0];
            std::string step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3;
            
            // e.g. AAA -> TAA -> TAG -> TGG
            stepCdn = refCdn; stepCdn[0] = altCdn[0];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> TGA -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> AGG -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> AGG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> TAG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            if (refAncestral) {
                return Nsum / 6;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. AAA <- TAA <- TGA <- TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- TAA <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA -> AGA -> TGA -> TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AGA <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // Get the average N for the twelve mutation paths:
                return Nsum / 12;
            }
        }
        
    }
    return N;
}

double calculateNtS(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral) {
    assert(refCdn.length() == altCdn.length());
    assert(diffNum >= 0 && diffNum <= 3);
    double NtS = 0;
    
    if (diffNum == 0) {
        return s_mapCodonToExpTs[refCdn];
    }
    
    if (diffNum == 1) {
        if (refAncestral) {
            return s_mapCodonToExpTs[refCdn];
        } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
            return (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[altCdn])/2;
        }
    } else {
        if (diffNum == 2) {
            // Find where the diffs are
            std::vector<std::string::size_type> diffPos;
            for (std::string::size_type i = 0; i != refCdn.length(); i++) {
                if (refCdn[i] != altCdn[i]) {
                    diffPos.push_back(i);
                }
            }
            
            // Then calculate N for the possible mutation paths:
            double Nsum = 0;
            // e.g. TAA -> TGA -> TGG
            std::string stepCdn = refCdn;
            stepCdn[diffPos[0]] = altCdn[diffPos[0]];
            Nsum = (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn])/2;
            
            // e.g. TAA -> TAG -> TGG
            stepCdn = refCdn;
            stepCdn[diffPos[1]] = altCdn[diffPos[1]];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn])/2);
            
            if (refAncestral) {
                return Nsum / 2;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. TAA <- TGA <- TGG
                stepCdn = altCdn;
                stepCdn[diffPos[0]] = refCdn[diffPos[0]];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn])/2);
                
                // e.g. TAA <- TAG <- TGG
                stepCdn = altCdn; // the reverse order of mutations
                stepCdn[diffPos[1]] = refCdn[diffPos[1]];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn])/2);
                // Get the average N for the four mutation paths:
                return Nsum / 4;
            }
        }
        
        // This could surely be simplified but I write out all the possible mutation paths explicitly
        // one could even have a lookup table for all the three letter pairs and what N scores they give
        if (diffNum == 3) {
            double Nsum = 0;
            
            // e.g. AAA -> TAA -> TGA -> TGG
            std::string stepCdn = refCdn; stepCdn[0] = altCdn[0];
            std::string step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3;
            
            // e.g. AAA -> TAA -> TAG -> TGG
            stepCdn = refCdn; stepCdn[0] = altCdn[0];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> TGA -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> AGG -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> AGG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> TAG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            if (refAncestral) {
                return Nsum / 6;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. AAA <- TAA <- TGA <- TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- TAA <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA -> AGA -> TGA -> TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AGA <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // Get the average N for the twelve mutation paths:
                return Nsum / 12;
            }
        }
        
    }
    return NtS;
}
