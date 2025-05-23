//
//  UtilsAnnotation.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef UtilsAnnotation_hpp
#define UtilsAnnotation_hpp

#include <stdio.h>
#include "UtilsGeneral.hpp"
#include "UtilsIUPAC.hpp"


inline std::string geneFromTranscript(const std::string& transcript) {
    std::vector<std::string> transcriptVec = split(transcript,'.');
    std::string geneName = transcriptVec[0] + "." + transcriptVec[1] + "." + transcriptVec[2] + "." + transcriptVec[3];
    return geneName;
}


class SimpleCoordsBed {
public:
    SimpleCoordsBed() {};

    SimpleCoordsBed(std::ifstream*& bedFile, std::map<int, string>& linearToGenomeCoordMap) {
        loadSimpleCoords(bedFile, linearToGenomeCoordMap);
    }
    
private:
    // Require: coordinates are sorted by scaffold
    std::map<string, std::vector<std::vector<string> > > loadSimpleCoords(std::ifstream*& bedFile, std::map<int,string>& linearToGenomeCoordMap) {
        std::map<string, std::vector<std::vector<string> > > coordsScaffoldMap;
        std::vector<std::vector<string> > coordsInScaffold;
        int linearPosition = 0; // Linear positions will be 0-indexed (will start at zero)
        string line;
        // Do the first line to find the name of the first scaffold
        getline(*bedFile, line);
        std::vector<string> bedVector = split(line, '\t');
        string currentScaffold = bedVector[0];
        int left = atoi(bedVector[1].c_str());
        int right = atoi(bedVector[2].c_str());
        int intervalLength = right - left;
        coordsInScaffold.push_back(bedVector);
        
        for (int i = 0; i < intervalLength; i++) {
            // Mapping will be 1-indexed, as in VCF files
            linearToGenomeCoordMap[linearPosition+i] = currentScaffold+"\t"+numToString(left+i+1);
        } linearPosition = linearPosition + intervalLength;
        
        //
        while (getline(*bedFile, line)) {
            std::vector<string> bedVector = split(line, '\t');
            string scaffold = bedVector[0];
            int left = atoi(bedVector[1].c_str());
            int right = atoi(bedVector[2].c_str());
            int intervalLength = right - left;
            
            for (int i = 0; i < intervalLength; i++) {
                // Mapping will be 1-indexed, as in VCF files
                int pos = linearPosition+i;
                if (pos % 1000000 == 0)
                    std::cerr << "Loaded and mapped " << pos/1000000 << "Mb" << std::endl;
                linearToGenomeCoordMap[pos] = currentScaffold+"\t"+numToString(left+i+1);
            } linearPosition = linearPosition + intervalLength;
            
            if (scaffold == currentScaffold) {
                coordsInScaffold.push_back(bedVector);
            } else {
                coordsScaffoldMap[currentScaffold] = coordsInScaffold;
                coordsInScaffold.clear();
                currentScaffold = scaffold;
            }
        }
        return coordsScaffoldMap;
    }
};

class LinkedCoordsBed {
public:
    LinkedCoordsBed() {};
    
    LinkedCoordsBed(std::ifstream*& bedFile) {
        elementCoordsVector = loadLinkedCoords(bedFile);
    }
    
    std::vector<std::vector<std::vector<string> > > elementCoordsVector;
    
    
    std::vector<std::vector<string> > getElementOuterBoundaries() {
        std::vector<std::vector<string> > allOuterBounds;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            //if (thisElement[0].size() < 4) {
            //    std::cerr << "thisElement[0].size()" << thisElement.size() << std::endl;
            //}
            std::vector<string> elementOuterBounds;
            elementOuterBounds.push_back(thisElement[0][0]);
            elementOuterBounds.push_back(thisElement[0][1]);
            elementOuterBounds.push_back(thisElement.back()[2]);
            allOuterBounds.push_back(elementOuterBounds);
        }
        return allOuterBounds;
    }
    
    std::vector<string> getElementNames() {
        std::vector<string> elementNames;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            //if (thisElement[0].size() < 4) {
            //    std::cerr << "thisElement[0].size()" << thisElement.size() << std::endl;
            //}
            elementNames.push_back(thisElement[0][3]);
        }
        return elementNames;
    }
    
    std::vector<double> getMeanPerElement(std::unordered_map<string, double>& posDxyMap) {
        std::vector<double> elementDxyValues;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            double elementDxyTotal = 0; int elementLengthTotal = 0;
            for (int j = 0; j < (int)thisElement.size(); j++) {
                std::vector<string> thisCoordVector = thisElement[j];
                string currentScaffold = thisCoordVector[0];
                int left = atoi(thisCoordVector[1].c_str());
                int right = atoi(thisCoordVector[2].c_str());
                int intervalLength = right - left;
                for (int k = 0; k < intervalLength; k++) {
                    string genomeLocation = currentScaffold+"\t"+numToString(left+k+1);
                    if (posDxyMap.count(genomeLocation) == 1) {
                        elementDxyTotal = elementDxyTotal + posDxyMap[genomeLocation];
                    }
                }
                elementLengthTotal = elementLengthTotal + intervalLength;
            }
            elementDxyValues.push_back(elementDxyTotal/elementLengthTotal);
        }
        return elementDxyValues;
    }
    
private:
    // Require: coordinates are sorted by scaffold
    std::vector<std::vector<std::vector<string> > > loadLinkedCoords(std::ifstream*& bedFile) {
        std::vector<std::vector<std::vector<string> > > elementCoordsVector;
        std::vector<std::vector<string> > thisElement;
        
        // Do the first line to find the name of the first scaffold and of the first element
        string line;
        getline(*bedFile, line);
        std::vector<string> bedVector = split(line, '\t');
        string currentElementName = bedVector[3];
        thisElement.push_back(bedVector);
        
        //
        while (getline(*bedFile, line)) {
            std::vector<string> bedVector = split(line, '\t');
            string element = bedVector[3];
            if (element == currentElementName) {
                thisElement.push_back(bedVector);
            } else {
                elementCoordsVector.push_back(thisElement);
                thisElement.clear();
                currentElementName = element;
                thisElement.push_back(bedVector);
            }
        }
        return elementCoordsVector;
    }
};



    
class Annotation {
public:
    Annotation() {};
    
    Annotation(const string geneFileName, bool includePartial) {
        if(!geneFileName.empty()) {
            std::ifstream* annotFile = new std::ifstream(geneFileName.c_str());
            annotationMap = loadAnnotationMap(annotFile, includePartial);
            getWgGeneTranscriptCounts();
            annotateGeneStartsEnds();
            initialised = true;
            
        }
    }
    
    std::map<std::string, std::vector<std::vector<std::string> > > annotationMap;
    std::map<string, int> geneTranscriptCounts;  // How many transcripts has each gene
    bool initialised = false;
    
    // For tracking as we go though the VCF
    string previousGene = "";
    string currentGene = "";
    string SNPcategory = "";
    bool bInGene = false;
    
    string getCategoryOfSNP(const string& SNPscaffold, const string& SNPlocus) {
        string SNPcategory = "other non-coding";
        std::vector<string> scaffoldTranscriptStartEnd = transcriptStartEndMap[SNPscaffold];
        string inGene = "";
        for (std::vector<std::vector<string> >::size_type i = 0; i != scaffoldTranscriptStartEnd.size(); i++) {
            std::vector<string> startEndVec = split(scaffoldTranscriptStartEnd[i], '\t');
            if (SNPlocus >= startEndVec[1] && SNPlocus <= startEndVec[2]) {
                string thisTranscript = startEndVec[0];
                int numDots = (int)std::count(thisTranscript.begin(), thisTranscript.end(), '.');
                SNPcategory = "intron";
                if (numDots == 4)  inGene = geneFromTranscript(thisTranscript);
                else inGene = thisTranscript;
                std::vector<string> exons = annotationMapTranscriptMap[SNPscaffold][thisTranscript];
                for (std::vector<string>::size_type j = 0; j != exons.size(); j++) {
                    std::vector<string> exonVec = split(exons[j], '\t');
                    if (SNPlocus >= exonVec[1] && SNPlocus <= exonVec[2]) {
                        SNPcategory = "exon";
                    }
                    break;
                }
                if (SNPcategory == "exon") { break; }
            } //else if (SNPcategory == "intron") { if (inGene != geneFromTranscript(startEndVec[0])) { break; } }
        }
        return SNPcategory;
    }
    
    void getSNPgeneDetails(const string& SNPscaffold, const int SNPlocus);
    
    bool bUpdateGene() {
        if (currentGene != "" && previousGene != "" && currentGene != previousGene) return true;
        else return false;
    }
    
    
    int getTranscriptCount(const std::string& geneOrTranscript) {
        int numTranscripts = 0;
        int numDots = (int)std::count(geneOrTranscript.begin(), geneOrTranscript.end(), '.');
        if (numDots < 4)
            numTranscripts = geneTranscriptCounts[geneOrTranscript];
        else if (numDots == 4) {
            //std::vector<std::string> transcriptVec = split(geneOrTranscript,'.');
            //std::string geneName = transcriptVec[0] + "." + transcriptVec[1] + "." + transcriptVec[2] + "." + transcriptVec[3];
            std::string geneName = geneFromTranscript(geneOrTranscript);
            numTranscripts = geneTranscriptCounts[geneName];
        } else {
            std::cerr << "NumDots: " << numDots << std::endl;
            std::cerr << "geneOrTranscript: " << geneOrTranscript << std::endl;
            assert(false);
        }
        return numTranscripts;
    }
    
private:
    std::map<std::string, std::vector<std::string> > transcriptStartEndMap; // Start and end positions of each transcript
    std::map<std::string, std::map<std::string, std::vector<std::string> > > annotationMapTranscriptMap; // Quickly find the exons of any transcript
    
    
    std::string getTranscriptName(const string& geneColumn, bool& bPartial) {
        string thisTranscriptName;
        std::vector<string> findIfPartial = split(geneColumn, '-');
        if (findIfPartial.size() == 2) {
            bPartial = true;
            thisTranscriptName = findIfPartial[1];
        } else if (findIfPartial.size() == 1) {
            thisTranscriptName = geneColumn;
        }
        return thisTranscriptName;
    }
    
    // Load up the annotation file
    std::map<string, std::vector<std::vector<string> > > loadAnnotationMap(std::ifstream*& geneFile, const bool bUsePartial) {
        std::map<string, std::vector<std::vector<string> > > annotationMap;
        std::vector<std::vector<string> > annotation;
        string line;
        std::vector<string> currentTranscript;
        getline(*geneFile, line);
        currentTranscript.push_back(line);
        std::vector<string> annotLineVec = split(line, '\t');
        std::vector<string> findIfPartial = split(annotLineVec[4], '-');
        string currentScaffold = annotLineVec[0];
        bool bThisTranscriptPartial; string currentTranscriptName = getTranscriptName(annotLineVec[4], bThisTranscriptPartial);
        while (getline(*geneFile, line)) {
            std::vector<string> annotLineVec = split(line, '\t');
            bool bCurrentLinePartial = false; string currentLineTranscriptName = getTranscriptName(annotLineVec[4], bCurrentLinePartial);
            if (bCurrentLinePartial)
                bThisTranscriptPartial = true;
            if (annotLineVec[0] == currentScaffold) {
                if (currentLineTranscriptName == currentTranscriptName) {
                    currentTranscript.push_back(line);
                } else {
                    if (!bThisTranscriptPartial || bUsePartial)
                        annotation.push_back(currentTranscript);
                    currentTranscript.clear();
                    currentTranscript.push_back(line);
                    currentTranscriptName = currentLineTranscriptName;
                    bThisTranscriptPartial = bCurrentLinePartial;
                }
            } else {
                if (!bThisTranscriptPartial || bUsePartial)
                    annotation.push_back(currentTranscript);
                annotationMap[currentScaffold] = annotation;
                annotation.clear();
                currentTranscript.clear();
                currentTranscript.push_back(line);
                currentScaffold = annotLineVec[0];
                currentTranscriptName = currentLineTranscriptName;
                bThisTranscriptPartial = bCurrentLinePartial;
            }
        }
        return annotationMap;
    }
    
    void getAnnotationPerGeneDetails(const std::vector<std::vector<string> >& scaffoldAnnotation) {
        int thisGeneTranscriptCount;
        string previousGeneName = "";
        for (std::vector<string>::size_type i = 0; i != scaffoldAnnotation.size(); i++) {
            std::vector<string> annotLineVec = split(scaffoldAnnotation[i][0], '\t');
            int numDots = (int)std::count(annotLineVec[4].begin(), annotLineVec[4].end(), '.');
            std::string geneName;
            if (numDots == 4)  geneName = geneFromTranscript(annotLineVec[4]);
            else geneName = annotLineVec[4];
            if (previousGeneName == geneName) {
                thisGeneTranscriptCount++;
            } else {
                geneTranscriptCounts[previousGeneName] = thisGeneTranscriptCount;
                thisGeneTranscriptCount = 1;
                previousGeneName = geneName;
            }
            if (i == (scaffoldAnnotation.size()-1)) {   // Add the final gene in the scaffold
                geneTranscriptCounts[previousGeneName] = thisGeneTranscriptCount;
            }
        }
    }
    
    
    void getWgGeneTranscriptCounts() {
        for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator it = annotationMap.begin(); it != annotationMap.end(); it++) {
            getAnnotationPerGeneDetails(it->second);
        }
    }
    
    void annotateGeneStartsEnds();
};


// Get the sequence for a particular individual, for a given region of the scaffold, as defined by the annotation vector
std::string getIndividualSequenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentIndividualWholeScaffoldSequence);
// Get the reference sequence for a given region of the scaffold defined by the annotation vector
std::string getReferenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentScaffoldReference);





class BedCoordinateFeatures {
public:
    BedCoordinateFeatures() : initialised(false) {};
    
    BedCoordinateFeatures(std::ifstream*& bedFile, int numValues = 0, bool leftCoordOneIndexed = false, bool debug = false) {
        loadBedFeatureMap(bedFile, leftCoordOneIndexed,numValues, debug);
        initialised = true;
    }
    
    bool initialised;
    map<string, vector< vector<int> > > BedFeatureMap;
    vector< map<std::string, vector<double>> > BedFeatureValueMaps;

    
    // Count how much of a genomic interval is covered by the bed features (returns the number of basepairs)
    // Left coordinates of features are 0-indexed (bed file)
    // Query is 1-indexed
    // query:              -----------------------
    // feature: 1)  ---                                   (f[1] <= start)
    //          2)                                  ---   (f[0] >= end)
    //          3)     ------                             (f[0] < start && f[1] <= end)
    //          4)              --------                  (f[0] >=start && f[1] <= end)
    //          5)                             -------    (f[0] >=start && f[1] > end)
    //          6)       -----------------------------    (f[0] < start && f[1] > end)
    int getNumBPinRegion(const string& scaffold, const int start, const int end) const {
        assert(start < end);
        try {
            std::vector<std::vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            /*for (std::map<std::string, std::vector<std::vector<int> > >::iterator i = accessibleGenomeMap.begin(); i != accessibleGenomeMap.end(); i++) {
             std::cerr << "There is scaffold: " << i->first << " in the map" << std::endl;
             } */
            
            
            // Binary search to find the first feature whose end coordinate is greater
            // or equal to the start of the region in question
            vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),start);
            int numBP = 0;
            if (itStart != featuresThisSc[1].end()) {  // if (start < f[1])    ---  excludind case 1)
                vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
                
                // Sum the lengths
                while (featuresThisSc[0][index] <= end && index < featuresThisSc[0].size()) { // if (f[0] >= end)   ---    ecluding case 2)
                    //std::cerr << "There are " << featuresThisSc[0].size() << " intervals" << std::endl;
                    //std::cerr << "Starts: " << std::endl;
                    //print_vector(featuresThisSc[0], std::cerr);
                    //std::cerr << "Ends: " << std::endl;
                    //print_vector(featuresThisSc[1], std::cerr);
                    // Now we know the overlap is > 0
                    //std::cerr << "start, end, featuresThisSc[0][index], featuresThisSc[1][index]: " << start << ", " << end << ", " << featuresThisSc[0][index] << ", " << featuresThisSc[1][index] << std::endl;
                    
                    if (featuresThisSc[0][index] < start && featuresThisSc[1][index] <= end)
                        numBP = numBP + (featuresThisSc[1][index] - start) + 1;
                    else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] <= end)
                        numBP = numBP + (featuresThisSc[1][index] - featuresThisSc[0][index]);
                    else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] > end)
                        numBP = numBP + (end - featuresThisSc[0][index]);
                    else if (featuresThisSc[0][index] < start && featuresThisSc[1][index] > end)
                        numBP = numBP + (end - start) + 1;
                    index++;
                }
            }
            return numBP;
        } catch (const std::out_of_range& oor) {
            //std::cerr << "No features on scaffold: " << scaffold << std::endl;
            return 0;
        }
    }
    
    std::vector<double> getMeanValuesInRegion(const string& scaffold, const int start, const int end) const {
        assert(start < end);
        try {
            vector<vector<int>> featuresThisSc = BedFeatureMap.at(scaffold);
            
            vector<vector<double>> valuesThisScaffold;
            for(int i = 0; i < BedFeatureValueMaps.size(); i++) {
                std::vector<double> thisValuesThisScaffold = BedFeatureValueMaps[i].at(scaffold);
                valuesThisScaffold.push_back(thisValuesThisScaffold);
            }
           // if (start > 116000) { exit (1); }
                                                                         
          
           // std::cerr << "Start: " << start << "; End: " << end << std::endl;
            //print_vector(featuresThisSc[0], std::cerr);
            //std::cerr << "Ends: " << std::endl;
            
            // Binary search to find the first feature whose end coordinate is greater
            // or equal to the start of the region in question
            vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),start);
           // std::cerr << "itStart: " << *itStart << std::endl;
            int numBPtotal = 0; int numBPthisFeature = 0;
            vector<double> sumPerBPvalues(BedFeatureValueMaps.size(),0);
            vector<double> meanValues(BedFeatureValueMaps.size(),NAN);
            if (itStart != featuresThisSc[1].end()) {  // if (start < f[1])    ---  excludind case 1)
                vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
                
               // std::cerr << "index\t" << index << std::endl;
                // Sum the lengths
                while (featuresThisSc[0][index] <= end && index < featuresThisSc[0].size()) { // if (f[0] >= end)   ---    ecluding case 2)
                    //std::cerr << "featuresThisSc[0][index]\t" << featuresThisSc[0][index] << std::endl;
                    //std::cerr << "BedFeatureValueMaps.size()\t" << BedFeatureValueMaps.size() << std::endl;
                    vector<double> valuesThisFeature(BedFeatureValueMaps.size(),0);
                    for(int i = 0; i < BedFeatureValueMaps.size(); i++) {
                        valuesThisFeature[i] = valuesThisScaffold[i][index];
                    }
                  //  std::cerr << "valuesThisFeature[0]\t" << valuesThisFeature[0] << std::endl;
                    if (featuresThisSc[0][index] < start && featuresThisSc[1][index] <= end) {
                        numBPthisFeature = (featuresThisSc[1][index] - start) + 1;
                      //  std::cerr << featuresThisSc[1][index] << " - " << start << " + 1" << std::endl;
                    } else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] <= end) {
                        numBPthisFeature = (featuresThisSc[1][index] - featuresThisSc[0][index]);
                      //  std::cerr << featuresThisSc[1][index] << " - " << featuresThisSc[0][index] << std::endl;
                    } else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] > end) {
                        numBPthisFeature = (end - featuresThisSc[0][index]);
                      //  std::cerr << end << " - " << featuresThisSc[0][index] << std::endl;
                    } else if (featuresThisSc[0][index] < start && featuresThisSc[1][index] > end) {
                        numBPthisFeature = (end - start) + 1;
                      //  std::cerr << end << " - " << start << " + 1" << std::endl;
                    }
                    
                   // std::cerr << "numBPthisFeature: " << numBPthisFeature << std::endl;
                   // std::cerr << "valuesThisFeature[0]: " << valuesThisFeature[0] << std::endl;
                    index++;
                    for(int i = 0; i < BedFeatureValueMaps.size(); i++) {
                        sumPerBPvalues[i] += (valuesThisFeature[i] * numBPthisFeature);
                    }
                    numBPtotal += numBPthisFeature;
                }
                // std::cerr << "sumPerBPvalues[0]: " << sumPerBPvalues[0] << std::endl;
                // std::cerr << "numBPtotal: " << numBPtotal << std::endl;
                for(int i = 0; i < BedFeatureValueMaps.size(); i++) {
                    meanValues[i] = (double)sumPerBPvalues[i]/numBPtotal;
                }
                // std::cerr << "meanValue: " << meanValues[0] << std::endl;
            }
            return meanValues;
        } catch (const std::out_of_range& oor) {
            //std::cerr << "No features on scaffold: " << scaffold << std::endl;
            vector<double> meanValues(NAN,BedFeatureValueMaps.size());
            return meanValues;
        }
    }
    
    
    vector< vector<string> > getFeaturesinRegion(const string& scaffold, const int start, const int end) {
        vector< vector<string> > allFeatures;
        try {
            vector< vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            //std::cerr << "There are " << aGThisSc[0].size() << " accessible intervals" << std::endl;
            /*for (std::map<std::string, std::vector<std::vector<int> > >::iterator i = accessibleGenomeMap.begin(); i != accessibleGenomeMap.end(); i++) {
             std::cerr << "There is scaffold: " << i->first << " in the map" << std::endl;
             } */
            
            
            // Binary search to find the first feature whose end coordinate is greater
            // or equal to the start of the region in question
            vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),start);
            vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
            // Now check that the start coordinate of the feature is within
            
            // Sum the lengths
            while (featuresThisSc[0][index] <= end) {
                vector<string> feature;
                feature.push_back(scaffold); feature.push_back(numToString(featuresThisSc[0][index]));
                feature.push_back(numToString(featuresThisSc[1][index]));
                index++;
            }
            return allFeatures;
        } catch (const std::out_of_range& oor) {
            std::cerr << "No features on scaffold: " << scaffold << std::endl;
            return allFeatures;
        }
    }
    
    
    
    // In the bed file, the left coordinate is zero indexed
    // bp is a 1-indexed coordinate
    bool findIfBPinBedFile(const string& scaffold, const int bp) {
        try {
            std::vector<std::vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            
            // Binary search to find the first element in the accessible genome annotation whose end coordinate is greater
            // or equal to the basepair in question
            std::vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),bp);
            std::vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
            
            // Then if the start coordinate is less than the coordinate of the basepair, the bp is covered by the bed file features
            if (featuresThisSc[0][index] < bp) {
                return true;
            } else {
                return false;
            }
        } catch (const std::out_of_range& oor) {
            return false;
        }
    }
    
    // Load up the file specifying the genomic coordinates (the bed file needs to be sorted by chromosome)
    void loadBedFeatureMap(std::ifstream*& bedFile, const int numValues, const bool leftCoordOneIndexed, const bool debug) {
        std::vector<std::vector<int> > BedFeaturesThisScaffold;
        std::vector<int> featureStarts;
        std::vector<int> featureEnds;
        std::vector <std::vector<double>> featureValues;
        string line;
        string previousScaffold = "";
        
        for(int i = 0; i < numValues; i++) {
            std::map<std::string, std::vector<double> > BedFeatureValueMapInit;
            BedFeatureValueMaps.push_back(BedFeatureValueMapInit);
            std::vector<double> featureValuesInit;
            featureValues.push_back(featureValuesInit);
        }
        
        while (getline(*bedFile, line)) {
            std::vector<string> currentFeature = split(line, '\t');
            
           // std::cerr << "line: " << line << std::endl;
           // std::cerr << "currentFeature.size(): " << currentFeature.size() << std::endl;
            
            if (currentFeature.size() != (3 + numValues)) {
                std::cerr << "ERROR: The file must have " << (3 + numValues) << "columns\n";
                exit(EXIT_FAILURE);
            }
            
            if (currentFeature[0] != previousScaffold && previousScaffold != "") {
                //std::cerr << "Loading scaffold: " << previousScaffold << std::endl;
                if (debug) {
                    std::cerr << "Loading scaffold: " << previousScaffold << std::endl;
                    std::cerr << "Starts: " << std::endl;
                    print_vector(featureStarts, std::cerr);
                    std::cerr << "Ends: " << std::endl;
                    print_vector(featureEnds, std::cerr);
                }
                BedFeaturesThisScaffold.push_back(featureStarts);
                BedFeaturesThisScaffold.push_back(featureEnds);
                BedFeatureMap[previousScaffold] = BedFeaturesThisScaffold;
                BedFeaturesThisScaffold.clear(); featureStarts.clear(); featureEnds.clear();
                for(int i = 0; i < numValues; i++) {
                    BedFeatureValueMaps[i][previousScaffold] = featureValues[i];
                    featureValues[i].clear();
                }
                //std::cerr << "Loading scaffold: " << currentFeature[0] << std::endl;
            }
            int subtractFromLeft = (leftCoordOneIndexed) ? 1 : 0;
            featureStarts.push_back(atoi(currentFeature[1].c_str()) - subtractFromLeft);
            featureEnds.push_back(atoi(currentFeature[2].c_str()));
            for(int i = 0; i < numValues; i++) { featureValues[i].push_back(stringToDouble(currentFeature[3+i])); }
            //std::cerr << "Loading scaffold: " << currentFeature[0] << std::endl;
            previousScaffold = currentFeature[0];
        }
        // Final line / final scaffold
        BedFeaturesThisScaffold.push_back(featureStarts);
        BedFeaturesThisScaffold.push_back(featureEnds);
        BedFeatureMap[previousScaffold] = BedFeaturesThisScaffold;
        // std::cerr << "Loading scaffold: " << previousScaffold << std::endl;
        for(int i = 0; i < numValues; i++) { BedFeatureValueMaps[i][previousScaffold] = featureValues[i]; }
    }
    
};



class AccessibleGenome : public BedCoordinateFeatures {
public:
    AccessibleGenome() {};
    
    AccessibleGenome(const string& bedFileName) {
        if (!bedFileName.empty()) {
            std::ifstream* accessibleGenomeBedFile = new std::ifstream(bedFileName);
            assertFileOpen(*accessibleGenomeBedFile, bedFileName);
            std::cerr << std::endl;
            std::cerr << "Loading the accessible genome annotation from the file: " << bedFileName << std::endl;
            loadBedFeatureMap(accessibleGenomeBedFile, false, false, false);
            std::cerr << "Done" << std::endl;
            initialised = true;
        }
    }
    
    int getAccessibleBPinRegion(const string& scaffold, const int start, const int end) const {
        return getNumBPinRegion(scaffold, start, end);
    }
    
    // Accessible genome is a bed file (left coordinate is zero indexed)
    // bp is a 1-indexed coordinate
    bool findIfBPaccessible(const string& scaffold, const int bp) {
        return findIfBPinBedFile(scaffold,bp);
    }
    
    string getAccessibleSeqForScaffold(const string& scaffold, const string& fullString) {
        std::vector<std::vector<int> > aGThisSc = BedFeatureMap[scaffold];
        
        string accessibleString; accessibleString.reserve(fullString.length()); accessibleString = "";
        for (std::vector<int>::size_type i = 0; i < aGThisSc[0].size(); i++) {
            accessibleString.append(fullString.substr(aGThisSc[0][i],aGThisSc[1][i]- aGThisSc[0][i]));
        }
        return accessibleString;
    }

};

class RecombinationMap : public BedCoordinateFeatures {
public:
    RecombinationMap() {};
    
    RecombinationMap(const string& bedFileName) {
        if (!bedFileName.empty()) {
            std::ifstream* RecombinationMapFile = new std::ifstream(bedFileName);
            std::cerr << std::endl;
            std::cerr << "Loading the recombination map" << std::endl;
            const int numValues = 1; const bool leftCoordOneIndexed = true;
            loadBedFeatureMap(RecombinationMapFile, numValues, leftCoordOneIndexed, false);
            std::cerr << "Done" << std::endl;
            initialised = true;
        }
    }
    
    double getMeanRecombinationRate(const string& scaffold, const int start, const int end) const {
        return getMeanValuesInRegion(scaffold, start, end)[0];
                                                                        
    }

};

class IntervalFile : public BedCoordinateFeatures {
public:
    IntervalFile() {};
    
    IntervalFile(const string& bedFileName, const bool oneIndexed = false) {
        if (!bedFileName.empty()) {
            std::ifstream* IntervalFileHandle = new std::ifstream(bedFileName);
            std::cerr << std::endl;
            std::cerr << "Loading the input file: " << bedFileName << std::endl;
            const int numValues = 1;
            loadBedFeatureMap(IntervalFileHandle, numValues, oneIndexed, false);
            std::cerr << "Done" << std::endl;
            initialised = true;
        }
    }
    
    std::vector<double> getMeanValuesForRegion(const string& scaffold, const int start, const int end) const {
        return getMeanValuesInRegion(scaffold, start, end);
                                                                        
    }

};





#endif /* UtilsAnnotation_hpp */
