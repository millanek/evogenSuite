//
//  UtilsAnnotation.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "UtilsAnnotation.hpp"


string getIndividualSequenceForThisRegion(const vector<string>& thisRegionAnnotation, const string& strand, const string& currentIndividualWholeScaffoldSequence) {
    string geneSeq = "";
    for (std::vector<std::string>::size_type j = 0; j != thisRegionAnnotation.size(); j++) {
        std::vector<string> annotLineVec = split(thisRegionAnnotation[j], '\t');
        geneSeq = geneSeq + currentIndividualWholeScaffoldSequence.substr(atoi(annotLineVec[1].c_str())-1,atoi(annotLineVec[2].c_str())-atoi(annotLineVec[1].c_str()) + 1);
    }
    if (strand == "-")
        geneSeq = reverseComplementIUPAC(geneSeq);
    
    return geneSeq;
}


string getReferenceForThisRegion(const vector<string>& thisRegionAnnotation, const string& strand, const string& currentScaffoldReference) {
    // std::cerr << "generating gene refseq" << std::endl;
    string refSeq = "";
    for (vector<string>::size_type j = 0; j != thisRegionAnnotation.size(); j++) {
        vector<string> annotLineVec = split(thisRegionAnnotation[j], '\t');
        refSeq = refSeq + currentScaffoldReference.substr(atoi(annotLineVec[1].c_str())-1,atoi(annotLineVec[2].c_str())-atoi(annotLineVec[1].c_str()) + 1);
    }
    if (strand == "-")
        refSeq = reverseComplementIUPAC(refSeq);
    // std::cerr << "Done" << std::endl;
    return refSeq;
}


void Annotation::getSNPgeneDetails(const string& chr, const int SNPlocus) {
    vector<string> scaffoldTranscriptStartEnd = transcriptStartEndMap[chr];
    
    string tempCurrentGene = ""; string tempSNPcategory = "nonCoding";
    for (vector< vector<string> >::size_type i = 0; i != scaffoldTranscriptStartEnd.size(); i++) {
        vector<string> startEndVec = split(scaffoldTranscriptStartEnd[i], '\t');
        string thisTranscript = startEndVec[0];
        int geneStart = atoi(startEndVec[1].c_str()); int geneEnd = atoi(startEndVec[2].c_str()); string strand = startEndVec[3];
        //if (SNPlocus == 20001) { print_vector_stream(startEndVec, std::cerr);}
        if (strand == "+") {
            if (SNPlocus >= (geneStart-3000) && SNPlocus < geneStart) {
                tempSNPcategory = "promoter"; tempCurrentGene = thisTranscript;
                break;
            }
        } else if (strand == "-") {
            if (SNPlocus > geneEnd && SNPlocus <= (geneEnd+3000)) {
                tempSNPcategory = "promoter"; tempCurrentGene = thisTranscript;
                break;
            }
        }
        if (SNPlocus >= geneStart && SNPlocus <= geneEnd) {
            //int numDots = (int)std::count(thisTranscript.begin(), thisTranscript.end(), '.');
            tempSNPcategory = "intron";
            //if (numDots == 4)  inGene = geneFromTranscript(thisTranscript);
            tempCurrentGene = thisTranscript;
            std::vector<string> exons = annotationMapTranscriptMap[chr][thisTranscript];
            // std::cerr << exons.size() << std::endl;
            for (std::vector<string>::size_type j = 0; j != exons.size(); j++) {
                std::vector<string> exonVec = split(exons[j], '\t');
                //if (SNPlocus == 20001) { print_vector_stream(exonVec, std::cerr); }
                if (SNPlocus >= atoi(exonVec[1].c_str()) && SNPlocus <= atoi(exonVec[2].c_str())) {
                    tempSNPcategory = "exon";
                    break;
                }
            }
            break;
        }
        //else if (SNPcategory == "intron") { if (inGene != geneFromTranscript(startEndVec[0])) { break; } }
    }
    
    currentGene = tempCurrentGene;
    SNPcategory = tempSNPcategory;
    
    // To set the first "previous gene"
    if (previousGene == "" && currentGene != "") previousGene = currentGene;
}



void Annotation::annotateGeneStartsEnds() {
    for (std::map<std::string, std::vector<std::vector<std::string> > >::iterator it = annotationMap.begin(); it != annotationMap.end(); it++) {
        std::vector<std::vector<std::string> > thisScaffoldAnnotation = it->second;
        std::vector<std::string> thisScaffoldTranscriptStartEnd;
        std::map<std::string, std::vector<std::string> > thisScaffoldTranscriptMap;
        for (std::vector<std::vector<std::string> >::size_type i = 0; i != thisScaffoldAnnotation.size(); i++) {
            std::vector<string> annotLineVec = split(thisScaffoldAnnotation[i][0], '\t');
            //if (i == 0) { print_vector_stream(annotLineVec, std::cerr); }
            string transcriptName = annotLineVec[4]; string strand = annotLineVec[3];
            string transcriptStart; string transcriptEnd;
            if (strand == "+") { transcriptStart = annotLineVec[1]; } else { transcriptEnd = annotLineVec[2]; }
            annotLineVec = split(thisScaffoldAnnotation[i].back(), '\t');
            if (strand == "+") { transcriptEnd = annotLineVec[2]; } else { transcriptStart = annotLineVec[1]; }
            // if (i == 0) { print_vector_stream(annotLineVec, std::cerr); }
            string transcriptStartEnd = transcriptName + "\t" + transcriptStart + "\t" + transcriptEnd + "\t" + strand;
            thisScaffoldTranscriptStartEnd.push_back(transcriptStartEnd);
            thisScaffoldTranscriptMap[transcriptName] = thisScaffoldAnnotation[i];
            //std::cerr << "Annotation processed: " << transcriptName << std::endl;
        }
        //std::cerr << "Annotation processed: " << it->first << std::endl;
        annotationMapTranscriptMap[it->first] = thisScaffoldTranscriptMap;
        transcriptStartEndMap[it->first] = thisScaffoldTranscriptStartEnd;
        thisScaffoldTranscriptMap.clear();
        thisScaffoldTranscriptStartEnd.clear();
        }
}
