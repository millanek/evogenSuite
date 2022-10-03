//
//  UtilsGeneral.hpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#ifndef UtilsGeneral_hpp
#define UtilsGeneral_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <getopt.h>
#include <cstdlib>
#include <stdexcept>
#include <math.h>
#include <random>
#include <unordered_map>
#include <ctime>
#include <deque>
#include "gzstream.hpp"

using std::string;
#define PROGRAM_BIN "evogenSuite"
#define PACKAGE_BUGREPORT "millanek@gmail.com"
#define GZIP_EXT ".gz"

#define LikelihoodsProbabilitiesAbsent 0
#define LikelihoodsProbabilitiesGP 1
#define LikelihoodsProbabilitiesGL 2
#define LikelihoodsProbabilitiesPL 3

#define AncestralAlleleMissing -1
#define AncestralAlleleRef 0
#define AncestralAlleleAlt 1


// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

// -----------------------------------------------------------------------
// --------------------- Workflow Info Printing --------------------------
// -----------------------------------------------------------------------

inline void reportProgessVCF(const int variantsProcessed, const std::clock_t startTime) {
    double durationOverall = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    std::cout << "Processed " << variantsProcessed << " variants in " << durationOverall << "secs" << std::endl;
}

inline void printMissingLikelihoodsWarning(const string& chr, const string& pos) {
    std::cerr << "WARNING: Could not find genotype likelihoods/probabilities (GP, PL, or GL fields) for variant at " << chr << " " << pos << std::endl;
    std::cerr << "WARNING: Did you really mean to use the -g option? Reverting to using called genotypes." << std::endl;
}

// -----------------------------------------------------------------------
// --------------------- Data Manipulation/Conversion --------------------
// -----------------------------------------------------------------------

void splitToDouble(const std::string &s, char delim, std::vector<double> &elems);
std::vector<double> splitToDouble(const std::string &s, char delim);
std::vector<std::string> split(const std::string &s, char delim);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}

class BadConversion : public std::runtime_error {
public:
    BadConversion(std::string const& s)
    : std::runtime_error(s)
    { }
};

inline double stringToDouble(std::string const& s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

// -----------------------------------------------------------------------
// --------------------- Vector Calculations -----------------------------
// -----------------------------------------------------------------------

template <class T> double vector_sum(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    return sum;
}

template <class T> double vector_average(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)vector.size();
    return average;
}

template <class T> double vector_average_withRegion(T vector, int regionLength) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)regionLength;
    return average;
}


template <typename T> std::map<int, int> tabulateVectorTemplate(T& vec) {
    std::vector<int> vecCopy(vec);
    std::sort(vecCopy.begin(), vecCopy.end());
    std::vector<int>::iterator it = std::unique(vecCopy.begin(), vecCopy.end());
    vecCopy.resize(std::distance(vecCopy.begin(), it));
    
    std::map<int, int>  table;
    //int pos = 0;
    for (std::vector<int>::size_type i = 0; i != vecCopy.size(); i++) {
        int mycount = std::count(vec.begin(), vec.end(), vecCopy[i]);
        table[vecCopy[i]] = mycount;
        //pos = pos + mycount;
    }
    return table;
}

// Does the same as R function table
std::map<int, int> tabulateVector(std::vector<int>& vec);

// -----------------------------------------------------------------------
// --------------------- File Manipulation -------------------------------
// -----------------------------------------------------------------------

std::string stripExtension(const std::string& filename);
std::string stripPath(const std::string& filename);
void assertGZOpen(gzstreambase& gh, const std::string& fn);
void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// -----------------------------------------------------------------------------
// --------------------- Printing and manipulating Vectors/Matrices ------------
// -----------------------------------------------------------------------------

// Print an arbitrary matrix (vector of vectors)
template <class T> void print_matrix(T matrix, std::ostream& outFile) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (j == (matrix[i].size()-1))
                outFile << matrix[i][j] << std::endl;
            else
                outFile << matrix[i][j] << "\t";
        }
    }
}

// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ostream& outFile, char delim = '\t', bool endLine = true) {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1)) {
            if (endLine) outFile << vector[i] << std::endl;
            else outFile << vector[i];
        } else {
            outFile << vector[i] << delim;
        }
    }
}

template <typename T> void reset_matrix_to_zero(std::vector<std::vector<T> >& m) {
    for (int i = 0; i < m.size(); i++) {
        std::fill(m[i].begin(), m[i].end(), 0);
    }
}

void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_rows, int m_columns = 0);
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_rows, int m_columns = 0);

void transformFromPhred(std::vector<double>& thisLikelihoods);
void transformFromGL(std::vector<double>& thisLikelihoods);
double getExpectedGenotype(const std::vector<double>& thisProbabilities);

bool testBiallelic(const std::string& altField);

void print80bpPerLineStdOut(std::ostream& outFile, std::string toPrint);
void print80bpPerLineFile(std::ofstream*& outFile, string toPrint);


// Reading in genome files
std::map<std::string,std::string> readMultiFastaToMap(const std::string& fileName);
std::string readMultiFastaToOneString(const string& fileName, int bytes = 50000000);

#endif /* UtilsGeneral_hpp */
