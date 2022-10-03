//
//  UtilsGeneral.cpp
//  evogenSuite
//
//  Created by Milan Malinsky on 20.09.22.
//

#include "UtilsGeneral.hpp"
 
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

void splitToDouble(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(stringToDouble(item));
    }
}

std::vector<double> splitToDouble(const std::string &s, char delim) {
    std::vector<double> elems;
    splitToDouble(s, delim, elems);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// Initialize a matrix
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_rows, int m_columns) {
    if (m_columns == 0)
        m_columns = m_rows;
    for (int i = 0; i < m_columns; i++) {
        std::vector<double> v(m_rows,0);
        m.push_back(v);
    }
}

// Initialize a matrix
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_rows, int m_columns) {
    if (m_columns == 0)
        m_columns = m_rows;
    for (int i = 0; i < m_columns; i++) {
        std::vector<int> v(m_rows,0);
        m.push_back(v);
    }
}

// Remove a single file extension from the filename
std::string stripPath(const std::string& filename)
{
    size_t slashPos = filename.find_last_of('/');
    if(slashPos == std::string::npos)
        return filename; // no path
    else
        return filename.substr(slashPos+1);
}


// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}


bool testBiallelic(const std::string& altField) {
    std::vector<std::string> altVector = split(altField, ',');
    if (altVector.size() == 1) { return true; }
    else { return false; }
}

bool testOverallReadDepth(const int maxReadDepth,const int minReadDepth, const std::string& infoField) {
    std::vector<std::string> info = split(infoField, ';');
    if (info[0] == "INDEL") {
        split(info[1], '=', info);
    } else {
        split(info[0], '=', info);
    }
    int DP = atoi((info.back()).c_str());
    if (DP <= maxReadDepth && DP >= minReadDepth) { return true; }
    else { return false; }
}




// Does the same as R function table
std::map<int, int> tabulateVector(std::vector<int>& vec) {
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



// Move all doubleton counts to the bottom left corner of the matrix
void rearrange_doubletons(std::vector<std::vector<int> >& doubletons){
    for (std::vector<std::vector<int> >::size_type i = 0; i < doubletons.size(); i++) {
        for (int j = 0; j < i; j++) {
            doubletons[i][j] = doubletons[i][j] + doubletons[j][i];
            doubletons[j][i] = 0;
        }
    }
}


std::string suffix(const std::string& seq, size_t len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}


// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    size_t suffix_length = sizeof(GZIP_EXT) - 1;
    
    // Assume files without an extension are not compressed
    if(filename.length() < suffix_length)
        return false;
    
    std::string extension = suffix(filename, suffix_length);
    return extension == GZIP_EXT;
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }
}

//
void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        ogzstream* pGZ = new ogzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
        assertFileOpen(*pWriter, filename);
        return pWriter;
    }
}



void print80bpPerLineStdOut(std::ostream& outStream, string toPrint) {
    string::size_type lines = toPrint.length() / 80;
    for (string::size_type j = 0; j <= lines; j++) {
        outStream << toPrint.substr(j*80,80) << std::endl;
    }
}

void print80bpPerLineFile(std::ofstream*& outFile, string toPrint) {
    string::size_type lines = toPrint.length() / 80;
    for (string::size_type j = 0; j <= lines; j++) {
        *outFile << toPrint.substr(j*80,80) << std::endl;
    }
}



std::map<string,string> readMultiFastaToMap(const string& fileName) {
    std::map<string, string> fastaSeqs;
    string line;
    std::ifstream* fastaFile = new std::ifstream(fileName.c_str());
    getline(*fastaFile, line);
    std::vector<std::string> fields = split(line, ' ');
    string currentScaffold = fields[0].substr(1,string::npos);
    fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
    while (getline(*fastaFile, line)) {
        if (line[0] != '>') {
            fastaSeqs[currentScaffold].append(line);
        } else {
            int l = (int)fastaSeqs[currentScaffold].length();
            std::cerr << currentScaffold << " length: " << l << std::endl;
            fastaSeqs[currentScaffold].shrink_to_fit();
            std::vector<std::string> fields = split(line, ' ');
            currentScaffold = fields[0].substr(1,string::npos);
            fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
        }
    }
    std::cerr << currentScaffold << " length: " << (int)fastaSeqs[currentScaffold].length() << std::endl;
    return fastaSeqs;
}

std::string readMultiFastaToOneString(const string& fileName, int bytes) {
    std::string fastaSeqs; fastaSeqs.reserve(bytes); fastaSeqs = "";
    string line;
    std::ifstream* fastaFile = new std::ifstream(fileName.c_str());
    getline(*fastaFile, line);
    string currentScaffold = line.substr(1,string::npos);
    while (getline(*fastaFile, line)) {
        if (line[0] != '>') {
            fastaSeqs.append(line);
        } else {
            continue;
        }
    }
    return fastaSeqs;
}


double getExpectedGenotype(const std::vector<double>& thisProbabilities) {
    double Egenotype = thisProbabilities[1] + 2*thisProbabilities[2];
    return Egenotype;
}

void transformFromPhred(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,-(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,-(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,-(thisLikelihoods[2]/10.0));
}

void transformFromGL(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,(thisLikelihoods[2]/10.0));
}

