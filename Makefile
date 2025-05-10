
CXXFLAGS=-std=c++17
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/evogenSuite

$(BIN)/evogenSuite: $(BIN) $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsStats.o $(BIN)/UtilsGeneral.o $(BIN)/gzstream.o $(BIN)/UtilsSetCounts.o $(BIN)/UtilsSetInfo.o $(BIN)/DistanceMatrix.o $(BIN)/DistanceToOutgroups.o $(BIN)/PaintWithFixedSites.o $(BIN)/PhysicalWindowAverages.o $(BIN)/CodingStats.o $(BIN)/UtilsCodingStats.o
	$(CXX) $(CXXFLAGS) $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsStats.o $(BIN)/UtilsGeneral.o $(BIN)/gzstream.o $(BIN)/UtilsSetCounts.o $(BIN)/UtilsSetInfo.o $(BIN)/DistanceMatrix.o $(BIN)/DistanceToOutgroups.o $(BIN)/PaintWithFixedSites.o $(BIN)/PhysicalWindowAverages.o $(BIN)/CodingStats.o $(BIN)/UtilsCodingStats.o -o $@ $(LDFLAGS)

$(BIN)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN):
	mkdir -p $@

clean:
	rm $(BIN)/*.o $(BIN)/evogenSuite
