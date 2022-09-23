
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/evogenSuite

$(BIN)/evogenSuite: $(BIN) $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsStats.o $(BIN)/UtilsGeneral.o $(BIN)/gzstream.o
	$(CXX) $(CXXFLAGS) $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsStats.o $(BIN)/UtilsGeneral.o $(BIN)/gzstream.o -o $@ $(LDFLAGS)

$(BIN)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN):
	mkdir -p $@

clean:
	rm $(BIN)/*.o $(BIN)/evogenSuite
