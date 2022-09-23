
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/evogenSuite

$(BIN)/evogenSuite: $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsStats.o $(BIN)/UtilsGeneral.o $(BIN)/gzstream.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(BIN)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN)/evogenSuite: $(BIN)/evogenSuite.o $(BIN)/AlleleFreq.o $(BIN)/Fst.o $(BIN)/PBS.o $(BIN)/UtilsAnnotation.o $(BIN)/UtilsIUPAC.o $(BIN)/UtilsGeneral.o $(BIN)/UtilsStats.o $(BIN)/gzstream.o | $(BIN)

$(BIN):
	mkdir -p $@

clean:
	rm $(BIN)/*.o $(BIN)/evogenSuite
