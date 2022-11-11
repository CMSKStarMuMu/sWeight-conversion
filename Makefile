#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore -lMinuit
ROOTCINT  := $(shell which rootcint)

#exe_files
EXECUTABLE1 := convert_sWeights
EXECUTABLE2 := merge_convert_sWeights
EXECUTABLE3 := merge2_convert_sWeights

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 


$(EXECUTABLE1): $(EXECUTABLE1).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE2): $(EXECUTABLE2).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE3): $(EXECUTABLE3).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS)

#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3)
