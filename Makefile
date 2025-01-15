STIM_DIR = ./Stim/src/
STIM_LIB_PATH = ./Stim/build/out/
STIM_INCLUDE = $(STIM_DIR)

CXX = g++
CXXFLAGS = -I$(STIM_INCLUDE) -std=c++23
LDFLAGS = -L$(STIM_LIB_PATH) -lstim

SRC = stim_test.cpp
EXEC = stim_test

$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXEC) $(LDFLAGS)

clean:
	rm -f $(EXEC)
