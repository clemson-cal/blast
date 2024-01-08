CXXFLAGS = -std=c++17 -Ofast -Wall

blast: src/blast.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -I vapor/include -lhdf5
