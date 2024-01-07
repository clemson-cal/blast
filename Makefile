CXXFLAGS = -std=c++17 -Ofast -Wall

shock: src/shock.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -I vapor/include -lhdf5
