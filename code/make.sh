mkdir -p bin
g++ src/runWatershed.cpp -I. -I./src -O3 -DNDEBUG -std=c++11 -o ./bin/rws -lboost_program_options
