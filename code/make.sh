mkdir -p bin
g++ src/runWatershedFull.cpp -I. -I./src -O3 -DNDEBUG -std=c++11 -o ./bin/runWatershedFullDouble -lboost_program_options
