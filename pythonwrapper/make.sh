g++ watershed/pywatershed.cpp -std=c++11 -Wall -I/usr/include/python2.7 -I../code/src -I../code -x c++ -lboost_python -lboost_numpy -lpython2.7 -v -fPIC -o PyWatershed.so -g -shared
