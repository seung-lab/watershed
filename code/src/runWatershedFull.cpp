//#pragma once
#include "watershed_full.hpp"
#include "watershed_config.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <memory>
#include <type_traits>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <chrono>
#include <stdlib.h>
#include <cstdlib>

namespace
{
    const int ERROR_COMMAND_LINE = 1;
    const int ERROR_UNHANDLED_EXCEPTION = 2;
}

int read_write_watershed(ws::Config &config)
{
    std::cout << "Reading in file " << config.inputFile 
        << " with x=" << config.xSize 
        << " y=" << config.ySize
        << " z=" << config.zSize << std::endl;
    affinity_graph_ptr<float> aff =
        read_affinity_graph<float>(config.inputFile, config.xSize, config.ySize, config.zSize);

    volume_ptr<uint32_t> seg;
    region_graph_ptr<uint32_t, float> mst;
    std::tie(seg, mst) = watershed_full<uint32_t, float>(aff, config);

    write_volume(config.outFileSegment, seg);

    write_region_graph_old_format(config.outFileDendPairs, config.outFileDendValues, *mst);
    return EXIT_SUCCESS;
}

/*
 * Run with args: inputFile outSegmentationFileName outFileDendPairs outFileDendValues xSize, ySize, zSize, lowv, highv, thold, lowt
 */
int main(int argc, char *argv[])
{
    ws::Config config = ws::Config();
    int status = ws::config::parseCmdLine(argc, argv, config);

    if (status != EXIT_SUCCESS) {
        return status;
    }
    if (config.inputFile.empty())
    {
        return ERROR_COMMAND_LINE;
    }

    return read_write_watershed(config);
}
