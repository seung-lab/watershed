//#pragma once
#include "felzenszwalb.hpp"
#include "agglomeration.hpp"
#include "region_graph.hpp"
#include "basic_watershed.hpp"
#include "limit_functions.hpp"
#include "merge_tree.hpp"
#include "types.hpp"
#include "utils.hpp"
#include "boost/program_options.hpp"

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

namespace po = boost::program_options;

namespace
{
    const int SUCCESS = 0;
    const int ERROR_COMMAND_LINE = 1;
    const int ERROR_UNHANDLED_EXCEPTION = 2;
}

struct Config
{
    std::string inputFile;
    std::string outFileSegment;
    std::string outFileDendPairs;
    std::string outFileDendValues;
    size_t xSize;
    size_t ySize;
    size_t zSize;
    float lowv;
    float highv;
	bool enableMerge;
    size_t thold;
    int lowt;
    std::string funcName;
    double funcArg1;
    double funcArg2;
    double funcArg3;
};

int parseCmdLine(int argc, char *argv[], Config &config)
{
    po::options_description genericOptions("Generic options");
    genericOptions.add_options()
        ("help", "print help messages");

    po::options_description inputOptions("Input options");
    inputOptions.add_options()
        ("inputFile", po::value<std::string>(&(config.inputFile))->required(),
            "Filename of affinity graph input file.\n"
            "  Data: \tAffinity graph values\n"
            "  Data Format: \t 4-D float array of [3][xSize][ySize][zSize] *need to confirm this* read in fortran_storage_order"
        ) // TODO fortran_storage_order should be a configuration option
        ("xSize", po::value<size_t>(&config.xSize)->required(),"Dimension in X Direction.")
        ("ySize", po::value<size_t>(&config.ySize)->required(),"Dimension in Y Direction.")
        ("zSize", po::value<size_t>(&config.zSize)->required(),"Dimension in Z Direction.");

    po::options_description watershedOptions("Watershed Options");
    watershedOptions.add_options()
        ("lowv", po::value<float>(&config.lowv)->default_value(0.3),"Minimum threshold for watershed.")
        ("highv", po::value<float>(&config.highv)->default_value(0.9),"Maximum threshold for watershed.")
        ("enableMerge", po::value<bool>(&config.enableMerge)->default_value(true),"Enable merge region step for single linkage clustering")
        ("lowt", po::value<size_t>(&config.thold)->default_value(256),"Minimum merge size")
        ("thold", po::value<size_t>(&config.thold)->default_value(256),"Maximum merge size (calculated from --func")
        ("funcName", po::value<std::string>(&config.funcName)->default_value("constant"),"Merge thresholding function.\n"
            "** NOT IMPLEMENTED Defaults to const_above_threshold(.3, thold)**\n\t"
            "Example inputs:\n"
            "  --funcName=const --funcArg1=.3 --funcArg2=1000\n"
            "  --funcName=linear --funcArg1=100\n"
            "  --funcName=square --funcArg1=100\n"
            "  --funcName=power --funcArg1=2 --funcArg2=5000"
        )
        ("funcArg1", po::value<double>(&config.funcArg1)->default_value(.3),"Argument 1 for merge thresholding function")
        ("funcArg2", po::value<double>(&config.funcArg2),"Argument 2 for merge thresholding function")
        ("funcArg3", po::value<double>(&config.funcArg3),"Argument 3 for merge thresholding function");

    po::options_description outputOptions("Watershed Options");
    outputOptions.add_options()
        ("outFileSegment", po::value<std::string>(&config.outFileSegment)->
                default_value("ws.segment.data.out"),
            "Filename of the segmentation output file.\n"
            "  Data: \tSegmentation Ids for each voxel\n"
            "  Data Format: \t3-Dimensional uint32_t array [xSize][ySize][zSize] "
        )
        ("outFileDendPairs", po::value<std::string>(&config.outFileDendPairs)->
                default_value("ws.dend_pairs"),
            "Filename of the MST Dendrogram pairs file. Sorted by weight first as specified in "
            "the dendValues file.\n"
            "  Data: \tGraph Edges for MST (weights in another file)\n"
            "  Data Format: \t2-D array of uint32_t pairs i.e. [[SegId1, SegId2][SegId1, SegId3]...]"
        )
        ("outFileDendValues", po::value<std::string>(&config.outFileDendValues)->
                default_value("ws.dend_values"),
            "Filename of the MST Dendrogram values file.\n"
            "  Data: \tGraph Edge Weight Probabilities\n"
            "  Data Format: \t1-D array of decreasing floats corresponding to the pair of SegIds "
            "specified in the dendPairs file"
        );

    po::variables_map variableMap;
    po::positional_options_description positionalOptions;
    positionalOptions.add("inputFile", 1)
        .add("xSize", 1).add("ySize", 1).add("zSize", 1)
        .add("lowv", 1).add("highv", 1).add("lowt", 1).add("thold", 1)
        .add("funcName", 1).add("funcArg1", 1).add("funcArg2", 1).add("funcArg3", 1)
        .add("outFileSegment", 1).add("outFileDendPairs", 1).add("outFileDendValues", 1);

    po::options_description cmdlineOptions;
    cmdlineOptions.add(genericOptions).add(inputOptions).add(watershedOptions).add(outputOptions);

    try
    {
        po::store(po::command_line_parser(argc,argv).options(cmdlineOptions).positional(positionalOptions).run(), variableMap); // possible to throw
        if(variableMap.count("help"))
        {
            std::cout << "Basic command to run watershed on affifinity graph with single linkage clustering. "
                "Options:\n" << cmdlineOptions << "\nUsage Examples:\n\t"
                "runWatershedFull ws.affinity.data 256 256 256 0.3 0.9 250 10 ws.segment.data ws.dend_pairs ws.dend_values"
                << std::endl;
            return SUCCESS;
        }

        po::notify(variableMap); // throws on input error
    }
    catch(po::error& e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << cmdlineOptions << std::endl; 
      return ERROR_COMMAND_LINE; 
    }
    return SUCCESS;
}

/*
 * Run with args: inputFile outSegmentationFileName outFileDendPairs outFileDendValues xSize, ySize, zSize, lowv, highv, thold, lowt
 */
int main(int argc, char *argv[])
{
    Config config = Config();
    parseCmdLine(argc, argv, config);

    if (config.inputFile.empty())
    {
        return ERROR_COMMAND_LINE;
    }


    std::cout << "Reading in file " << config.inputFile 
        << " with x=" << config.xSize 
        << " y=" << config.ySize
        << " z=" << config.zSize << std::endl;
    affinity_graph_ptr<float> aff =
        read_affinity_graph<float>(config.inputFile, config.xSize, config.ySize, config.zSize);

    volume_ptr<uint32_t> seg;
    std::vector<std::size_t> counts;

    std::cout << "Running watershed with lowv=" << config.lowv << " highv=" << config.highv;
    std::tie(seg , counts) = watershed<uint32_t>(aff, config.lowv, config.highv);

    std::cout << "Getting region graph";
    auto rg = get_region_graph(aff, seg , counts.size()-1);

    if (config.enableMerge)
    {
        std::cout << "Merging with funcName=" << config.funcName 
            << "(" << config.funcArg1 << "," << config.funcArg2 << "," << config.funcArg3 << ")"
            << " thold=" << config.thold << " lowt= " << config.lowt;
        merge_segments_with_function (seg, rg, counts, const_above_threshold(config.funcArg1, config.thold), config.lowt);
    }
    else
    {
        std::cout << "Merging regions disabled" << std::endl;
    }

    auto mt = get_merge_tree(*rg, counts.size()-1);

    write_volume(config.outFileSegment, seg);

    write_region_graph_old_format(config.outFileDendPairs, config.outFileDendValues, *mt);
    return 0;
}
