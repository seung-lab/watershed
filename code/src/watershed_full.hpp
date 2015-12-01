#pragma once

#include "agglomeration.hpp"
#include "region_graph.hpp"
#include "basic_watershed.hpp"
#include "limit_functions.hpp"
#include "merge_tree.hpp"
#include "types.hpp"
#include "watershed_config.hpp"
#include <cstddef>
#include <iostream>
#include <map>
std::tuple<volume_ptr<uint32_t>, region_graph_ptr<uint32_t, float>> watershed_full(
        const affinity_graph_ptr<float>& aff, const ws::Config& config)
{
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

    region_graph_ptr<uint32_t, float> mt = get_merge_tree(*rg, counts.size()-1);

    return std::make_tuple(seg, mt);
}

