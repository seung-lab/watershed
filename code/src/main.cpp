//#pragma once
#include "felzenszwalb.hpp"
#include "agglomeration.hpp"
#include "region_graph.hpp"
#include "basic_watershed.hpp"
#include "limit_functions.hpp"
#include "merge_tree.hpp"
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

std::size_t get_rand_idx( std::vector<uint32_t>& v)
{
    if ( v.size() <= 1 )
    {
        return 0;
    }

    std::size_t r = 0;

    for ( int i = 0; i < v.size() - 1; ++i )
    {
        for ( int j = i+1; j < v.size(); ++j )
            r += v[i]*v[j];
    }

    return r;
}

double square_sum( std::vector<uint32_t>& v)
{
    double r = 0;

    for ( const auto& x: v )
    {
        r += x*x;
    }

    return r;
}

template< typename ID >
void fill_void( ID* arr, std::size_t len )
{
    ID maxi = 0;
    for ( std::size_t i = 0; i < len; ++i )
    {
        maxi = std::max(maxi, arr[i]);
    }

    for ( std::size_t i = 0; i < len; ++i )
    {
        if ( arr[i] == 0 ) arr[i] = ++maxi;
    }
}

std::pair<double,double>
compare_volumes( volume_ref<uint32_t>& gt,
                 volume_ref<uint32_t>& ws,
                 std::size_t size )
{
    std::map<uint32_t, std::map<uint32_t, uint32_t>> map;
    std::map<uint32_t, std::map<uint32_t, uint32_t>> invmap;
    std::map<uint32_t, uint32_t> setg, setw;

    double rand_split = 0;
    double rand_merge = 0;

    double t_sq = 0;
    double s_sq = 0;

    double total = 0;
    std::map<uint32_t, std::map<uint32_t, std::size_t>> p_ij;

    std::map<uint32_t, std::size_t> s_i, t_j;

    for ( std::ptrdiff_t z = 28; z < size-28; ++z )
        for ( std::ptrdiff_t y = 28; y < size-28; ++y )
            for ( std::ptrdiff_t x = 16; x < size-16; ++x )
            {
                uint32_t wsv = ws[x][y][z];
                uint32_t gtv = gt[x][y][z];

                if ( gtv )
                {
                    ++total;

                    ++p_ij[gtv][wsv];
                    ++s_i[wsv];
                    ++t_j[gtv];
                }
            }

    double sum_p_ij = 0;
    for ( auto& a: p_ij )
    {
        for ( auto& b: a.second )
        {
            sum_p_ij += b.second * b.second;
        }
    }

    double sum_t_k = 0;
    for ( auto& a: t_j )
    {
        sum_t_k += a.second * a.second;
    }


    double sum_s_k = 0;
    for ( auto& a: s_i )
    {
        sum_s_k += a.second * a.second;
    }


    //std::cout << sum_p_ij << "\n";
    std::cout << "Rand Split: " << (sum_p_ij/sum_t_k) << "\n";
    std::cout << "Rand Merge: " << (sum_p_ij/sum_s_k) << "\n";
    std::cout << "Rand alpha: " << (sum_p_ij*2/(sum_t_k+sum_s_k)) << "\n";


    return std::make_pair(sum_p_ij/sum_t_k,
                          sum_p_ij/sum_s_k);
}

std::vector<double> reduce( const std::vector<double>& v )
{
    std::vector<double> ret;
    ret.push_back(v[0]);
    ret.push_back(v[1]);

    for ( std::size_t i = 2; i + 4 < v.size(); i += 2 )
    {
        double oldx = v[i-2];
        double oldy = v[i-1];

        double x = v[i];
        double y = v[i+1];

        double nextx = v[i+2];
        double nexty = v[i+3];

        if ( std::abs(oldx-x) > 0.0001 ||
             std::abs(oldy-y) > 0.0001 ||
             std::abs(nextx-x) > 0.0001 ||
             std::abs(nexty-y) > 0.0001 )
        {
            ret.push_back(x);
            ret.push_back(y);
        }
    }

    ret.push_back(v[v.size()-2]);
    ret.push_back(v[v.size()-1]);

    return ret;
}



////////////////////////////////////////////////////////////////////////////
//
// Felzenszwalb implementation

void process_felzenszwalb( volume_ptr<uint32_t> gt_ptr,
                           affinity_graph_ptr<float> aff,
                           double k,
                           const std::string& fname )
{
    auto seg = felzenszwalb<uint32_t>(aff, k);
    write_volume(fname + ".dat", seg);

    auto prc = reduce(felzenszwalb_err(gt_ptr, aff, k));
    write_to_file(fname + "_pr.dat", prc.data(), prc.size());
}

////////////////////////////////////////////////////////////////////////////
//
// Const above threshold

void process_const_above_thold( volume_ptr<uint32_t> gt_ptr,
                                affinity_graph_ptr<float> aff,
                                float thold,
                                std::size_t sz,
                                const std::string& fname )
{
    volume_ptr<uint32_t>     segg  ;
    std::vector<std::size_t> counts;

    {
        std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);
        auto rg = get_region_graph(aff, segg, counts.size()-1);

        merge_segments_with_function
            (segg, rg, counts,
             const_above_threshold(thold, sz), 100);

        write_volume(fname + ".dat", segg);
    }

    {
        std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);
        auto rg = get_region_graph(aff, segg, counts.size()-1);

        auto r = reduce(merge_segments_with_function_err
                        (segg, gt_ptr, rg, counts,
                         const_above_threshold(thold, sz), 100));

        write_to_file(fname + "_pr.dat", r.data(), r.size());
    }
}


////////////////////////////////////////////////////////////////////////////
//
// Const above threshold

void process_square( volume_ptr<uint32_t> gt_ptr,
                     affinity_graph_ptr<float> aff,
                     std::size_t atlow,
                     std::size_t at1,
                     const std::string& fname )
{
    volume_ptr<uint32_t>     segg  ;
    std::vector<std::size_t> counts;

    {
        std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);
        auto rg = get_region_graph(aff, segg, counts.size()-1);

        merge_segments_with_function
            (segg, rg, counts,
             square_fn(0.3, atlow, at1), 100);

        write_volume(fname + ".dat", segg);
    }

    {
        std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);
        auto rg = get_region_graph(aff, segg, counts.size()-1);

        auto r = reduce(merge_segments_with_function_err
                        (segg, gt_ptr, rg, counts,
                         square_fn(0.3, atlow, at1), 100));

        write_to_file(fname + "_pr.dat", r.data(), r.size());
    }
}





////////////////////////////////////////////////////////////////////////////
//
// Felzenszwalb implementation

// void process_square_fn( volume_ptr<uint32_t> gt_ptr,
//                         affinity_graph_ptr<float> aff,
//                         std::size_t k,
//                         const std::string& fname )
// {
//     volume_ptr<uint32_t>     segg  ;
//     std::vector<std::size_t> counts;

//     {
//         std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);
//         auto rg = get_region_graph(aff, segg, counts.size()-1);

//         auto r = reduce(merge_segments_with_function_err
//                         (segg, gt_ptr, rg, counts,
//                          square_fn(0.3, 50000, 250), 100));

//         write_to_file("./experiments/square/50000_pr.dat",
//                       r.data(), r.size());

//         write_volume("./experiments/square/50000.dat", segg);
//     }

//     auto seg = felzenszwalb<uint32_t>(aff, k);
//     write_volume(fname + ".dat", seg);

//     auto prc = reduce(felzenszwalb_err(gt_ptr, aff, k));
//     write_to_file(fname + "_pr.dat", prc.data(), prc.size());
// }

int main()
{
    // load the ground truth and the affinity graph

    volume_ptr<uint32_t> gt_ptr =
        read_volume<uint32_t>("./data/gt.in", 256);


    affinity_graph_ptr<float> aff =
        read_affinity_graph<float>("./data/ws_test_256.raw",
                                   256, 256, 256);


    if ( 1 )
    {
        volume_ptr<uint32_t>     seg   ;
        std::vector<std::size_t> counts;

        std::tie(seg, counts) = watershed<uint32_t>(aff, -1, 2);
        write_volume("./experiments/watershed/basic.out", seg);

        std::tie(seg, counts) = watershed<uint32_t>(aff, 0.1, 0.99);
        write_volume("./experiments/watershed/minmax.out", seg);

        // {
        //     std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);

        //     auto rg = get_region_graph(aff, segg, counts.size()-1);
        //     merge_segments_with_function(segg, rg,
        //                                  counts, limit_fn2, 100);

        //     write_volume("./experiments/voutall.out", segg);
        // }

        //return 0;
    }


    //
    // Linear
    //
    if ( 1 )
    {
        std::vector<double> r;
        for ( std::size_t thold = 200; thold <= 100000; thold += 100 )
        {
            if ( thold > 1000 ) thold += 900;
            if ( thold > 10000 ) thold += 9000;

            //double k = static_cast<double>(thold) / 1000;

            std::cout << "THOLD: " << thold << "\n";

            volume_ptr<uint32_t>     seg   ;
            std::vector<std::size_t> counts;

            {
                std::tie(seg , counts) = watershed<uint32_t>(aff, 0.3, 0.99);
                auto rg = get_region_graph(aff, seg , counts.size()-1);

                merge_segments_with_function
                    (seg, rg, counts,
                     linear(thold), 10);

                write_volume("experiments/linear/"
                             + std::to_string(thold) + ".dat", seg);

                auto x = compare_volumes(*gt_ptr, *seg, 256);
                r.push_back(x.first);
                r.push_back(x.second);
            }
            write_to_file("experiments/linear.dat", r.data(), r.size());
        }

        //return 0;
    }

    //
    // Square
    //
    if ( 1 )
    {
        std::vector<double> r;
        for ( std::size_t thold = 200; thold <= 100000; thold += 100 )
        {
            if ( thold > 1000 ) thold += 900;
            if ( thold > 10000 ) thold += 9000;

            std::cout << "THOLD: " << thold << "\n";

            volume_ptr<uint32_t>     seg   ;
            std::vector<std::size_t> counts;

            {
                std::tie(seg , counts) = watershed<uint32_t>(aff, 0.3, 0.99);
                auto rg = get_region_graph(aff, seg , counts.size()-1);

                merge_segments_with_function
                    (seg, rg, counts,
                     square(thold), 10);

                write_volume("experiments/square/"
                             + std::to_string(thold) + ".dat", seg);

                auto x = compare_volumes(*gt_ptr, *seg, 256);
                r.push_back(x.first);
                r.push_back(x.second);
            }
            write_to_file("experiments/square.dat", r.data(), r.size());
        }

//        return 0;
    }

    //
    // Felzenszwalb implementation
    //
    if ( 1 )
    {
        std::vector<double> r;
        for ( std::size_t thold = 100; thold <= 50000; thold += 100 )
        {
            if ( thold > 1000 ) thold += 900;
            if ( thold > 10000 ) thold += 9000;

            double k = static_cast<double>(thold) / 1000;

            auto seg = felzenszwalb<uint32_t>(aff, k);
            write_volume("experiments/felzenszwalb/"
                         + std::to_string(k) + ".dat", seg);

            auto x = compare_volumes(*gt_ptr, *seg, 256);

            r.push_back(x.first);
            r.push_back(x.second);
        }
        write_to_file("experiments/felzenszwalb.dat", r.data(), r.size());
    }

    //
    // simple thold fn
    //
    if ( 1 )
    {


        std::vector<double> r;
        for ( std::size_t thold = 100; thold <= 50000; thold += 100 )
        {
            if ( thold > 1000 ) thold += 900;
            if ( thold > 10000 ) thold += 9000;

            std::cout << "THOLD: " << thold << "\n";

            volume_ptr<uint32_t>     seg   ;
            std::vector<std::size_t> counts;

            {
                std::tie(seg , counts) = watershed<uint32_t>(aff, 0.3, 0.99);
                auto rg = get_region_graph(aff, seg , counts.size()-1);

                merge_segments_with_function
                    (seg, rg, counts,
                     const_above_threshold(0.3, thold), 100);

                write_volume("experiments/threshold/"
                             + std::to_string(thold) + ".dat", seg);

                auto x = compare_volumes(*gt_ptr, *seg, 256);
                r.push_back(x.first);
                r.push_back(x.second);
            }
        }
        write_to_file("experiments/threshold.dat", r.data(), r.size());
    }

    return 0;

    volume_ptr<uint32_t>     segg  ;
    std::vector<std::size_t> counts;

    std::tie(segg, counts) = watershed<uint32_t>(aff, -1, 2);

    write_volume("voutraw.out", segg);

    for ( float low = 0.01; low < 0.051; low += 0.01 )
    {
        for ( float high = 0.998; high > 0.989; high -= 0.002 )
        {
//            std::tie(segg, counts) = watershed<uint32_t>(aff, low, high);
//            write_volume("vout." + std::to_string(low) + "." +
//                                 std::to_string(high) + ".out", segg);
        }
    }

    std::tie(segg, counts) = watershed<uint32_t>(aff, 0.5, 2);

    write_volume("voutmax.out", segg);

    std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);

    write_volume("voutminmax.out", segg);


//    return 0;

    // auto rg = get_region_graph(aff, segg, counts.size()-1);

    // //yet_another_watershed(segg, rg, counts, 0.3);

    // //write_volume("voutanouther.out", segg);


    // merge_segments_with_function(segg, rg, counts, limit_fn3, 100);

    // write_volume("voutdo.out", segg);




    auto rg = get_region_graph(aff, segg, counts.size()-1);

    //yet_another_watershed(segg, rg, counts, 0.3);

    //write_volume("voutanouther.out", segg);

    auto r = merge_segments_with_function_err(segg, gt_ptr, rg,
                                              counts, limit_fn2, 100);

    write_to_file("./experiments/custom/precision_recall.dat",
                  r.data(), r.size());

    write_volume("voutall4x.out", segg);

    return 0;

    write_region_graph("voutall.rg", *rg);

    auto mt = get_merge_tree(*rg, counts.size()-1);

    write_region_graph("voutall.mt", *mt);

    return 0;


}
