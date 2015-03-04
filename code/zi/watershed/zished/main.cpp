//#pragma once

#include "agglomeration.hpp"
#include "region_graph.hpp"
#include "basic_watershed.hpp"
#include "limit_functions.hpp"
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


template< typename ID, typename F >
inline region_graph_ptr<ID,F>
get_merge_tree( const region_graph<ID,F>& rg, std::size_t max_segid )
{
    zi::disjoint_sets<ID>      sets(max_segid+1);
    std::vector<std::list<ID>> edges(max_segid+1);
    region_graph_ptr<ID,F>     mt_ptr( new region_graph<ID,F> );

    for ( const auto& e: rg )
    {
        ID v1 = std::get<1>(e);
        ID v2 = std::get<2>(e);

        ID s1 = sets.find_set(v1);
        ID s2 = sets.find_set(v2);

        if ( s1 != s2 && s1 && s2 )
        {
            mt_ptr->push_back(e);
            sets.join(s1, s2);

            edges[v1].push_back(v2);
            edges[v2].push_back(v1);

        }
    }

    std::vector<ID> order(max_segid+1);
    ID curr = 0;

    for ( ID i = 0; i <= max_segid; ++i )
    {
        if ( order[i] == 0 )
        {
            std::deque<ID> queue;
            queue.push_back(i);
            order[i] = ++curr;

            while ( queue.size() )
            {
                ID x = queue.front();
                queue.pop_front();

                for ( auto& y: edges[x] )
                {
                    if ( order[y] == 0 )
                    {
                        order[y] = ++curr;
                        queue.push_back(y);
                    }
                }
            }
        }
    }

    for ( auto& e: *mt_ptr )
    {
        if ( order[std::get<2>(e)] < order[std::get<1>(e)] )
        {
            std::swap(std::get<2>(e), std::get<1>(e));
        }
    }

    return mt_ptr;
}



template< typename ID, typename F, typename K >
inline volume_ptr<ID> other_alg( const affinity_graph_ptr<F>& aff_ptr,
                                 const K& kconst )
{
    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    auto dummy          = get_dummy_segmentation<ID>(xdim, ydim, zdim);
    volume<ID>&    seg  = *std::get<0>(dummy);

    region_graph_ptr<ID,F> rg_ptr = get_region_graph(aff_ptr,
                                                     std::get<0>(dummy), total);
    region_graph<ID,F>& rg  = *rg_ptr;

    zi::disjoint_sets<ID> sets(total+1);

    F k = kconst;

    std::vector<F> int_val(total+1);
    std::vector<F> sizes(total+1);

    std::fill_n(sizes.begin(), total+1, 1);

    for ( auto& it: rg )
    {
        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        F  w  = static_cast<F>(1) - std::get<0>(it);

        if ( s1 != s2 )
        {
            F v1 = int_val[s1] + k / sizes[s1];
            F v2 = int_val[s2] + k / sizes[s2];

            if ( w < std::min(v1, v2) )
            {
                sizes[s1] += sizes[s2];
                sizes[s2] = 0;

                int_val[s1] = w;
                int_val[s2] = 0;

                ID s = sets.join(s1,s2);

                std::swap(sizes[s1], sizes[s]);
                std::swap(int_val[s1], int_val[s]);
            }
        }
    }

    std::cout << "Done with merging" << std::endl;

    std::vector<ID> remaps(sizes.size());

    ID next_id = 1;

    for ( ID id = 0; id < sizes.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) )
        {
            remaps[s] = next_id;
            ++next_id;
        }
    }

    ID* seg_raw = seg.data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps[sets.find_set(seg_raw[idx])];
    }

    std::cout << "Done with remapping, total: " << (next_id-1) << std::endl;

    region_graph<ID,F> new_rg;

    for ( auto& it: rg )
    {
        ID s1 = remaps[sets.find_set(std::get<1>(it))];
        ID s2 = remaps[sets.find_set(std::get<2>(it))];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            new_rg.emplace_back(std::get<0>(it), mm.first, mm.second);
        }
    }

    rg.swap(new_rg);

    std::cout << "Done with updating the region graph, size: "
              << rg.size() << std::endl;

    return std::get<0>(dummy);
}

template< typename ID,
          typename F,
          typename L,
          class = typename std::enable_if<is_numeric<F>::value>::type,
          class = typename std::enable_if<std::is_integral<ID>::value>::type,
          class = typename std::enable_if<std::is_convertible<L,F>::value>::type >
inline void yet_another_watershed( const volume_ptr<ID>& seg_ptr,
                                   const region_graph_ptr<ID,F> rg_ptr,
                                   std::vector<std::size_t>& counts,
                                   const L& lowl)
{
    F low = static_cast<F>(lowl);

    std::vector<std::size_t> new_counts({0});
    std::vector<ID>          remaps(counts.size());

    zi::disjoint_sets<ID>    sets(counts.size());
    std::vector<F>           maxs(counts.size());


    region_graph<ID,F>& rg  = *rg_ptr;

    ID next_id = 1;
    ID merged = 0;


    for ( auto& it: rg )
    {
        if ( std::get<0>(it) <= low )
        {
            break;
        }

        ID s1 = std::get<1>(it);
        ID s2 = std::get<2>(it);

        F f = std::get<0>(it);

        if ( s1 && s2 )
        {
            if ( (remaps[s1] == 0) || (remaps[s2] == 0) )
            {
                if ( remaps[s1] == 0 )
                {
                    std::swap(s1,s2);
                }

                if ( remaps[s1] == 0 )
                {
                    maxs[next_id] = f;
                    remaps[s1] = remaps[s2] = next_id;
                    new_counts.push_back(counts[s1]+counts[s2]);
                    ++next_id;
                }
                else
                {
                    ID actual = sets.find_set(remaps[s1]);
                    remaps[s2] = remaps[s1];
                    new_counts[actual] += counts[s2];
                }
            }
            else
            {
                ID a1 = sets.find_set(remaps[s1]);
                ID a2 = sets.find_set(remaps[s2]);

                if ( 0 && a1 != a2 && ((maxs[a1]==f)||(maxs[a2]==f)) )
                {
                    ++merged;
                    new_counts[a1] += new_counts[a2];
                    new_counts[a2] = 0;
                    maxs[a1] = std::max(maxs[a1],maxs[a2]);
                    maxs[a2] = 0;
                    ID a = sets.join(a1,a2);
                    std::swap(new_counts[a], new_counts[a1]);
                    std::swap(maxs[a], maxs[a1]);
                }
            }
        }
    }

    next_id -= merged;

    std::vector<ID> remaps2(counts.size());

    next_id = 1;

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(remaps[id]);
        if ( s && (remaps2[s]==0) )
        {
            remaps2[s] = next_id;
            new_counts[next_id] = new_counts[s];
            ++next_id;
        }
    }

    new_counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps2[remaps[seg_raw[idx]]];
    }

    region_graph<ID,F> new_rg;

    for ( auto& it: rg )
    {
        ID s1 = remaps2[remaps[std::get<1>(it)]];
        ID s2 = remaps2[remaps[std::get<2>(it)]];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            new_rg.emplace_back(std::get<0>(it), mm.first, mm.second);
        }
    }

    rg.swap(new_rg);

    counts.swap(new_counts);

    std::cout << "New count: " << counts.size() << std::endl;

    std::cout << "Done with updating the region graph, size: "
              << rg.size() << std::endl;
}


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

void compare_volumes( const std::string& gtf, const std::string& wsf,
                      std::size_t size )
{
    volume_ptr<uint32_t> gt_ptr = read_volume<uint32_t>(gtf, size);
    volume_ptr<uint32_t> ws_ptr = read_volume<uint32_t>(wsf, size);
    //volume_ptr<uint32_t> ws_ptr = get_dummy_segmentation<uint32_t>(size,size,size);

    //fill_void(gt_ptr->data(), size*size*size);
    //fill_void(ws_ptr->data(), size*size*size);


    volume<uint32_t>& gt = *gt_ptr;
    volume<uint32_t>& ws = *ws_ptr;

    std::map<uint32_t, std::map<uint32_t, uint32_t>> map;
    std::map<uint32_t, std::map<uint32_t, uint32_t>> invmap;
    std::map<uint32_t, uint32_t> setg, setw;

    double rand_split = 0;
    double rand_merge = 0;

    double t_sq = 0;
    double s_sq = 0;

    double total = 0;
    std::map<uint32_t, std::map<uint32_t, double>> p_ij;

    std::map<uint32_t, double> s_i, t_j;

    for ( std::ptrdiff_t z = 28; z < size-28; ++z )
        for ( std::ptrdiff_t y = 28; y < size-28; ++y )
            for ( std::ptrdiff_t x = 16; x < size-16; ++x )
            {
                uint32_t wsv = ws[x][y][z];
                uint32_t gtv = gt[x][y][z];

                if ( gtv )
                {
                    total += 1;

                    p_ij[gtv][wsv] += 1;

                    s_i[wsv] += 1;
                    t_j[gtv] += 1;
                }
            }

    double sum_p_ij = 0;
    for ( auto& a: p_ij )
    {
        for ( auto& b: a.second )
        {
            sum_p_ij += (b.second / total) * (b.second / total);
        }
    }

    double sum_t_k = 0;
    for ( auto& a: t_j )
    {
        sum_t_k += (a.second / total) * (a.second / total);
    }


    double sum_s_k = 0;
    for ( auto& a: s_i )
    {
        sum_s_k += (a.second / total) * (a.second / total);
    }

    //std::cout << sum_p_ij << "\n";
    std::cout << "Rand Split: " << (sum_p_ij/sum_t_k) << "\n";
    std::cout << "Rand Merge: " << (sum_p_ij/sum_s_k) << "\n";
    std::cout << "Rand alpha: " << (sum_p_ij*2/(sum_t_k+sum_s_k)) << "\n";


    // for ( auto& a: invmap )
    // {
    //     std::vector<uint32_t> v;
    //     double sum = 0;
    //     for ( auto& b: a.second )
    //     {
    //         v.push_back(b.second);
    //         sum += b.second;
    //     }
    //     rand_split += square_sum(v);
    //     t_sq += sum;
    // }


    // for ( auto& a: map )
    // {
    //     if ( a.second.size() > 1 )
    //     {
    //         std::vector<uint32_t> v;
    //         double sum;
    //         for ( auto& b: a.second )
    //         {
    //             v.push_back(b.second);
    //             sum += b.second;
    //         }
    //         rand_merge += square_sum(v) / sum;
    //     }
    // }

    // std::cout << "Total gt: " << setg.size() << "\n";
    // std::cout << "Total ws: " << setw.size() << "\n";

    // std::cout << "Rand Split: " << rand_split << "\n";
    // std::cout << "Rand Merge: " << rand_merge << "\n";

}



int main()
{
    compare_volumes("../../../data/gt.in", "../voutall4x.out", 256);

    //return 0;

    affinity_graph_ptr<float> aff = read_affinity_graph<float>("../../../data/ws_test_256.raw",
                                                               256, 256, 256);

    // for ( float f = 0.1; f < 1.05; f += 0.1 )
    // {

    //     auto rrr = other_alg<uint32_t>(aff, f);
    //     write_volume("voutoth" + std::to_string(f) + ".out", rrr);
    //     compare_volumes("../../data/gt.in", "./voutoth" + std::to_string(f) + ".out", 256);
    // }

    // return 0;

    // volume_ptr<uint32_t> sptr = get_dummy_segmentation<uint32_t>(256, 256, 256);
    // std::vector<std::size_t> cnts(256*256*256+1);
    // std::fill_n(cnts.begin(), 256*256*256+1, 1);


    // auto rgf = get_region_graph(aff, sptr, cnts.size()-1);

    // merge_segments_with_function(sptr, rgf, cnts, limit_fn2, 100);

    // write_volume("voutall4dir.out", sptr);

    // return 0;




    std::cout << "Multiplied" << std::endl;


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

    volume_ptr<uint32_t> gt_ptr = read_volume<uint32_t>("../../../data/gt.in",
                                                        256);

    merge_segments_with_function_err(segg, gt_ptr, rg, counts, limit_fn4, 100);

    write_volume("voutall4x.out", segg);

    return 0;

    write_region_graph("voutall.rg", *rg);

    auto mt = get_merge_tree(*rg, counts.size()-1);

    write_region_graph("voutall.mt", *mt);

    return 0;


}
