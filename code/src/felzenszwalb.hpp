#pragma once

#include "types.hpp"
#include "utils.hpp"
#include "region_graph.hpp"

#include <zi/disjoint_sets/disjoint_sets.hpp>


template< typename ID, typename F, typename K >
inline volume_ptr<ID> felzenszwalb( const affinity_graph_ptr<F>& aff_ptr,
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
    std::vector<std::size_t>& counts = std::get<1>(dummy);

    for ( auto& it: rg )
    {
        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        F  w  = static_cast<F>(1) - std::get<0>(it);

        if ( s1 != s2 )
        {
            F v1 = int_val[s1] + k / counts[s1];
            F v2 = int_val[s2] + k / counts[s2];

            if ( w < std::min(v1, v2) )
            {
                counts[s1] += counts[s2];
                counts[s2] = 0;

                int_val[s1] = w;
                int_val[s2] = 0;

                ID s = sets.join(s1,s2);

                std::swap(counts[s1], counts[s]);
                std::swap(int_val[s1], int_val[s]);
            }
        }
    }

    std::cout << "Done with merging" << std::endl;

    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

    for ( ID id = 0; id < counts.size(); ++id )
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



template< typename ID, typename F, typename K >
inline volume_ptr<ID> felzenszwalbws( const affinity_graph_ptr<F>& aff_ptr,
                                      volume_ptr<ID> seg_ptr,
                                      region_graph_ptr<ID,F> rg_ptr,
                                      std::vector<std::size_t>& counts,
                                      const K& kconst )
{

    //auto dummy          = get_dummy_segmentation<ID>(xdim, ydim, zdim);
    volume<ID>&    seg  = *seg_ptr; //std::get<0>(dummy);

    //region_graph_ptr<ID,F> rg_ptr = get_region_graph(aff_ptr,
    //                                                   std::get<0>(dummy), total);
    region_graph<ID,F>& rg  = *rg_ptr;

    zi::disjoint_sets<ID> sets(counts.size());

    F k = kconst;

    std::vector<F> int_val(counts.size());
    //std::vector<std::size_t>& counts = std::get<1>(dummy);

    for ( auto& it: rg )
    {
        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        F  w  = static_cast<F>(1) - std::get<0>(it);

        if ( s1 != s2 )
        {
            F v1 = int_val[s1] + k / counts[s1];
            F v2 = int_val[s2] + k / counts[s2];

            if ( w < std::min(v1, v2) )
            {
                counts[s1] += counts[s2];
                counts[s2] = 0;

                int_val[s1] = w;
                int_val[s2] = 0;

                ID s = sets.join(s1,s2);

                std::swap(counts[s1], counts[s]);
                std::swap(int_val[s1], int_val[s]);
            }
        }
    }

    std::cout << "Done with merging" << std::endl;

    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) )
        {
            remaps[s] = next_id;
            ++next_id;
        }
    }

    ID* seg_raw = seg.data();

    for ( std::ptrdiff_t idx = 0; idx < seg.num_elements(); ++idx )
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

    return seg_ptr;
}





template< typename ID, typename F, typename K >
inline std::vector<double>
felzenszwalb_err( const volume_ptr<ID>& gt_ptr,
                  const affinity_graph_ptr<F>& aff_ptr,
                  const K& kconst )
{
    std::vector<double> ret;

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    auto dummy          = get_dummy_segmentation<ID>(xdim, ydim, zdim);
    volume<ID>&    seg  = *std::get<0>(dummy);
    volume<ID>&    gt   = *gt_ptr;

    region_graph_ptr<ID,F> rg_ptr = get_region_graph(aff_ptr,
                                                     std::get<0>(dummy), total);
    region_graph<ID,F>& rg  = *rg_ptr;

    zi::disjoint_sets<ID> sets(total+1);

    std::size_t tot = 0;

    std::vector<std::size_t> s_i(total+1);
    std::map<ID,std::size_t> t_j;

    std::vector<std::map<ID,std::size_t>> p_ij(total+1);


    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                uint32_t sgv = seg[x][y][z];
                uint32_t gtv =  gt[x][y][z];

                if ( gtv )
                {
                    tot += 1;

                    ++p_ij[sgv][gtv];

                    ++s_i[sgv];
                    ++t_j[gtv];
                }
            }

    double sum_p_ij = 0;
    for ( auto& a: p_ij )
    {
        for ( auto& b: a )
        {
            sum_p_ij += static_cast<double>(b.second) * b.second;
        }
    }

    double sum_t_k = 0;
    for ( auto& a: t_j )
    {
        sum_t_k += static_cast<double>(a.second) * a.second;
    }


    double sum_s_k = 0;
    for ( auto& a: s_i )
    {
        sum_s_k += static_cast<double>(a) * a;
    }

    ret.push_back(sum_p_ij/sum_t_k);
    ret.push_back(sum_p_ij/sum_s_k);

    F k = kconst;

    std::vector<F> int_val(total+1);
    std::vector<std::size_t>& counts = std::get<1>(dummy);

    for ( auto& it: rg )
    {
        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        F  w  = static_cast<F>(1) - std::get<0>(it);

        if ( s1 != s2 )
        {
            F v1 = int_val[s1] + k / counts[s1];
            F v2 = int_val[s2] + k / counts[s2];

            if ( w < std::min(v1, v2) )
            {
                counts[s1] += counts[s2];
                counts[s2] = 0;

                int_val[s1] = w;
                int_val[s2] = 0;

                sum_s_k -= static_cast<double>(s_i[s1]) * s_i[s1];
                sum_s_k -= static_cast<double>(s_i[s2]) * s_i[s2];

                s_i[s1] += s_i[s2];
                s_i[s2]  = 0;

                sum_s_k += static_cast<double>(s_i[s1]) * s_i[s1];


                for ( auto& b: p_ij[s1] )
                {
                    sum_p_ij -= static_cast<double>(b.second) * b.second;
                }

                for ( auto& b: p_ij[s2] )
                {
                    sum_p_ij -= static_cast<double>(b.second) * b.second;
                    p_ij[s1][b.first] += b.second;
                }

                for ( auto& b: p_ij[s1] )
                {
                    sum_p_ij += static_cast<double>(b.second) * b.second;
                }

                p_ij[s2].clear();

                ret.push_back(sum_p_ij/sum_t_k);
                ret.push_back(sum_p_ij/sum_s_k);

                ID s = sets.join(s1,s2);

                std::swap(counts[s1], counts[s]);
                std::swap(int_val[s1], int_val[s]);
                std::swap(s_i[s], s_i[s1]);
                std::swap(p_ij[s], p_ij[s1]);            }
        }
    }

    std::cout << "Done with merging" << std::endl;

    return ret;
}
