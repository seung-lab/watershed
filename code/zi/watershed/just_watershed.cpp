//#pragma once

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

#include <zi/disjoint_sets/disjoint_sets.hpp>





template< typename ID, typename F >
inline bool write_region_graph_to_file( const std::string& fname,
                                        const region_graph<ID,F>& rg )
{
    std::ofstream f(fname.c_str(), (std::ios::out | std::ios::binary) );

    if ( !f )
    {
        return false;
    }

    F* data = new F[rg.size() * 3];

    std::size_t idx = 0;

    for ( const auto& e: rg )
    {
        data[idx++] = static_cast<F>(std::get<1>(e));
        data[idx++] = static_cast<F>(std::get<2>(e));
        data[idx++] = static_cast<F>(std::get<0>(e));
    }

    f.write( reinterpret_cast<char*>(data), rg.size() * 3 * sizeof(F));

    return true;
}

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

template< typename ID, typename F >
inline region_graph_ptr<ID,F>
get_region_graph( const affinity_graph_ptr<F>& aff_ptr,
                  const volume_ptr<ID> seg_ptr,
                  std::size_t max_segid)
{

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    volume<ID>& seg = *seg_ptr;
    affinity_graph<F> aff = *aff_ptr;

    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );

    region_graph<ID,F>& rg = *rg_ptr;

    std::vector<std::map<ID,F>> edges(max_segid+1);

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][0]);
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][1]);
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][2]);
                }
            }

    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
    {
        for ( const auto& p: edges[id1] )
        {
            rg.emplace_back(p.second, id1, p.first);
        }
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::stable_sort(std::begin(rg), std::end(rg),
                     std::greater<std::tuple<F,ID,ID>>());

    std::cout << "Sorted" << std::endl;

    return rg_ptr;
}

template< typename ID >
inline volume_ptr<ID> get_dummy_segmentation( std::size_t xdim,
                                              std::size_t ydim,
                                              std::size_t zdim )
{

    volume_ptr<ID> seg_ptr( new volume<ID>(boost::extents[xdim][ydim][zdim],
                                           boost::fortran_storage_order()));
    volume<ID>& seg = *seg_ptr;

    for ( ID i = 0; i < xdim*ydim*zdim; ++i )
    {
        seg.data()[i] = i+1;
    }

    return seg_ptr;
}

template< typename AG >
inline void reverse_affinity( AG& ag )
{
    for ( auto& a: ag )
    {
        std::get<0>(a) = static_cast<float>(1) - std::get<0>(a);
    }
}

template< typename ID, typename F, typename L, typename M >
inline void merge_segments( const volume_ptr<ID>& seg_ptr,
                            const region_graph_ptr<ID,F> rg_ptr,
                            std::vector<std::size_t>& counts,
                            const L& tholds,
                            const M& lowt )
{
    zi::disjoint_sets<ID> sets(counts.size());

    typename region_graph<ID,F>::iterator rit = rg_ptr->begin();

    region_graph<ID,F>& rg  = *rg_ptr;

    for ( auto& it: tholds )
    {
        std::size_t size = static_cast<std::size_t>(it.first);
        F           thld = static_cast<F>(it.second);

        while ( (rit != rg.end()) && ( std::get<0>(*rit) > thld) )
        {
            ID s1 = sets.find_set(std::get<1>(*rit));
            ID s2 = sets.find_set(std::get<2>(*rit));

            if ( s1 != s2 && s1 && s2 )
            {
                if ( (counts[s1] < size) || (counts[s2] < size) )
                {
                    counts[s1] += counts[s2];
                    counts[s2]  = 0;
                    ID s = sets.join(s1,s2);
                    std::swap(counts[s], counts[s1]);
                }
            }
            ++rit;
        }
    }

    std::cout << "Done with merging" << std::endl;

    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

    std::size_t low = static_cast<std::size_t>(lowt);

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) && (counts[s] >= low) )
        {
            remaps[s] = next_id;
            counts[next_id] = counts[s];
            ++next_id;
        }
    }

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

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
}


template< typename ID, typename F, typename K >
inline volume_ptr<ID> other_alg( const affinity_graph_ptr<F>& aff_ptr,
                                 const K& kconst )
{
    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    volume_ptr<ID> seg_ptr = get_dummy_segmentation<ID>(xdim, ydim, zdim);
    volume<ID>&    seg     = *seg_ptr;

    region_graph_ptr<ID,F> rg_ptr = get_region_graph(aff_ptr, seg_ptr, total);
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

    ID* seg_raw = seg_ptr->data();

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

    return seg_ptr;
}


template< class C >
struct is_numeric:
    std::integral_constant<bool,
                           std::is_integral<C>::value ||
                           std::is_floating_point<C>::value> {};

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


template< typename ID, typename F, typename FN, typename M >
inline void merge_segments_with_function( const volume_ptr<ID>& seg_ptr,
                                          const region_graph_ptr<ID,F> rg_ptr,
                                          std::vector<std::size_t>& counts,
                                          const FN& func,
                                          const M& lowt )
{
    zi::disjoint_sets<ID> sets(counts.size());

    region_graph<ID,F>& rg  = *rg_ptr;

    for ( auto& it: rg )
    {
        std::size_t size = func(std::get<0>(it));

        if ( size == 0 )
        {
            break;
        }

        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        if ( s1 != s2 && s1 && s2 )
        {
            if ( (counts[s1] < size) || (counts[s2] < size) )
            {
                counts[s1] += counts[s2];
                counts[s2]  = 0;
                ID s = sets.join(s1,s2);
                std::swap(counts[s], counts[s1]);
            }
        }
    }

    std::cout << "Done with merging" << std::endl;

    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

    std::size_t low = static_cast<std::size_t>(lowt);

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) && (counts[s] >= low) )
        {
            remaps[s] = next_id;
            counts[next_id] = counts[s];
            ++next_id;
        }
    }

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

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
}


template< typename ID, typename F, typename L, typename H >
inline std::pair<volume_ptr<ID>, ID>
simple_watershed( const affinity_graph_ptr<F>& aff_ptr,
                  const L& lowv,
                  const H& highv,
                  std::vector<std::size_t>& counts )
{
    using float_type = F;
    typedef watershed_traits<ID>    traits    ;

    float_type low  = static_cast<float_type>(lowv);
    float_type high = static_cast<float_type>(highv);

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t size = xdim * ydim * zdim;

    volume_ptr<ID> seg_ptr( new volume<ID>(boost::extents[xdim][ydim][zdim],
                                           boost::fortran_storage_order()));
    affinity_graph<F>& aff = *aff_ptr;
    volume<ID>&        seg = *seg_ptr;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                ID& id = seg[x][y][z] = 0;

                F negx = (x>0) ? aff[x][y][z][0] : low;
                F negy = (y>0) ? aff[x][y][z][1] : low;
                F negz = (z>0) ? aff[x][y][z][2] : low;
                F posx = (x<(xdim-1)) ? aff[x+1][y][z][0] : low;
                F posy = (y<(ydim-1)) ? aff[x][y+1][z][1] : low;
                F posz = (z<(zdim-1)) ? aff[x][y][z+1][2] : low;

                F m = std::max({negx,negy,negz,posx,posy,posz});

                if ( m > low )
                {
                    if ( negx == m || negx >= high ) { id |= 0x01; }
                    if ( negy == m || negy >= high ) { id |= 0x02; }
                    if ( negz == m || negz >= high ) { id |= 0x04; }
                    if ( posx == m || posx >= high ) { id |= 0x08; }
                    if ( posy == m || posy >= high ) { id |= 0x10; }
                    if ( posz == m || posz >= high ) { id |= 0x20; }
                }
            }


    const std::ptrdiff_t dir[6] = { -1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim };
    const ID dirmask[6]  = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
    const ID idirmask[6] = { 0x08, 0x10, 0x20, 0x01, 0x02, 0x04 };

    // get plato escapers

    std::vector<std::ptrdiff_t> bfs;

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( !(seg_raw[idx+dir[d]] & idirmask[d]) )
                {
                    seg_raw[idx] |= 0x40;
                    bfs.push_back(idx);
                    d = 6; // break;
                }
            }
        }
    }

    std::size_t bfs_index = 0;

    while ( bfs_index < bfs.size() )
    {
        std::ptrdiff_t idx = bfs[bfs_index];

        ID to_set = 0;

        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( seg_raw[idx+dir[d]] & idirmask[d] )
                {
                    if ( !( seg_raw[idx+dir[d]] & 0x40 ) )
                    {
                        bfs.push_back(idx+dir[d]);
                        seg_raw[idx+dir[d]] |= 0x40;
                    }
                }
                else
                {
                    to_set = dirmask[d];
                }
            }
        }
        seg_raw[idx] = to_set;
        ++bfs_index;
    }

    bfs.clear();

    //std::vector<std::size_t> counts({0});

    counts.resize(1);
    counts[0] = 0;

    ID next_id = static_cast<ID>(1);

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        if ( seg_raw[idx] == 0 )
        {
            seg_raw[idx] |= traits::high_bit;
            ++counts[0];
        }

        if ( !( seg_raw[idx] & traits::high_bit ) && seg_raw[idx] )
        {
            bfs.push_back(idx);
            bfs_index = 0;
            seg_raw[idx] |= 0x40;

            while ( bfs_index < bfs.size() )
            {
                std::ptrdiff_t me = bfs[bfs_index];

                for ( std::ptrdiff_t d = 0; d < 6; ++d )
                {
                    if ( seg_raw[me] & dirmask[d] )
                    {
                        std::ptrdiff_t him = me + dir[d];
                        if ( seg_raw[him] & traits::high_bit )
                        {
                            counts[ seg_raw[him] & ~traits::high_bit ]
                                += bfs.size();

                            for ( auto& it: bfs )
                            {
                                seg_raw[it] = seg_raw[him];
                            }

                            bfs.clear();
                            d = 6; // break
                        }
                        else if ( !( seg_raw[him] & 0x40 ) )
                        {
                            seg_raw[him] |= 0x40;
                            bfs.push_back( him );

                        }
                    }
                }
                ++bfs_index;
            }

            if ( bfs.size() )
            {
                counts.push_back( bfs.size() );
                for ( auto& it: bfs )
                {
                    seg_raw[it] = traits::high_bit | next_id;
                }
                ++next_id;
                bfs.clear();
            }
        }
    }

    std::cout << "found: " << (next_id-1) << " components\n";

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        seg_raw[idx] &= traits::mask;
    }

    return std::make_pair(seg_ptr, next_id-1);
}


template< typename ID, typename F, typename L, typename H >
inline std::tuple< volume_ptr<ID>, std::vector<std::size_t> >
watershed( const affinity_graph_ptr<F>& aff_ptr, const L& lowv, const H& highv )
{
    using affinity_t = F;
    using id_t       = ID;
    using traits     = watershed_traits<id_t>;

    affinity_t low  = static_cast<affinity_t>(lowv);
    affinity_t high = static_cast<affinity_t>(highv);

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t size = xdim * ydim * zdim;

    std::tuple< volume_ptr<id_t>, std::vector<std::size_t> > result
        ( volume_ptr<id_t>( new volume<id_t>(boost::extents[xdim][ydim][zdim],
                                           boost::fortran_storage_order())),
          std::vector<std::size_t>(1) );

    auto& counts = std::get<1>(result);
    counts[0] = 0;

    affinity_graph<F>& aff = *aff_ptr;
    volume<id_t>&      seg = *std::get<0>(result);

    id_t* seg_raw = seg.data();

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                id_t& id = seg[x][y][z] = 0;

                F negx = (x>0) ? aff[x][y][z][0] : low;
                F negy = (y>0) ? aff[x][y][z][1] : low;
                F negz = (z>0) ? aff[x][y][z][2] : low;
                F posx = (x<(xdim-1)) ? aff[x+1][y][z][0] : low;
                F posy = (y<(ydim-1)) ? aff[x][y+1][z][1] : low;
                F posz = (z<(zdim-1)) ? aff[x][y][z+1][2] : low;

                F m = std::max({negx,negy,negz,posx,posy,posz});

                if ( m > low )
                {
                    if ( negx == m || negx >= high ) { id |= 0x01; }
                    if ( negy == m || negy >= high ) { id |= 0x02; }
                    if ( negz == m || negz >= high ) { id |= 0x04; }
                    if ( posx == m || posx >= high ) { id |= 0x08; }
                    if ( posy == m || posy >= high ) { id |= 0x10; }
                    if ( posz == m || posz >= high ) { id |= 0x20; }
                }
            }


    const std::ptrdiff_t dir[6] = { -1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim };
    const id_t dirmask[6]  = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
    const id_t idirmask[6] = { 0x08, 0x10, 0x20, 0x01, 0x02, 0x04 };

    // get plato corners

    std::vector<std::ptrdiff_t> bfs;

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( !(seg_raw[idx+dir[d]] & idirmask[d]) )
                {
                    seg_raw[idx] |= 0x40;
                    bfs.push_back(idx);
                    d = 6; // break;
                }
            }
        }
    }

    // divide the plateaus

    std::size_t bfs_index = 0;

    while ( bfs_index < bfs.size() )
    {
        std::ptrdiff_t idx = bfs[bfs_index];

        id_t to_set = 0;

        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( seg_raw[idx+dir[d]] & idirmask[d] )
                {
                    if ( !( seg_raw[idx+dir[d]] & 0x40 ) )
                    {
                        bfs.push_back(idx+dir[d]);
                        seg_raw[idx+dir[d]] |= 0x40;
                    }
                }
                else
                {
                    to_set = dirmask[d];
                }
            }
        }
        seg_raw[idx] = to_set;
        ++bfs_index;
    }

    bfs.clear();

    // main watershed logic

    id_t next_id = 1;

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        if ( seg_raw[idx] == 0 )
        {
            seg_raw[idx] |= traits::high_bit;
            ++counts[0];
        }

        if ( !( seg_raw[idx] & traits::high_bit ) && seg_raw[idx] )
        {
            bfs.push_back(idx);
            bfs_index = 0;
            seg_raw[idx] |= 0x40;

            while ( bfs_index < bfs.size() )
            {
                std::ptrdiff_t me = bfs[bfs_index];

                for ( std::ptrdiff_t d = 0; d < 6; ++d )
                {
                    if ( seg_raw[me] & dirmask[d] )
                    {
                        std::ptrdiff_t him = me + dir[d];
                        if ( seg_raw[him] & traits::high_bit )
                        {
                            counts[ seg_raw[him] & ~traits::high_bit ]
                                += bfs.size();

                            for ( auto& it: bfs )
                            {
                                seg_raw[it] = seg_raw[him];
                            }

                            bfs.clear();
                            d = 6; // break
                        }
                        else if ( !( seg_raw[him] & 0x40 ) )
                        {
                            seg_raw[him] |= 0x40;
                            bfs.push_back( him );

                        }
                    }
                }
                ++bfs_index;
            }

            if ( bfs.size() )
            {
                counts.push_back( bfs.size() );
                for ( auto& it: bfs )
                {
                    seg_raw[it] = traits::high_bit | next_id;
                }
                ++next_id;
                bfs.clear();
            }
        }
    }

    std::cout << "found: " << (next_id-1) << " components\n";

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        seg_raw[idx] &= traits::mask;
    }

    return result;
}



template< typename T >
affinity_graph_ptr<T> mult_aff( const affinity_graph_ptr<T>& aff, int n )
{

    std::ptrdiff_t xdim = aff->shape()[0];
    std::ptrdiff_t ydim = aff->shape()[1];
    std::ptrdiff_t zdim = aff->shape()[2];

    affinity_graph_ptr<T> r(new affinity_graph<T>
                            (boost::extents[xdim*n][ydim*n][zdim*n][3],
                             boost::fortran_storage_order()));

    for ( std::ptrdiff_t z = 0; z < n; ++z )
        for ( std::ptrdiff_t y = 0; y < n; ++y )
            for ( std::ptrdiff_t x = 0; x < n; ++x )
            {

                typedef boost::multi_array_types::index_range range;

                (*r)[boost::indices[range(x*xdim,(x+1)*xdim)]
                     [range(y*ydim,(y+1)*ydim)]
                     [range(z*zdim,(z+1)*zdim)]
                     [range(0,3)]] = *aff;

            }

    return r;

}

std::size_t limit_fn( float v )
{
    //return 50;

    // size threshold based on affinity
    if ( v > 1 )
    {
        return 2000;
    }

    if ( v < 0.3 )
    {
        return 0;
    }

    if ( v < 0.5 )
    {
        return 50;
    }

    v *= 10;


    return static_cast<std::size_t>(50+v*v*v);
}


std::size_t limit_fn2( float v )
{
    if ( v < 0.3 ) return 0;

    if ( v > 0.98 ) return 5000;
    if ( v > 0.97 ) return 2000;
    if ( v > 0.96 ) return 1500;
    if ( v > 0.95 ) return 500;

    if ( v > 0.5 ) return 250;
    //return 500;
    return 100;

    // size threshold based on affinity
    if ( v > 1 )
    {
        return 2000;
    }

    if ( v < 0.3 )
    {
        return 0;
    }

    if ( v < 0.5 )
    {
        return 150;
    }

    v *= 10;


    return static_cast<std::size_t>(50+v*v*v);
}


std::size_t limit_fn4( float v )
{
    if ( v < 0.3 ) return 0;

    return 250 + 50000.0 * (v-0.3)*(v-0.3) / 0.7 / 0.7;


    if ( v > 0.98 ) return 5000;
    if ( v > 0.97 ) return 2000;
    if ( v > 0.96 ) return 1500;
    if ( v > 0.95 ) return 500;

    if ( v > 0.5 ) return 250;
    //return 500;
    return 100;

    // size threshold based on affinity
    if ( v > 1 )
    {
        return 2000;
    }

    if ( v < 0.3 )
    {
        return 0;
    }

    if ( v < 0.5 )
    {
        return 150;
    }

    v *= 10;


    return static_cast<std::size_t>(50+v*v*v);
}



std::size_t limit_fn3( float v )
{
    if ( v < 0.3 ) return 0;
    return 100;

    if ( v > 0.98 ) return 5000;
    if ( v > 0.97 ) return 2000;
    if ( v > 0.96 ) return 1500;
    if ( v > 0.95 ) return 500;

    if ( v > 0.5 ) return 250;
    //return 500;
    return 100;

    // size threshold based on affinity
    if ( v > 1 )
    {
        return 2000;
    }

    if ( v < 0.3 )
    {
        return 0;
    }

    if ( v < 0.5 )
    {
        return 150;
    }

    v *= 10;


    return static_cast<std::size_t>(50+v*v*v);
}


template< typename ID, typename F >
inline volume_ptr<ID>
get_colormap( const volume_ptr<ID>& orig,
              const region_graph<ID,F>& rg,
              std::size_t max_segid )
{

    std::ptrdiff_t xdim = orig->shape()[0];
    std::ptrdiff_t ydim = orig->shape()[1];
    std::ptrdiff_t zdim = orig->shape()[2];

    std::ptrdiff_t size = xdim * ydim * zdim;

    volume_ptr<ID> seg_ptr( new volume<ID>(boost::extents[xdim][ydim][zdim],
                                           boost::fortran_storage_order()));

    volume<ID>& seg = *seg_ptr;
    seg = *orig;

    std::vector<ID>            colormap(max_segid+1);
    std::vector<std::list<ID>> edges(max_segid+1);
    std::vector<std::set<ID>>  taken(max_segid+1);

    for ( const auto& e: rg )
    {
        ID v1 = std::get<1>(e);
        ID v2 = std::get<2>(e);

        if ( v1 != v2 && v1 && v2 )
        {
            edges[v1].push_back(v2);
            edges[v2].push_back(v1);

            //std::cout << edges[v1].size() << ' ' << edges[v2].size() << '\n';
        }
    }

    return seg_ptr;
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
    compare_volumes("../../data/gt.in", "./voutall4x.out", 256);

    return 0;

    affinity_graph_ptr<float> aff = read_affinity_graph_from_file<float>("../../data/ws_test_256.raw",
                                                                         256, 256, 256);

    // for ( float f = 0.1; f < 1.05; f += 0.1 )
    // {

    //     auto rrr = other_alg<uint32_t>(aff, f);
    //     write_volume_to_file("voutoth" + std::to_string(f) + ".out", rrr);
    //     compare_volumes("../../data/gt.in", "./voutoth" + std::to_string(f) + ".out", 256);
    // }

    // return 0;

    // volume_ptr<uint32_t> sptr = get_dummy_segmentation<uint32_t>(256, 256, 256);
    // std::vector<std::size_t> cnts(256*256*256+1);
    // std::fill_n(cnts.begin(), 256*256*256+1, 1);


    // auto rgf = get_region_graph(aff, sptr, cnts.size()-1);

    // merge_segments_with_function(sptr, rgf, cnts, limit_fn2, 100);

    // write_volume_to_file("voutall4dir.out", sptr);

    // return 0;




    std::cout << "Multiplied" << std::endl;


    volume_ptr<uint32_t>     segg  ;
    std::vector<std::size_t> counts;

    auto seg = simple_watershed<uint32_t>(aff, -0.00001, 1.2, counts);

    std::tie(segg, counts) = watershed<uint32_t>(aff, -1, 2);

    write_volume_to_file("voutraw.out", segg);

    for ( float low = 0.01; low < 0.051; low += 0.01 )
    {
        for ( float high = 0.998; high > 0.989; high -= 0.002 )
        {
//            std::tie(segg, counts) = watershed<uint32_t>(aff, low, high);
//            write_volume_to_file("vout." + std::to_string(low) + "." +
//                                 std::to_string(high) + ".out", segg);
        }
    }

    std::tie(segg, counts) = watershed<uint32_t>(aff, 0.5, 2);

    write_volume_to_file("voutmax.out", segg);

    std::tie(segg, counts) = watershed<uint32_t>(aff, 0.3, 0.99);

    write_volume_to_file("voutminmax.out", segg);


//    return 0;

    // auto rg = get_region_graph(aff, segg, counts.size()-1);

    // //yet_another_watershed(segg, rg, counts, 0.3);

    // //write_volume_to_file("voutanouther.out", segg);


    // merge_segments_with_function(segg, rg, counts, limit_fn3, 100);

    // write_volume_to_file("voutdo.out", segg);




    auto rg = get_region_graph(aff, segg, counts.size()-1);

    //yet_another_watershed(segg, rg, counts, 0.3);

    //write_volume_to_file("voutanouther.out", segg);


    merge_segments_with_function(segg, rg, counts, limit_fn4, 100);

    write_volume_to_file("voutall4x.out", segg);

    return 0;

    write_region_graph_to_file("voutall.rg", *rg);

    auto mt = get_merge_tree(*rg, counts.size()-1);

    write_region_graph_to_file("voutall.mt", *mt);

    auto cmap = get_colormap(segg, *rg, counts.size()-1);

    return 0;

    // single likage clusteroing

    //std::cout << (std::chrono::system_clock::now() - start)
    //          << " simple ws done" << std::endl;




    //std::cout << (std::chrono::system_clock::now() - start)
    //          << " region graph created" << std::endl;

    std::list<std::pair<std::size_t, float>> rules;

    // rules.emplace_back(1500,0.98);
    // rules.emplace_back(500,0.95);
    // rules.emplace_back(250,0.91);
    // rules.emplace_back(100,0.89);
    // rules.emplace_back(25,0.0);


    //merge_segments_with_function(seg.first, rg, counts, limit_fn, 25);

    yet_another_watershed(seg.first, rg, counts, 0.1);
    yet_another_watershed(seg.first, rg, counts, 0.2);
    yet_another_watershed(seg.first, rg, counts, 0.3);
    yet_another_watershed(seg.first, rg, counts, 0.4);
    //yet_another_watershed(seg.first, rg, counts);

    //merge_segments_with_function(seg.first, rg, counts, limit_fn, 25);

    //std::cout << (std::chrono::system_clock::now() - start)
    //          << " segments merged" << std::endl;


    std::cout << "here" << std::endl;



    //int x;
    //std::cin >> x;

    //while (1);

}
