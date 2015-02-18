//#pragma once

#include <boost/multi_array.hpp>
#include <boost/multi_array/types.hpp>
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
#include <vector>
#include <chrono>
#include <set>

#include <zi/disjoint_sets/disjoint_sets.hpp>

template < typename T >
void clear_container( T& c )
{
    T n;
    n.swap(c);
}

template < typename T > struct watershed_traits;

template <> struct watershed_traits<uint16_t>
{
    static const uint16_t high_bit = 0x8000;
    static const uint16_t mask     = 0x7FFF;
    static const uint16_t visited  = 0x1000;
    static const uint16_t dir_mask = 0x007F;
};

template <> struct watershed_traits<uint32_t>
{
    static const uint32_t high_bit = 0x80000000;
    static const uint32_t mask     = 0x7FFFFFFF;
    static const uint32_t visited  = 0x00001000;
    static const uint32_t dir_mask = 0x0000007F;
};

template <> struct watershed_traits<uint64_t>
{
    static const uint64_t high_bit = 0x8000000000000000LL;
    static const uint64_t mask     = 0x7FFFFFFFFFFFFFFFLL;
    static const uint64_t visited  = 0x0000000000001000LL;
    static const uint64_t dir_mask = 0x000000000000007FLL;
};

template <> struct watershed_traits<int16_t>
: watershed_traits<uint16_t> {};

template <> struct watershed_traits<int32_t>
: watershed_traits<uint32_t> {};

template <> struct watershed_traits<int64_t>
: watershed_traits<uint64_t> {};


template < typename T >
using volume = boost::multi_array<T,3>;

template < typename T >
using affinity_graph = boost::multi_array<T,4>;

template < typename T >
using volume_ptr = std::shared_ptr<volume<T>>;

template < typename T >
using affinity_graph_ptr = std::shared_ptr<affinity_graph<T>>;


template < typename T >
inline bool read_from_file( const std::string& fname, T* data, std::size_t n )
{
    FILE* f = std::fopen(fname.c_str(), "rbXS");

    if ( !f )
    {
        return false;
    }

    std::size_t nread = std::fread(data, sizeof(T), n, f);

    std::fclose(f);

    return nread == n;
}

template < typename T >
inline affinity_graph_ptr<T>
read_affinity_graph_from_file( const std::string& fname,
                               std::size_t xsize,
                               std::size_t ysize,
                               std::size_t zsize )
{
    affinity_graph_ptr<T> aff(new affinity_graph<T>
                              (boost::extents[xsize][ysize][zsize][3],
                               boost::fortran_storage_order()));

    if ( !read_from_file(fname, aff->data(), xsize*ysize*zsize*3) )
    {
        throw 0;
    }

    return aff;
}


template < typename T >
inline bool
write_volume_to_file( const std::string& fname,
                      const volume_ptr<T>& vol )
{
    std::ofstream f(fname.c_str(), (std::ios::out | std::ios::binary) );

    if ( !f )
    {
        return false;
    }

    f.write( reinterpret_cast<char*>(vol->data()),
             vol->shape()[0] * vol->shape()[1] * vol->shape()[2] * sizeof(T));

    return true;
}

template< typename ID, typename F >
using region_graph = std::vector<std::tuple<F,ID,ID>>;

template< typename ID, typename F >
using region_graph_ptr = std::shared_ptr<region_graph<ID,F>>;


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
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z]
                     && (seg[x][y][z] != seg[x-1][y][z] ) )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][0]);
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z]
                     && (seg[x][y][z] != seg[x][y-1][z] ) )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][1]);
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1]
                     && (seg[x][y][z] != seg[x][y][z-1] ) )
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

    std::sort(std::begin(rg), std::end(rg),
              std::greater<std::tuple<F,ID,ID>>());

    std::cout << "Sorted" << std::endl;

    return rg_ptr;
}


template< typename ID, typename F >
inline region_graph_ptr<ID,F>
get_avg_region_graph( const affinity_graph_ptr<F>& aff_ptr,
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
    std::vector<std::map<ID,F>> numed(max_segid+1);


    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z]
                     && (seg[x][y][z] != seg[x-1][y][z] ) )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr += aff[x][y][z][0];
                    numed[mm.first][mm.second] += 1;
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z]
                     && (seg[x][y][z] != seg[x][y-1][z] ) )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr += aff[x][y][z][1];
                    numed[mm.first][mm.second] += 1;
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1]
                     && (seg[x][y][z] != seg[x][y][z-1] ) )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                    F& curr = edges[mm.first][mm.second];
                    curr += aff[x][y][z][2];
                    numed[mm.first][mm.second] += 1;
                }
            }

    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
    {
        for ( const auto& p: edges[id1] )
        {
            rg.emplace_back(p.second / numed[id1][p.first], id1, p.first);
        }
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::sort(std::begin(rg), std::end(rg),
              std::greater<std::tuple<F,ID,ID>>());

    std::cout << "Sorted" << std::endl;

    return rg_ptr;
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

    std::set<std::pair<ID,ID>> pairs;

    for ( auto& it: rg )
    {
        ID s1 = remaps[sets.find_set(std::get<1>(it))];
        ID s2 = remaps[sets.find_set(std::get<2>(it))];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            if ( pairs.count(mm) == 0 )
            {
                new_rg.emplace_back(std::get<0>(it), mm.first, mm.second);
                pairs.insert(mm);
            }
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

    if ( v < 0.2 )
    {
        return 0;
    }

    if ( v < 0.5 ) return 200;
    if ( v < 0.8 ) return 500;
    if ( v < 0.95) return 1000;

    return 2000;

    // if ( v < 0.5 )
    // {
    //     return 50;
    // }

    // v *= 10;


    // return static_cast<std::size_t>(50+v*v*v);
}


std::size_t limit_fn2( float v )
{
    if ( v < 0.95 ) return 0;
    return 10000;
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


int main()
{

    //auto start = std::chrono::system_clock::now();


    //std::cout << std::hex << watershed_traits<uint32_t>::high_bit << std::endl;

    affinity_graph_ptr<float> aff = read_affinity_graph_from_file<float>("../../data/ws_test_250.raw",
                                                                         250, 250, 250);

    //std::cout << (std::chrono::system_clock::now() - start)
    //          << " data loaded" << std::endl;


    //std::fill_n(aff->data(), 160*160*160*3, 0);

    //  (*aff)[10][10][10][0] = 0.5;
    // (*aff)[10][10][150][0] = 0.5;
    // (*aff)[10][100][100][0] = 0.5;


    //aff = mult_aff(aff, 3);

    std::cout << "Multiplied" << std::endl;


    volume_ptr<uint32_t>     segg  ;
    std::vector<std::size_t> counts;

    auto seg = simple_watershed<uint32_t>(aff, 0.1, 0.99, counts);

    std::tie(segg, counts) = watershed<uint32_t>(aff, 0.1, 0.99);



    //return 0;

    auto rg = get_region_graph(aff, segg, counts.size()-1);

    //merge_segments_with_function(seg.first, rg, counts, limit_fn2, 0);

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

    //yet_another_watershed(segg, rg, counts, 0.2);

    //rg = get_avg_region_graph(aff, segg, counts.size()-1);

    //yet_another_watershed(segg, rg, counts, 0.2);

    //rg = get_avg_region_graph(aff, segg, counts.size()-1);

    //yet_another_watershed(segg, rg, counts, 0.2);

    //rg = get_avg_region_graph(aff, segg, counts.size()-1);

    //yet_another_watershed(segg, rg, counts, 0.2);

    //rg = get_avg_region_graph(aff, segg, counts.size()-1);

    merge_segments_with_function(segg, rg, counts, limit_fn, 300);


    rg = get_region_graph(aff, segg, counts.size()-1);

    // rg = get_avg_region_graph(aff, segg, counts.size()-1);

    // yet_another_watershed(segg, rg, counts, 0.0);
    // rg = get_avg_region_graph(aff, segg, counts.size()-1);

    // yet_another_watershed(segg, rg, counts, 0.0);
    // rg = get_avg_region_graph(aff, segg, counts.size()-1);

    // yet_another_watershed(segg, rg, counts, 0.0);
    // rg = get_avg_region_graph(aff, segg, counts.size()-1);

    // yet_another_watershed(segg, rg, counts, 0.0);

    //yet_another_watershed(segg, rg, counts, 0.0);

    //yet_another_watershed(segg, rg, counts, 0.0);

    write_volume_to_file("voutraw.out", segg);

    //yet_another_watershed(seg.first, rg, counts);

    //merge_segments_with_function(seg.first, rg, counts, limit_fn, 25);

    //std::cout << (std::chrono::system_clock::now() - start)
    //          << " segments merged" << std::endl;


    std::cout << "here" << std::endl;



    //int x;
    //std::cin >> x;

    //while (1);

}
