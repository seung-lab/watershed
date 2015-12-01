#pragma once

#include "types.hpp"
#include <zi/disjoint_sets/disjoint_sets.hpp>

#include <queue>
#include <list>
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
