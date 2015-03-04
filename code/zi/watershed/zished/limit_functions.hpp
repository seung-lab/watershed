#pragma once

#include <cstddef>

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
