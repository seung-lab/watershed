#define NPY_DEPRECATED_API NPY_1_7_API_VERSION
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numpy.hpp>
#include <boost/multi_array.hpp>
#include <boost/multi_array/types.hpp>
#include <boost/multi_array/extent_gen.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <vector>
#include <cstdlib>
#include <string>
#include <memory>
#include <cstdint>
#include <assert.h>
#include <tuple>
#include <chrono>

#include "basic_watershed.hpp"
#include "region_graph.hpp"
#include "limit_functions.hpp"
#include "agglomeration.hpp"
#include "merge_tree.hpp"
#include "watershed_full.hpp"
#include "watershed_config.hpp"

namespace b = boost;
namespace bp = boost::python;
namespace bdma = boost::detail::multi_array;
namespace np = boost::numpy;
namespace chrono = std::chrono;

//Loops over a bp::dict, and returns a vector of key strings
// used to define loop in pyopt_to_znnopt
std::vector<std::string> extract_dict_keys( bp::dict const & layer_dict )
{
    std::vector<std::string> res;
    std::size_t num_keys = bp::extract<std::size_t>( layer_dict.attr("__len__")() );

    for (std::size_t i=0; i < num_keys; i++ )
    {
        //fetching the keys through the bp interface...one at a time...
        res.push_back( bp::extract<std::string>(layer_dict.attr("keys")()[i]) );
    }

    return res;
}

/**
 * Calls watershed on the input affinity graph
 * @param np_affinity_graph this is an affinity graph in regular c order (c-z-y-x)
 * @param lowv this is the low threshold when running watershed
 * @param highv this is hte high threshold when running watershed
 * @return a tuple of the segmentation volume (ndarray) and the counts(list) per segmentation.
 *  WARNING: segmentation volume is returned in FORTRAN ORDER (x-y-z)
 *  TODO investigate to change watershed to use c order instead!
 */
bp::tuple py_watershed(np::ndarray& np_affinity_graph, float lowv, float highv)
{
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    std::size_t channel_size = np_affinity_graph.shape(0);
    std::size_t z_size = np_affinity_graph.shape(1);
    std::size_t y_size = np_affinity_graph.shape(2);
    std::size_t x_size = np_affinity_graph.shape(3);
    std::cout << "Begin Converting input data to boost multi_array" << std::endl;
    affinity_graph_ptr<float> b_affinity_graph_ptr(new b::multi_array_ref<float, 4> (
                reinterpret_cast<float *>(np_affinity_graph.get_data()),
                               boost::extents[x_size][y_size][z_size][channel_size],
                               boost::fortran_storage_order()));
    std::cout << "Finish Converting affinity graph to boost multi_array in " <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Original dimensions :: " << std::endl <<
        channel_size << "," << z_size << "," << y_size << "," << x_size <<
        " dtype= " << bp::extract<char const *> (bp::str(np_affinity_graph.get_dtype())) <<
        std::endl;

    std::cout << "Begin calling watershed" << std::endl;

    // Run watershed here!
    std::tuple<volume_ptr<uint32_t>, std::vector<std::size_t>> tuple = watershed<uint32_t>(b_affinity_graph_ptr, lowv, highv);
    std::cout << "Finished calling watershed in " << 
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting data back to ndarray";
    b::multi_array<uint32_t, 3> b_segmentation_volume = *(std::get<0>(tuple));
    std::vector<size_t> b_segmentation_counts = std::get<1>(tuple);

    np::ndarray np_segmentation_volume = np::from_data(
            reinterpret_cast<void *> (b_segmentation_volume.data()),
            np::dtype::get_builtin<uint32_t>(),
            bp::make_tuple(x_size, y_size, z_size),
            bp::make_tuple(y_size * z_size * sizeof(uint32_t), 
                z_size * sizeof(uint32_t),
                sizeof(uint32_t)),
            bp::object());
    std::cout << "Finished converting segmentation volume back to ndarray in " <<
        chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
         << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    bp::list bp_segmentation_counts;
    for (size_t i = 0; i < b_segmentation_counts.size(); ++i)
    {
        bp_segmentation_counts.append(b_segmentation_counts[i]);
    }
    std::cout << "Finished converting segmentation counts back to ndarray in " <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    // COPY the data here so that the memory doesn't get deallocated
    return bp::make_tuple(np_segmentation_volume.copy(), bp_segmentation_counts);
}


/**
 * Generates a region graph of the segmentation data.
 * @param np_affinity_graph this is an affinity graph in regular c order (c-z-y-x)
 * @param np_volume this is the segmentation labels of the volume in FORTRAN order (x-y-z-c)
 * @param max_segid this is the maximum value of the segmentation label
 * @return the region graph as a list of tuples of weights and nodes (Weight, Edge1Id, Edge2Id)
 */
bp::list py_region_graph(
        np::ndarray& np_affinity_graph, np::ndarray& np_volume, std::size_t max_segid)
{
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    std::size_t channel_size = np_affinity_graph.shape(0);
    std::size_t z_size = np_affinity_graph.shape(1);
    std::size_t y_size = np_affinity_graph.shape(2);
    std::size_t x_size = np_affinity_graph.shape(3);
    std::cout << "Begin Converting input data to boost multi_array"  << std::endl;

    affinity_graph_ptr<float> b_affinityGraph_ptr(new b::multi_array_ref<float, 4> (
                reinterpret_cast<float *>(np_affinity_graph.get_data()),
                               boost::extents[x_size][y_size][z_size][channel_size],
                               boost::fortran_storage_order()));
    std::cout << "Finished converting affinity_graph to boost multi_array "  <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    volume_ptr<uint32_t> b_segmentation_ptr(new b::multi_array_ref<uint32_t, 3> (
                reinterpret_cast<uint32_t *>(np_volume.get_data()),
                               boost::extents[x_size][y_size][z_size]
                               ));
    std::cout << "Finished converting segmentation_volume to boost multi_array in " <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin calling region_graph" << std::endl;
    region_graph_ptr<uint32_t, float> region_graph_ptr = get_region_graph(b_affinityGraph_ptr, b_segmentation_ptr, max_segid);
    std::cout << "Finished calling region_graph in " << 
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting output data in " << std::endl;
    bp::list region_graph_edges;
    for (auto e : *region_graph_ptr )
    {
        float weight = std::get<0>(e);
        uint32_t node1 = std::get<1>(e);
        uint32_t node2 = std::get<2>(e);
        bp::tuple edge = bp::make_tuple(weight, node1, node2);
        region_graph_edges.append(edge);
    }
    std::cout << "Finished converting output region_graph_edges in " <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();
    return region_graph_edges;
}

/**
 * Merge segments based on function currently uses const_above_threshold()
 * IN PLACE MODIFICATION of segmentation volume and region_graph_edges
 * @param np_segmentation_volume this is the segmentation labels of the 
 *  volume in FORTRAN order (x-y-z)
 * @param bp_region_graph_edges list of tuples as edges of the region graph in (weight , id, id)
 * @param bp_counts python list of counts per segmentation label
 * @param lowt lower size threshold for merge
 * @param min_mult coefficient for upper size threshold thold
 * @param thold upper size threshold for merge
 */
void py_merge_segments(
        np::ndarray& np_segmentation_volume, bp::list& bp_region_graph_edges,
        bp::list &bp_counts, size_t lowt, float min_mult, size_t thold)
{
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    std::size_t x_size = np_segmentation_volume.shape(0);
    std::size_t y_size = np_segmentation_volume.shape(1);
    std::size_t z_size = np_segmentation_volume.shape(2);

    std::cout << "Begin Converting input data to boost multi_array"  << std::endl;
    volume_ptr<uint32_t> b_segmentation_ptr(new b::multi_array_ref<uint32_t, 3> (
                reinterpret_cast<uint32_t *>(np_segmentation_volume.get_data()),
                               boost::extents[x_size][y_size][z_size]
                               ));
    std::cout << "Finished converting segmentation input volume to boost multi_array "  <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting list of boost tuples to vector of tuples " 
        << bp::len(bp_region_graph_edges) << std::endl;
    region_graph_ptr<uint32_t, float> temp_region_graph_ptr(new region_graph<uint32_t, float>);

    for (int i=0; i < bp::len(bp_region_graph_edges); i++)
    {
        bp::tuple region_graph_edge = bp::extract<bp::tuple>(bp_region_graph_edges[i]);
        temp_region_graph_ptr->emplace_back(
                bp::extract<float>(region_graph_edge[0]),
                bp::extract<uint32_t>(region_graph_edge[1]),
                bp::extract<uint32_t>(region_graph_edge[2]));
    }
    std::cout << "Finished converting region_graph_edges to boost multi_array in "  <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::vector<size_t> b_counts;
    for (int i=0; i < bp::len(bp_counts); i++)
    {
        b_counts.emplace_back(bp::extract<uint32_t>(bp_counts[i]));
    }
    std::cout << "Finished converting bp_counts to vector in "  <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    // CALL MERGE_SEGMENTS_WITH_FUNCTION
    std::cout << "Begin calling merge segments with function"  << std::endl;
    merge_segments_with_function(b_segmentation_ptr, temp_region_graph_ptr,
            b_counts, const_above_threshold(min_mult, thold), lowt);
    std::cout << "Finished calling merge segments with function in "  <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting region_graph back to ndarray" << std::endl;
    while (bp::len(bp_region_graph_edges) > 0)
    {
        bp_region_graph_edges.pop();
    }
    std::cout << "Finished clearing original region graph list in " <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();
    for (auto e : *temp_region_graph_ptr )
    {
        float weight = std::get<0>(e);
        uint32_t node1 = std::get<1>(e);
        uint32_t node2 = std::get<2>(e);
        bp::tuple edge = bp::make_tuple(weight, node1, node2);
        bp_region_graph_edges.append(edge);
    }
    std::cout << "Finished converting output region_graph_edges in " <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    while (bp::len(bp_counts) > 0)
    {
        bp_counts.pop();
    }
    std::cout << "Finished clearing original counts graph list in " <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    for (size_t i = 0; i < b_counts.size(); ++i )
    {
        bp_counts.append(b_counts[i]);
    }
}


/**
 * Merges a tree to maximal or minimal spanning tree depending on order of inputs
 * @param region_graph_edges the region graph as a list of edge weights and nodes (Weight, Edge1Id, Edge2Id)
 * @param max_segid the maximum segmentid in the segmentation volume
 */
bp::list py_merge_tree(bp::list region_graph_edges, std::size_t max_segid) {
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin Converting input data to boost multi_array"  << std::endl;
    region_graph_ptr<uint32_t, float> temp_region_graph_ptr(new region_graph<uint32_t, float>);
    for (int i=0; i < bp::len(region_graph_edges); i++)
    {
        bp::tuple region_graph_edge = bp::extract<bp::tuple>(region_graph_edges[i]);
        temp_region_graph_ptr->emplace_back(
                bp::extract<float>(region_graph_edge[0]),
                bp::extract<uint32_t>(region_graph_edge[1]),
                bp::extract<uint32_t>(region_graph_edge[2]));
    }
    std::cout << "Finished converting region_graph boost multi_array "  <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    // CALL MERGE_TREE
    std::cout << "Begin call merge_tree to MST" << std::endl;
    region_graph_ptr<uint32_t, float> merged_tree = get_merge_tree(*temp_region_graph_ptr, max_segid);
    std::cout << "Finished calling merge_tree " << 
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting merge_tree back to ndarray" << std::endl;
    bp::list merge_tree_edges;
    for (auto e : *temp_region_graph_ptr )
    {
        float weight = std::get<0>(e);
        uint32_t node1 = std::get<1>(e);
        uint32_t node2 = std::get<2>(e);
        bp::tuple edge = bp::make_tuple(weight, node1, node2);
        merge_tree_edges.append(edge);
    }
    std::cout << "Finished converting output merge_tree in " <<
       chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    return merge_tree_edges;
}

bp::tuple py_watershed_full(np::ndarray& np_affinity_graph, bp::dict config)
{
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    std::size_t channel_size = np_affinity_graph.shape(0);
    std::size_t z_size = np_affinity_graph.shape(1);
    std::size_t y_size = np_affinity_graph.shape(2);
    std::size_t x_size = np_affinity_graph.shape(3);
    std::cout << "Begin Converting input data to boost multi_array"  << std::endl;

    affinity_graph_ptr<float> b_affinityGraph_ptr(new b::multi_array_ref<float, 4> (
                reinterpret_cast<float *>(np_affinity_graph.get_data()),
                               boost::extents[x_size][y_size][z_size][channel_size],
                               boost::fortran_storage_order()));
    std::cout << "Finished converting affinity_graph to boost multi_array "  <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    // Merge the input configuration into the default configuration
    ws::Config watershed_config = ws::Config();
    watershed_config.xSize = x_size;
    watershed_config.ySize = y_size;
    watershed_config.zSize = z_size;

    // leverage the command line parser to set up watershed config defaults
    std::vector<std::string> keys = extract_dict_keys(config);
    std::vector<std::string> argv_builder;
    for (size_t k=0; k < keys.size(); k++ )
    {
        argv_builder.push_back("--" + keys[k]);
        argv_builder.push_back(bp::extract<std::string>(bp::str(config.get(keys[k]))));
    }
    char *fake_argv[argv_builder.size()];
    for (size_t i = 0; i < argv_builder.size(); ++i)
    {
        fake_argv[i] = const_cast<char *>(argv_builder[i].c_str());
    }

    ws::config::parseCmdLine(argv_builder.size(), fake_argv, watershed_config, false, false);

    volume_ptr<uint32_t> b_segmentation_volume;
    region_graph_ptr<uint32_t, float> mst;

    // CALL Watershed full
    std::cout << "Begin call merge_tree to MST" << std::endl;
    std::tie(b_segmentation_volume, mst) =
        watershed_full<uint32_t, float>(b_affinityGraph_ptr, watershed_config);
    std::cout << "Finished calling watershed full" << 
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    std::cout << "Begin converting outputs back to python" << std::endl;
    np::ndarray np_segmentation_volume = np::from_data(
            reinterpret_cast<void *> (b_segmentation_volume->data()),
            np::dtype::get_builtin<uint32_t>(),
            bp::make_tuple(x_size, y_size, z_size),
            bp::make_tuple(y_size * z_size * sizeof(uint32_t), 
                z_size * sizeof(uint32_t),
                sizeof(uint32_t)),
            bp::object());

    std::cout << "Finished converting segmentation volume back to ndarray in " <<
        chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - t1).count() 
         << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    bp::list merge_tree_edges;
    for (auto e : *mst )
    {
        float weight = std::get<0>(e);
        uint32_t node1 = std::get<1>(e);
        uint32_t node2 = std::get<2>(e);
        bp::tuple edge = bp::make_tuple(weight, node1, node2);
        merge_tree_edges.append(edge);
    }
    std::cout << "Finished converting MST" <<
       chrono::duration_cast<chrono::milliseconds>(
               chrono::high_resolution_clock::now() - t1).count() 
       << " ms" << std::endl;
    t1 = chrono::high_resolution_clock::now();

    return bp::make_tuple(np_segmentation_volume, merge_tree_edges);
}


BOOST_PYTHON_MODULE(PyWatershed) {
    Py_Initialize();
    boost::numpy::initialize();
    bp::def("watershed", py_watershed);
    bp::def("regionGraph", py_region_graph);
    bp::def("mergeSegments", py_merge_segments);
    bp::def("mergeTree", py_merge_tree);
    bp::def("watershedFull", py_watershed_full);
}
/*
 *
 *int main(int argc, char *argv[])
 *{
 *
 *    Py_Initialize();
 *    boost::numpy::initialize();
 *
 *    try {
 *        std::cout << "BEGIN" << std::endl;
 *        np::ndarray data = np::zeros(bp::make_tuple(72), np::dtype::get_builtin<float>());
 *
 *        for (uint8_t i = 0; i < 72; ++i) {
 *            data[i] = i+1;
 *        }
 *        data = data.reshape(bp::make_tuple(3, 4, 2, 3));
 *
 *        std::cout << "Reshaped Array:\n" << bp::extract<char const*> (bp::str(data)) << std::endl;
 *
 *        // run watershed
 *        bp::tuple tuple = py_watershed(data, .3, .9);
 *
 *        std::cout << "Show results:" << std::endl;
 *        np::ndarray segVolume = bp::extract<np::ndarray>(tuple[0]);
 *        bp::list segCounts = bp::extract<bp::list>(tuple[1]);
 *
 *        std::cout << "returned segmentation :: " << std::endl << bp::extract < char const * > (bp::str(segVolume)) << std::endl;
 *        std::cout << "returned counts :: " << std::endl << bp::extract < char const * > (bp::str(segCounts)) << std::endl;
 *
 *        // get region graph
 *        bp::list region_graph = py_region_graph(data, segVolume, bp::len(segCounts)-1);
 *
 *        std::cout << "returned region counts :: " << std::endl << bp::extract < char const * > (bp::str(region_graph)) << std::endl;
 *
 *        // merge segments
 *        py_merge_segments(segVolume, region_graph, segCounts, 0, .3, 256);
 *
 *        std::cout << "returned region graph after merge :: " << std::endl << bp::extract < char const * > (bp::str(region_graph)) << std::endl;
 *
 *        // get MST
 *        bp::list merge_tree = py_merge_tree(region_graph, segCounts.shape(0) - 1);
 *
 *        std::cout << "returned MST after merge :: " << std::endl << bp::extract < char const * > (bp::str(merge_tree)) << std::endl;
 *    }
 *    catch (bp::error_already_set)
 *    {
 *        PyObject *ptype, *pvalue, *ptraceback;
 *        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
 *
 *        bp::handle<> hType(ptype);
 *        bp::object extype(hType);
 *        bp::handle<> hTraceback(ptraceback);
 *        bp::object traceback(hTraceback);
 *
 *        //Extract error message
 *        std::string strErrorMessage = bp::extract<std::string>(pvalue);
 *
 *        //Extract line number (top entry of call stack)
 *        // if you want to extract another levels of call stack
 *        // also process traceback.attr("tb_next") recurently
 *        long lineno = bp::extract<long> (traceback.attr("tb_lineno"));
 *        std::string filename = bp::extract<std::string>(traceback.attr("tb_frame").attr("f_code").attr("co_filename"));
 *        std::string funcname = bp::extract<std::string>(traceback.attr("tb_frame").attr("f_code").attr("co_name"));
 *        std::cout << strErrorMessage << std::endl;
 *        std::cout << filename << std::endl;
 *        std::cout << funcname << std::endl;
 *    }
 *    return 0;
 *}
 *
 */
