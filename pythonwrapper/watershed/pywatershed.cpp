#define NPY_DEPRECATED_API NPY_1_7_API_VERSION
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numpy.hpp>
#include <boost/multi_array.hpp>
#include <boost/multi_array/types.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/return_value_policy.hpp>


#include <string>
#include <memory>
#include <cstdint>
#include <assert.h>
#include <tuple>

#include "basic_watershed.hpp"
#include "region_graph.hpp"
#include "agglomeration.hpp"
#include "merge_tree.hpp"

namespace b = boost;
namespace bp = boost::python;
namespace np = boost::numpy;

template <typename TAffinity>
bp::tuple pyWatershed(np::ndarray& npAffinityGraph, TAffinity lowv, TAffinity highv)
{
	std::cout<< lowv << " HIGH IS " << highv << std::endl;
	std::cout << "Original array :: " << std::endl <<
		bp::extract<char const *>(bp::str(npAffinityGraph)) << std::endl;
	 std::size_t channelSize = npAffinityGraph.shape(0);
	 std::size_t zSize = npAffinityGraph.shape(1);
	 std::size_t ySize = npAffinityGraph.shape(2);
	 std::size_t xSize = npAffinityGraph.shape(3);
	std::cout << "Original dimensions :: " << std::endl <<
		channelSize << "," << zSize << "," << ySize << "," << xSize <<
		" dtype= " << bp::extract<char const *> (bp::str(npAffinityGraph.get_dtype())) <<
		std::endl;

	for (int c = 0; c < channelSize; ++c) {
		for (int z = 0; z < zSize; ++z) {
			for (int y = 0; y < ySize; ++y) {
				for (int x = 0; x < xSize; ++x) {
	/*
	 *for (int x = 0; x < xSize; ++x) {
	 *    for (int y = 0; y < ySize; ++y) {
	 *        for (int z = 0; z < zSize; ++z) {
	 *            for (int c = 0; c < channelSize; ++c) {
	 */
					std::cout << bp::extract<char const *>(bp::str(npAffinityGraph[c][z][y][x])) << std::endl;
				}
			}
		}
	}
	std::cout << "np ^^ done" << std::endl;

    affinity_graph_ptr<float> bAffinityGraph_ptr(new b::multi_array_ref<float, 4> (
				reinterpret_cast<float *>(npAffinityGraph.get_data()),
							   boost::extents[xSize][ySize][zSize][channelSize],
							   boost::fortran_storage_order()));
	for (int c = 0; c < channelSize; ++c) {
		for (int z = 0; z < zSize; ++z) {
			for (int y = 0; y < ySize; ++y) {
				for (int x = 0; x < xSize; ++x) {
					std::cout << x << "," << y << "," << z << "," << c << " = "
						<< unsigned((*bAffinityGraph_ptr)[x][y][z][c]) << std::endl;
				}
			}
		}
	}


	
	std::tuple<volume_ptr<uint32_t>, std::vector<std::size_t>> tuple = watershed<uint32_t>(bAffinityGraph_ptr, lowv, highv);

	b::multi_array_ref<uint32_t, 3> bSegmentationVolume = *(std::get<0>(tuple));
	np::ndarray npSegmentationVolume = np::from_data(
			reinterpret_cast<void *> (bSegmentationVolume.data()),
			np::dtype::get_builtin<uint32_t>(),
			bp::make_tuple(xSize, ySize, zSize, channelSize),
			bp::make_tuple(ySize * zSize * channelSize * sizeof(uint32_t), 
				zSize * channelSize * sizeof(uint32_t),
				channelSize * sizeof(uint32_t),
				sizeof(uint32_t)),
			bp::object());

	std::vector<size_t> segmentationCounts = std::get<1>(tuple);
	np::ndarray npSegmentationCounts = np::from_data(
			reinterpret_cast<void *> (segmentationCounts.data()),
			np::dtype::get_builtin<size_t>(),
			bp::make_tuple(segmentationCounts.size()),
			bp::make_tuple(sizeof(size_t)),
			bp::object());

	int data[] = { 1 , 2, 3 ,4,5};
	np::ndarray dummyarr  = np::from_data(
			data, np::dtype::get_builtin<int>(),
			bp::make_tuple(5), bp::make_tuple(sizeof(int)), bp::object());

	return bp::make_tuple(npSegmentationVolume, npSegmentationCounts, dummyarr);
};

	bp::object own;
np::ndarray blah() {
	int data[] = {1,2,3,4,5} ;
	bp::tuple shape = bp::make_tuple(5) ;
	bp::tuple stride = bp::make_tuple(sizeof(int)) ;
	np::dtype dt = np::dtype::get_builtin<int>();
	np::ndarray data_ex = np::from_data(data,dt, shape,stride,own);
	std::cout << "Single dimensional array ::" << std::endl << bp::extract < char const * > (bp::str(data_ex)) << std::endl;

	return data_ex;
}

//void (*pws)(np::ndarray, float, float) = pyWatershed;

BOOST_PYTHON_MODULE(PyWatershed) {
	Py_Initialize();
	boost::numpy::initialize();
	bp::def("watershed", pyWatershed<float>);
	bp::def("watershed", pyWatershed<double>);
	bp::def("blah", blah);//, bp::return_value_policy<bp::manage_new_object>());
}
/*
 *int main(int argc, char *argv[])
 *{
 *
 *    Py_Initialize();
 *    boost::numpy::initialize();
 *    std::cout << "BEGIN" << std::endl;
 *    np::ndarray data = np::zeros(bp::make_tuple(72), np::dtype::get_builtin<float>());
 *
 *    for (uint8_t i = 0; i < 72; ++i) {
 *        data[i] = i+1;
 *    }
 *    data = data.reshape(bp::make_tuple(3, 4, 2, 3));
 *
 *    std::cout << "Reshaped Array:\n" << bp::extract<char const*> (bp::str(data)) << std::endl;
 *
 *    pyWatershed(data, .3, .4);
 *    return 0;
 *};
 */
