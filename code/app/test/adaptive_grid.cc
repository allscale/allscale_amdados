#ifndef AMDADOS_PLAIN_MPI
// This implementation is based on Allscale API, no MPI at all.

#include <gtest/gtest.h>
#include "allscale/utils/assert.h"
#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/utils/printer/vectors.h"
#include "allscale/utils/printer/pairs.h"
#include "allscale/utils/static_grid.h"
#include "allscale/utils/string_utils.h"
#include <ctime>
#include <cstdlib>
#include <iostream>

namespace allscale {
namespace api {
namespace user {
namespace data {

// Stuff from the file:  "data_item_test.inl"
namespace {

template<typename Region>
void checkParameters(const Region& a, const Region& b) {

	using allscale::api::core::is_region;
	using allscale::api::core::isSubRegion;

	// very first check: syntax
	static_assert(is_region<Region>::value, "Provided region does not fit region concept!");

	// check input parameters
	ASSERT_FALSE(a.empty()) << "Require first parameter to be non-empty: " << a;
	ASSERT_FALSE(b.empty()) << "Require second parameter to be non-empty: " << b;
	ASSERT_NE(a,b) << "Requires parameters to be not equivalent.";

	ASSERT_FALSE(isSubRegion(a,b));
	ASSERT_FALSE(isSubRegion(b,a));

	// compute intersection
	Region c = Region::intersect(a,b);
	ASSERT_FALSE(c.empty());
}


template<typename Region>
void testRegion(const Region& a, const Region& b) {

	using allscale::api::core::isSubRegion;

	// check whether parameters are two distinct set with a non-empty common sub-set
	checkParameters(a,b);

	// compute intersection
	Region c = Region::intersect(a,b);
	Region d = Region::merge(a,b);

	// so, we have:
	//		region a and b -- neither a subset of the other
	//		region c - a non-empty sub-set of region a and b
	//      region d - the union of a and b


	// test semantic
	Region e;
	EXPECT_TRUE(e.empty());


	// check sub-region relation

	EXPECT_TRUE(isSubRegion(e,e));
	EXPECT_TRUE(isSubRegion(a,a));
	EXPECT_TRUE(isSubRegion(b,b));
	EXPECT_TRUE(isSubRegion(c,c));
	EXPECT_TRUE(isSubRegion(d,d));

	EXPECT_TRUE(isSubRegion(e,a));
	EXPECT_TRUE(isSubRegion(e,b));
	EXPECT_TRUE(isSubRegion(e,c));
	EXPECT_TRUE(isSubRegion(e,d));

	EXPECT_FALSE(isSubRegion(a,e));
	EXPECT_FALSE(isSubRegion(b,e));
	EXPECT_FALSE(isSubRegion(c,e));
	EXPECT_FALSE(isSubRegion(d,e));

	EXPECT_FALSE(isSubRegion(a,c));
	EXPECT_FALSE(isSubRegion(b,c));
	EXPECT_FALSE(isSubRegion(d,c));

	EXPECT_FALSE(isSubRegion(d,a));
	EXPECT_FALSE(isSubRegion(d,b));

	EXPECT_TRUE(isSubRegion(c,a));
	EXPECT_TRUE(isSubRegion(c,b));
	EXPECT_TRUE(isSubRegion(a,d));
	EXPECT_TRUE(isSubRegion(b,d));


	// check equivalencis
	EXPECT_EQ(e, Region::difference(e,e));
	EXPECT_EQ(e, Region::difference(a,a));
	EXPECT_EQ(e, Region::difference(a,d));
	EXPECT_EQ(e, Region::difference(b,b));
	EXPECT_EQ(e, Region::difference(b,d));

	EXPECT_EQ(Region::difference(b,a), Region::difference(b,c)) << "a=" << a << "\nb=" << b << "\nc=" << c;
	EXPECT_EQ(Region::difference(a,b), Region::difference(a,c)) << "a=" << a << "\nb=" << b << "\nc=" << c;

}


template<typename Fragment>
void testFragment(const typename Fragment::shared_data_type& shared, const typename Fragment::region_type& a, const typename Fragment::region_type& b) {

	using Region = typename Fragment::region_type;

	// check whether parameters are two distinct set with a non-empty common sub-set
	checkParameters(a,b);

	// compute empty region, intersection, and union of a and b
	Region e;
	Region c = Region::intersect(a,b);
	Region d = Region::merge(a,b);


	// create a few fragments
	Fragment empty(shared,e);
	EXPECT_EQ(e, empty.getCoveredRegion());

	Fragment fA(shared,a);
	EXPECT_EQ(a, fA.getCoveredRegion());

	Fragment fB(shared,b);
	EXPECT_EQ(b, fB.getCoveredRegion());

	// create an empty fragment
	Fragment tmp(shared,e);
	EXPECT_EQ(e, tmp.getCoveredRegion());

	// resize fragment to size of c
	tmp.resize(c);
	EXPECT_EQ(c, tmp.getCoveredRegion());

	// load c-share for fA into tmp
	tmp.insert(fA,c);
	EXPECT_EQ(c, tmp.getCoveredRegion());

	// resize tmp to d
	tmp.resize(d);
	EXPECT_EQ(d, tmp.getCoveredRegion());

	// load fA and fB into tmp
	tmp.insert(fA, Region::difference(a,c));
	tmp.insert(fB, Region::difference(b,a));

	auto facada = fA.mask();
}


template<typename Fragment>
void testFragment(const typename Fragment::region_type& a, const typename Fragment::region_type& b) {
	testFragment<Fragment>(typename Fragment::shared_data_type(), a, b);
}


template<typename Fragment, typename Region>
utils::Archive extract(const Fragment& fragment, const Region& region) {
	utils::ArchiveWriter writer;
	fragment.extract(writer,region);
	return std::move(writer).toArchive();
}

template<typename Fragment>
void insert(Fragment& fragment, const utils::Archive& archive) {
	utils::ArchiveReader reader(archive);
	fragment.insert(reader);
}

} // anonymous namespace



    TEST(MyAdaptiveGridCell, MyRefinementCoarsening) {
        const int NUM_LAYERS = 3;
        using FourLayerCellConfig = CellConfig<2, layers<layer<1,1>,
                                                         layer<8,8>,
                                                         layer<2,2> > >;

        AdaptiveGridCell<int, FourLayerCellConfig> cell;
//        FourLayerCellConfig config;
        std::cout << FourLayerCellConfig().dims << std::endl;
        std::cout << FourLayerCellConfig().num_layers << std::endl;
        assert_true((NUM_LAYERS == FourLayerCellConfig().num_layers));
//        std::exit(0);

        // Function prints all the layers.
        auto PrintAllLayers = [&](
                AdaptiveGridCell<int, FourLayerCellConfig> & cell)
        {
            unsigned int old = cell.getActiveLayer();
            for (unsigned int i = 0; i < NUM_LAYERS; ++i) {
                cell.setActiveLayer(i);
                const auto sz = cell.getActiveLayerSize();
                std::cout << "layer: " << i << ", size: " << sz << std::endl;
                for (unsigned int y = 0; y < sz.y; ++y) {
                    for (unsigned int x = 0; x < sz.x; ++x) {
                        std::cout << std::setw(4) << cell[{x,y}];
                    }
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
            std::cout << "---------------------------------------------------";
            std::cout << std::endl << std::endl << std::flush;
            cell.setActiveLayer(old);
        };

        // Initialize. Important!
        for (unsigned int i = 0; i < NUM_LAYERS; ++i) {
            cell.setActiveLayer(i);
            //cell.forAllActiveNodes([=](int& element) { element = i; });
            const auto sz = cell.getActiveLayerSize();
            std::cout << "layer: " << i << ", size: " << sz << std::endl;
            unsigned int count = 0;
            for (unsigned int y = 0; y < sz.y; ++y) {
                for (unsigned int x = 0; x < sz.x; ++x) {
                    assert_true(count == x + y*sz.x);
                    cell[{x,y}] = ++count;
                }
            }
            assert_true(count == sz.x * sz.y);
        }
        PrintAllLayers(cell);

        // Show an alternative way to access the layer.
        auto & cellGrid = cell.getLayer<0>();
        cellGrid.forEach([&](auto& element) { element = 77; });
        PrintAllLayers(cell);

        // Test different ways to access layer's elements.
        for (unsigned int i = 0; i < NUM_LAYERS; ++i) {
            cell.setActiveLayer(0);
            unsigned long count = 0;
            const auto sz = cell.getActiveLayerSize();
            cell.forAllActiveNodes(
                [&](const allscale::utils::Vector<long,2> & pos, int & cur) {
                    const long x = pos.x;
                    const long y = pos.y;
                    assert_true(count == x*sz.y + y);   // y is faster!
                    cell[{x,y}] = ++count;
                    assert_true((cell[{x,y}] == cur));
            });
            assert_true(count == sz.x * sz.y);
        }
        std::cout << "#########################################" << std::endl;
        PrintAllLayers(cell);
        std::cout << "#########################################" << std::endl;

        for (int i = 0; i < NUM_LAYERS - 1; ++i) {
            cell.setActiveLayer(i);
            cell.coarsen([](const int & element) { return element; });
            PrintAllLayers(cell);
        }
        for (int i = NUM_LAYERS - 1; i > 0; --i) {
            cell.setActiveLayer(i);
            cell.refine([](const int & element) { return element; });
            PrintAllLayers(cell);
        }

//        cell.setActiveLayer(1);
//        cell.coarsen([](const int & element) { return element; });
//        PrintAllLayers(cell);
//
//        cell.setActiveLayer(2);
//        cell.coarsen([](const int & element) { return element; });
//        PrintAllLayers(cell);
//
//        cell.setActiveLayer(3);
//        cell.refine([](const int & element) { return element; });
//        PrintAllLayers(cell);
//
//        cell.setActiveLayer(2);
//        cell.refine([](const int & element) { return element; });
//        PrintAllLayers(cell);
//
//        cell.setActiveLayer(1);
//        cell.refine([](const int & element) { return element; });
//        PrintAllLayers(cell);

        std::srand(std::time(nullptr));
        for (int test = 0; test < 10; ++test) {
            for (int l = 0; l < 4; ++l) {
                cell.setActiveLayer(l);
                const auto sz = cell.getActiveLayerSize();
                const int W = static_cast<int>(sz.x);
                const int H = static_cast<int>(sz.y);
                if (!((W > 1) && (H > 1)))
                    continue;

                bool ok = false;
                std::vector<int> boundary;
                std::vector<int> left( static_cast<size_t>(H)),
                                 right(static_cast<size_t>(H));
                std::vector<int> bottom(static_cast<size_t>(W)),
                                 top(   static_cast<size_t>(W));

                for (int y = 0; y < H; ++y) {
                    left [static_cast<size_t>(y)] = cell[{0  ,y}] = std::rand();
                    right[static_cast<size_t>(y)] = cell[{W-1,y}] = std::rand();
                }
                boundary = cell.getBoundary(Left);
                ok = std::equal(left.begin(), left.end(), boundary.begin());
                assert_true(ok);
                boundary = cell.getBoundary(Right);
                ok = std::equal(right.begin(), right.end(), boundary.begin());
                assert_true(ok);

                for (int x = 0; x < W; ++x) {
                    bottom[static_cast<size_t>(x)] = cell[{x,  0}] = std::rand();
                    top   [static_cast<size_t>(x)] = cell[{x,H-1}] = std::rand();
                }
                boundary = cell.getBoundary(Down);
                ok = std::equal(bottom.begin(), bottom.end(), boundary.begin());
                assert_true(ok);
                boundary = cell.getBoundary(Up);
                ok = std::equal(top.begin(), top.end(), boundary.begin());
                assert_true(ok);

                std::generate(  left.begin(),   left.end(), std::rand);
                std::generate( right.begin(),  right.end(), std::rand);
                std::generate(bottom.begin(), bottom.end(), std::rand);
                std::generate(   top.begin(),    top.end(), std::rand);
                cell.setBoundary(Left, left);
                cell.setBoundary(Right, right);
                for (int y = 0; y < H; ++y) {
                    ok = (left [static_cast<size_t>(y)] == cell[{0  ,y}]);
                    assert_true(ok);
                    ok = (right[static_cast<size_t>(y)] == cell[{W-1,y}]);
                    assert_true(ok);
                }
                for (int x = 0; x < W; ++x) {
                    ok = (bottom[static_cast<size_t>(x)] = cell[{x,  0}]);
                    assert_true(ok);
                    ok = (top   [static_cast<size_t>(x)] = cell[{x,H-1}]);
                    assert_true(ok);
                }
            }
        }
    }

} // end namespace data
} // end namespace user
} // end namespace api
} // end namespace allscale

#endif  // AMDADOS_PLAIN_MPI
