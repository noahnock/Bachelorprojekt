//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include <iostream>
#include <vector>
#include "./Rtree/Rtree.h"
#include "Rtree/RtreeFileReader.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pointGeo = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
using polygonGeo = bg::model::polygon<pointGeo>;
using boxGeo = bg::model::box<pointGeo>;
using rTreeValue = std::pair<boxGeo, long long>;
using multiBox = std::vector<rTreeValue>;

class BoostRtree {
private:
    bgi::rtree<rTreeValue , bgi::rstar<16>> rtree_;

public:
    void createEntries(const multiBox& entries) {
        std::cout << "Creating the entries in the R-tree..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();

        rtree_ = bgi::rtree<rTreeValue, bgi::rstar<16>>(entries.begin(), entries.end());

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << "Created the R-tree with " << entries.size() << " entries in " << duration.count() / 1000000.0 << " seconds" << std::endl;
    }

    multiBox SearchInTree(boxGeo query) {
        std::cout << "Searching in the R-tree..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        multiBox results;
        std::for_each(rtree_.qbegin(bgi::intersects(query)), rtree_.qend(), [&](rTreeValue const& result) {
            results.push_back(result);
        });

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << results.size() << " elements were found in " << duration.count() / 1000000.0 << " seconds" << std::endl;
        return results;
    }
};

void searchInBoostRtree() {
    BoostRtree tree;
    multiBoxGeo rectangles = FileReaderWithoutIndex::LoadEntries("../switzerland_raw/converted_data_100k.boundingbox.tmp");
    multiBox entries;
    for (RTreeValue entry : rectangles) {
        rTreeValue newEntry = std::make_pair(entry.box, entry.id);
        entries.push_back(newEntry);
    }
    tree.createEntries(entries);
    tree.SearchInTree(BasicGeometry::CreateBoundingBox(9.88657, 47.38431, 9.88671, 47.6088));
}

void Benchmark() {
    searchInBoostRtree();
}
