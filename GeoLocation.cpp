// Copyright 2023, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Noah Nock (noah.v.nock@gmail.com)

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include "Rtree/Rtree.h"
#include <filesystem>
#include <regex>
#include "Rtree/ExternalSorting.cpp"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pointGeo = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
using polygonGeo = bg::model::polygon<pointGeo>;
using multiPolygonGeo = bg::model::multi_polygon<polygonGeo>;
using boxGeo = bg::model::box<pointGeo>;
using rTreeValue = std::pair<boxGeo, long long>;
using multiBoxGeo = std::vector<rTreeValue>;
using boundingBoxValues = std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>;
using bg::append;
using bg::make;

class GeoLocation {
private:
    bgi::rtree<rTreeValue , bgi::rstar<16>> rtree;

public:

    static multiBoxGeo loadEntries(const std::string& path) {
        std::ifstream infile(path);

        multiBoxGeo boxes;
        unsigned int lineCount = 0;

        std::cout << "Loading in the entries..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();

        for (std::string line; std::getline(infile, line);)
        {

            if (line != "minX minY,maxX maxY,id") {
                lineCount++;

                std::string mins = line.substr(0, line.find(','));
                //std::string maxs = line.substr(line.find(',') + 1, line.length() - (mins.length() + 1));
                std::string maxs = line.substr(line.find(',') + 1, line.find(',', line.find(',') + 1)- line.find(',') - 1);
                unsigned int id = std::stod(line.substr(line.find(',', line.find(',') + 1) + 1, line.length() - (line.find(',', line.find(',') + 1) + 1)));

                double minX = std::stod(mins.substr(0, mins.find(' ')));
                double minY = std::stod(mins.substr(mins.find(' ') + 1, mins.length() - (mins.find(' ') + 1)));
                double maxX = std::stod(maxs.substr(0, maxs.find(' ')));
                double maxY = std::stod(maxs.substr(maxs.find(' ') + 1, maxs.length() - (maxs.find(' ') + 1)));

                boxGeo box = createBoundingBox(minX, minY, maxX, maxY);
                rTreeValue boxWithId = std::make_pair(box, id);
                boxes.push_back(boxWithId);

                if (lineCount % 5000000 == 0) {
                    std::cout << "Reached entry number " << lineCount << std::endl;
                }
            }

        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << "Loaded " << lineCount << " entries in " << duration.count() / 1000000.0 << " seconds" << std::endl;

        return boxes;
    }

    void createEntries(const multiBoxGeo& entries) {
        std::cout << "Creating the entries in the R-Tree..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();

        rtree = bgi::rtree<rTreeValue, bgi::rstar<16>>(entries.begin(), entries.end());

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << "Created the R-Tree with " << entries.size() << " entries in " << duration.count() / 1000000.0 << " seconds" << std::endl;
    }

    void addEntry(rTreeValue value) {
        rtree.insert(value);
    }

    unsigned int removeEntry(rTreeValue value) {
        return rtree.remove(value);
    }

    multiBoxGeo SearchInTree(boxGeo query) {
        multiBoxGeo results;
        std::for_each(rtree.qbegin(bgi::intersects(query)), rtree.qend(), [&](rTreeValue const& result) {
            results.push_back(result);
            //std::cout << bg::wkt(query) << " intersects with " << bg::wkt(result.first) << ", which has the Id: " << result.second << std::endl;
        });

        return results;
    }

    static boxGeo createBoundingBox(double pointOneX, double pointOneY, double pointTwoX, double pointTwoY) {
        return make<boxGeo>(make<pointGeo>(pointOneX, pointOneY), make<pointGeo>(pointTwoX, pointTwoY));
    }

    static boundingBoxValues convertWKTToBoundingBoxes(const std::string& wkt) {

        std::size_t posMPolStart = wkt.find("MULTIPOLYGON(((");
        std::size_t posMPolEnd = wkt.find(")))") + 2;

        if (posMPolStart == std::string::npos || posMPolEnd == std::string::npos) {
            return boundingBoxValues {};
        }

        std::string newWkt = wkt.substr(posMPolStart, posMPolEnd - posMPolStart + 1);

        multiPolygonGeo mPolygon;
        bg::read_wkt(newWkt, mPolygon);

        boundingBoxValues boxes(mPolygon.size());

        for (size_t i = 0; i < mPolygon.size(); i++) {
            polygonGeo polygon = mPolygon[i];

            double minX = -1;
            double maxX = -1;
            double minY = -1;
            double maxY = -1;

            for (auto& p : polygon.outer()) {
                double x = bg::get<0>(p);
                double y = bg::get<1>(p);

                if (x < minX || minX == -1) {
                    minX = x;
                }

                if (x > maxX) {
                    maxX = x;
                }

                if (y < minY || minY == -1) {
                    minY = y;
                }

                if (y > maxY) {
                    maxY = y;
                }
            }

            boxes[i] = std::make_pair(std::make_pair(minX, minY), std::make_pair(maxX, maxY));
        }

        return boxes;
    }

    static void convertWKTFile(const std::string& fileFromPath, const std::string& fileToPath) {
        std::ifstream infile(fileFromPath);
        size_t lineCount = 0;
        size_t id = 1;

        std::cout << "Starting to convert the file" << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        std::ofstream outfile;
        outfile.open(fileToPath, std::ios_base::app);
        outfile << "minX minY,maxX maxY,id" << std::endl;
        for (std::string line; std::getline(infile, line); )
        {
            lineCount++;
            boundingBoxValues boxes = convertWKTToBoundingBoxes(line);
            if (!boxes.empty()) {
                for (const auto & box : boxes) {
                    double minX = box.first.first;
                    double minY = box.first.second;
                    double maxX = box.second.first;
                    double maxY = box.second.second;
                    outfile << minX << " " << minY << "," << maxX << " " << maxY << "," << id << std::endl;
                    id++;
                }
            }
            if (lineCount % 100000 == 0) {
                std::cout << "Reached line number " << lineCount << std::endl;
            }
        }
        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << "Converted " << lineCount << " lines in " << duration.count() / 1000000.0 << " seconds" << std::endl;
    }

    static void centerOrdering(multiBoxGeo& boxes, size_t dim) {
        if (dim == 0) {
            // order by centerX
            auto sortRuleLambda = [] (rTreeValue b1, rTreeValue b2) -> bool
            {
                double center1 = (b1.first.min_corner().get<0>() + b1.first.max_corner().get<0>()) / 2;
                double center2 = (b2.first.min_corner().get<0>() + b2.first.max_corner().get<0>()) / 2;
                return center1 < center2;
            };

            std::sort(boxes.begin(), boxes.end(), sortRuleLambda);
        } else if (dim == 1) {
            // order by centerY
            auto sortRuleLambda = [] (rTreeValue b1, rTreeValue b2) -> bool
            {
                double center1 = (b1.first.min_corner().get<1>() + b1.first.max_corner().get<1>()) / 2;
                double center2 = (b2.first.min_corner().get<1>() + b2.first.max_corner().get<1>()) / 2;
                return center1 < center2;
            };

            std::sort(boxes.begin(), boxes.end(), sortRuleLambda);
        }
    }

    static double costFunctionTGS(boxGeo& b0, boxGeo& b1) {
        /*
         * To accelerate the algorithm, make sure that both boxes are of similar size.
         */

        double areaB0 = (b0.max_corner().get<0>() - b0.min_corner().get<0>()) * (b0.max_corner().get<1>() - b0.min_corner().get<1>());
        double areaB1 = (b1.max_corner().get<0>() - b1.min_corner().get<0>()) * (b1.max_corner().get<1>() - b1.min_corner().get<1>());
        double difference = std::abs(areaB0 - areaB1);

        return difference;
    }

    static bool pointWithinBox(pointGeo point, boxGeo box) {
        return point.get<0>() >= box.min_corner().get<0>() && point.get<0>() <= box.max_corner().get<0>() && point.get<1>() >= box.min_corner().get<1>() && point.get<1>() <= box.max_corner().get<1>();
    }

    static std::vector<multiBoxGeo> TGSRecursive(const std::vector<multiBoxGeo>& orderedInputRectangles, size_t S, size_t M) {  // https://dl.acm.org/doi/pdf/10.1145/288692.288723
        /*
         * inputRectangles needs to be pre-sorted with centerOrdering for both d0 and d1
         */

        size_t n = orderedInputRectangles[0].size();

        if (n <= S) {
            // stop condition
            return std::vector<multiBoxGeo> { orderedInputRectangles[0] };
        }

        double bestCost = -1;
        size_t bestDim = 0;
        unsigned long bestI = 1;
        boxGeo bestB0 = createBoundingBox(0, 0, 0, 0);

        for (size_t dim = 0; dim <= 1; dim++) {
            for (size_t i = 1; i < n / M; i++) {
                // calculate B0 and B1
                double minXB0 = 0;
                double minYB0 = 0;
                double maxXB0 = 1;
                double maxYB0 = 1;

                double minXB1 = 0;
                double minYB1 = 0;
                double maxXB1 = 1;
                double maxYB1 = 1;

                if (dim == 0) {
                    minXB0 = (orderedInputRectangles[dim][0].first.min_corner().get<0>() + orderedInputRectangles[dim][0].first.max_corner().get<0>()) / 2;
                    maxXB0 = (orderedInputRectangles[dim][i * S].first.min_corner().get<0>() + orderedInputRectangles[dim][i * S].first.max_corner().get<0>()) / 2;

                    minXB1 = (orderedInputRectangles[dim][i * S + 1].first.min_corner().get<0>() + orderedInputRectangles[dim][i * S + 1].first.max_corner().get<0>()) / 2;
                    maxXB1 = (orderedInputRectangles[dim][orderedInputRectangles[dim].size() - 1].first.min_corner().get<0>() + orderedInputRectangles[dim][orderedInputRectangles[dim].size() - 1].first.max_corner().get<0>()) / 2;
                } else {
                    minYB0 = (orderedInputRectangles[dim][0].first.min_corner().get<1>() + orderedInputRectangles[dim][0].first.max_corner().get<1>()) / 2;
                    maxYB0 = (orderedInputRectangles[dim][i * S].first.min_corner().get<1>() + orderedInputRectangles[dim][i * S].first.max_corner().get<1>()) / 2;

                    minYB1 = (orderedInputRectangles[dim][i * S + 1].first.min_corner().get<1>() + orderedInputRectangles[dim][i * S + 1].first.max_corner().get<1>()) / 2;
                    maxYB1 = (orderedInputRectangles[dim][orderedInputRectangles[dim].size() - 1].first.min_corner().get<1>() + orderedInputRectangles[dim][orderedInputRectangles[dim].size() - 1].first.max_corner().get<1>()) / 2;
                }

                boxGeo b0 = createBoundingBox(minXB0, minYB0, maxXB0, maxYB0);
                boxGeo b1 = createBoundingBox(minXB1, minYB1, maxXB1, maxYB1);


                double cost = costFunctionTGS(b0, b1);

                if (bestCost == -1 || cost < bestCost) {
                    bestCost = cost;
                    bestDim = dim;
                    bestI = i;
                    bestB0 = b0;
                }
            }
        }

        // split the rectangles
        multiBoxGeo s0BestDim(orderedInputRectangles[bestDim].begin(), orderedInputRectangles[bestDim].begin() + bestI * S);
        multiBoxGeo s1BestDim(orderedInputRectangles[bestDim].begin() + bestI * S, orderedInputRectangles[bestDim].end());

        multiBoxGeo s0OtherDim;
        multiBoxGeo s1OtherDim;

        size_t otherDim = 1 - bestDim;

        for (rTreeValue box : orderedInputRectangles[otherDim]) {
            //if (elementBelongsToFirstBox(box, bestB0, bestB1, s0BestDim, bestDim)) {
            pointGeo boxCenter;
            if (bestDim == 0) {
                boxCenter = make<pointGeo>((box.first.min_corner().get<0>() + box.first.max_corner().get<0>()) / 2, 0.5);
            } else {
                boxCenter = make<pointGeo>(0.5, (box.first.min_corner().get<1>() + box.first.max_corner().get<1>()) / 2);
            }
            if (pointWithinBox(boxCenter, bestB0) && s0OtherDim.size() < bestI * S) {
                s0OtherDim.push_back(box);
            } else {
                s1OtherDim.push_back(box);
            }
        }

        std::vector<multiBoxGeo> s0;
        std::vector<multiBoxGeo> s1;

        if (bestDim == 0) {
            s0.push_back(s0BestDim);
            s0.push_back(s0OtherDim);
            s1.push_back(s1BestDim);
            s1.push_back(s1OtherDim);
        } else {
            s0.push_back(s0OtherDim);
            s0.push_back(s0BestDim);
            s1.push_back(s1OtherDim);
            s1.push_back(s1BestDim);
        }

        // recursion
        std::vector<multiBoxGeo> result0 = TGSRecursive(s0, S, M);
        std::vector<multiBoxGeo> result1 = TGSRecursive(s1, S, M);

        std::vector<multiBoxGeo> result;
        result.insert(result.begin(), result0.begin(), result0.end());
        result.insert(result.end(), result1.begin(), result1.end());

        return result;
    }

    static std::vector<multiBoxGeo> TGS(multiBoxGeo& inputRectangles, size_t S, size_t M) {
        std::vector<multiBoxGeo> orderedInputRectangles;

        multiBoxGeo RectanglesD0(inputRectangles);
        multiBoxGeo RectanglesD1(inputRectangles);

        centerOrdering(RectanglesD0, 0);
        centerOrdering(RectanglesD1, 1);

        orderedInputRectangles.push_back(RectanglesD0);
        orderedInputRectangles.push_back(RectanglesD1);

        return TGSRecursive(orderedInputRectangles, S, M);
    }

    static void saveTGSResult(const std::vector<multiBoxGeo>& result, const std::string& folder_path) {
        std::filesystem::create_directory(folder_path);

        std::ofstream baseFile;
        baseFile.open(folder_path + "/base.csv");
        baseFile << "minX minY,maxX maxY,id" << std::endl;

        size_t index = 1;
        for (const multiBoxGeo& boxes : result) {
            if (!boxes.empty()) {
                std::ofstream currentFile;
                currentFile.open(folder_path + "/file_" + std::to_string(index) + ".csv");
                currentFile << "minX minY,maxX maxY,id" << std::endl;

                double globalMinX = -1;
                double globalMinY = -1;
                double globalMaxX = -1;
                double globalMaxY = -1;

                for (rTreeValue box : boxes) {
                    double minX = box.first.min_corner().get<0>();
                    double minY = box.first.min_corner().get<1>();
                    double maxX = box.first.max_corner().get<0>();
                    double maxY = box.first.max_corner().get<1>();
                    currentFile << minX << " " << minY << "," << maxX << " " << maxY << "," << box.second << std::endl;

                    if (globalMinX == -1 || minX < globalMinX) {
                        globalMinX = minX;
                    }
                    if (globalMinY == -1 || minY < globalMinY) {
                        globalMinY = minY;
                    }
                    if (maxX > globalMaxX) {
                        globalMaxX = maxX;
                    }
                    if (maxY > globalMaxY) {
                        globalMaxY = maxY;
                    }
                }

                baseFile << globalMinX << " " << globalMinY << "," << globalMaxX << " " << globalMaxY << "," << index << std::endl;

                index++;
            }
        }
    }

    void Search(boxGeo query, const std::string& folder) {
        multiBoxGeo baseEntries = loadEntries(folder + "/base.csv");
        createEntries(baseEntries);

        multiBoxGeo baseResults = SearchInTree(query);
        for (rTreeValue baseResult : baseResults) {
            unsigned int id = baseResult.second;

            multiBoxGeo entries = loadEntries(folder + "/file_" + std::to_string(id) + ".csv");

            for (const rTreeValue& result : entries) {
                if (bg::intersects(result.first, query)) {
                    std::cout << result.first.min_corner().get<0>() << " " << result.first.min_corner().get<1>() << "," << result.first.max_corner().get<0>()
                              << " " << result.first.max_corner().get<1>() << "," << result.second << std::endl;
                }
            }

            //createEntries(entries);

            /*multiBoxGeo results = SearchInTree(query);
            for (rTreeValue result : results) {
                std::cout << result.first.min_corner().get<0>() << " " << result.first.min_corner().get<1>() << "," << result.first.max_corner().get<0>()
                        << " " << result.first.max_corner().get<1>() << "," << result.second << std::endl;
            }*/
        }

    }

    /*static void TestConverter(Rtree& tree, const std::string& path) {
        std::ifstream infile(path);
        tree.OpenConversion("../conversion_test_old");
        size_t lineCount = 0;

        for (std::string line; std::getline(infile, line);)
        {
            tree.ConvertWordToRtreeEntry(line, 0, lineCount);
            lineCount++;
        }

        tree.CloseConversion();
    }*/
};

/*void CheckForDuplicateIds(multiBoxGeo& boxes) {
    std::vector<unsigned int> cache;
    int i = 0;
    for (rTreeValue box : boxes) {
        bool found = false;
        for (unsigned int id : cache) {
            if (id == box.second) {
                std::cout << id << std::endl;
                i++;
                found = true;
            }
        }
        cache.push_back(box.second);
    }
    std::cout << "Found " << i << " duplicates" << std::endl;
}*/

void showCase() {
    GeoLocation test;
    //multiPolygonGeo m;
    //bg::read_wkt("MULTIPOLYGON(((6.0859436 50.7580645,6.0859881 50.7582046,6.0860374 50.7581984,6.0860460 50.7582255,6.0862678 50.7581973,6.0862594 50.7581708,6.0866540 50.7581207,6.0867201 50.7581123,6.0867354 50.7581603,6.0868121 50.7584019,6.0872680 50.7583439,6.0871550 50.7579880,6.0870223 50.7580049,6.0869902 50.7580090,6.0869831 50.7579865,6.0867916 50.7573839,6.0866235 50.7574053,6.0866417 50.7574623,6.0867314 50.7577444,6.0867755 50.7578832,6.0866545 50.7578986,6.0867048 50.7580568,6.0866710 50.7580611,6.0865367 50.7580782,6.0865118 50.7580000,6.0864775 50.7578917,6.0864619 50.7578427,6.0862767 50.7578662,6.0861611 50.7578809,6.0859437 50.7579086,6.0859913 50.7580585,6.0859436 50.7580645),(6.0868632 50.7581593,6.0870314 50.7581379,6.0870580 50.7582210,6.0868892 50.7582423,6.0868632 50.7581593)))", m);
    //std::cout << bg::wkt(test.convertWKTToBoundingBoxes("MULTIPOLYGON(((6.0859436 50.7580645,6.0859881 50.7582046,6.0860374 50.7581984,6.0860460 50.7582255,6.0862678 50.7581973,6.0862594 50.7581708,6.0866540 50.7581207,6.0867201 50.7581123,6.0867354 50.7581603,6.0868121 50.7584019,6.0872680 50.7583439,6.0871550 50.7579880,6.0870223 50.7580049,6.0869902 50.7580090,6.0869831 50.7579865,6.0867916 50.7573839,6.0866235 50.7574053,6.0866417 50.7574623,6.0867314 50.7577444,6.0867755 50.7578832,6.0866545 50.7578986,6.0867048 50.7580568,6.0866710 50.7580611,6.0865367 50.7580782,6.0865118 50.7580000,6.0864775 50.7578917,6.0864619 50.7578427,6.0862767 50.7578662,6.0861611 50.7578809,6.0859437 50.7579086,6.0859913 50.7580585,6.0859436 50.7580645),(6.0868632 50.7581593,6.0870314 50.7581379,6.0870580 50.7582210,6.0868892 50.7582423,6.0868632 50.7581593)))"));
    //std::cout << test.convertWKTToBoundingBoxes("MULTIPOLYGON(((1.0 1.0, 4.0 4.0), (2.0 2.0, 3.0 3.0)))").size();

    //std::cout << bg::wkt(test.convertWKTToBoundingBoxes("\"MULTIPOLYGON(((10.4254953 50.5721282,10.4256401 50.5721333,10.4256470 50.5720548,10.4255021 50.5720497,10.4254953 50.5721282)))\"^^<http://www.opengis.net/ont/geosparql#wktLiteral>")[0]);
    //test.convertWKTToBoundingBoxes("rtdfgdfgdrgrdg");

    /*test.createEntries({std::make_pair(GeoLocation::createBoundingBox(1.0, 1.0, 2.0, 2.0), 1),
                        std::make_pair(GeoLocation::createBoundingBox(3.0, 3.0, 4.0, 4.0), 2)});
    test.addEntry(std::make_pair(GeoLocation::createBoundingBox(2.0, 2.0, 3.0, 3.0), 3));
    test.Search(GeoLocation::createBoundingBox(2.9, 2.9, 5.0, 5.0));*/

    //test.convertWKTFile("../osm-germany_6tbDDu.tsv", "../germany_data_tidy.csv");
    //multiBoxGeo boxes = test.loadEntries("../data100k.csv");
    //test.createEntries(boxes);

    /*test.centerOrdering(boxes, 1);

    for (auto & box : boxes) {
        std::cout << bg::wkt(box.first) << std::endl;
        //std::cout << box.first.min_corner().get<0>() << std::endl;
    }*/

    //std::vector<multiBoxGeo> result = test.TGS(boxes, 316, 316);
    //test.saveTGSResult(result, "../100k");
    //std::cout << "finished";


    //auto startTime = std::chrono::high_resolution_clock::now();
    //test.Search(test.createBoundingBox(5.9204, 50.9949, 5.92056, 50.995), "../100k");
    //multiBoxGeo boxes = test.loadEntries("../germany_data_tidy.csv");
    //std::vector<multiBoxGeo> result = test.TGS(boxes, 6633, 6633);
    //test.saveTGSResult(result, "../germany");

    // new tests
    //multiBoxGeo boxes = test.loadEntries("../germany_data_tidy.csv");
    auto startTime = std::chrono::high_resolution_clock::now();

    //Rtree tree = Rtree();
    /*Node test_ = Node(7, test.createBoundingBox(0, 0, 1, 1));
    Node test_2 = Node(2, test.createBoundingBox(2, 2, 3, 3));
    Node test__2 = Node(42, test.createBoundingBox(123, 456, 789, 369));
    test_.AddChild(test_2);
    test_.AddChild(test__2);
    SaveNode(test_, false, "../test.bin");*/
    //Node test_3 = loadNode("../test.bin");
    /*multiBoxGeo boxes = tree.LoadEntries("../switzerland_raw");
    auto stopTime_ = std::chrono::high_resolution_clock::now();
    auto duration_ = std::chrono::duration_cast<std::chrono::microseconds>(stopTime_ - startTime_);
    std::cout << "Loaded in " << duration_.count() / 1000000.0 << " seconds" << std::endl;
    auto startTime = std::chrono::high_resolution_clock::now();
    tree.BuildTree(boxes, 16, "../switzerland_raw/rtree_build"); // took 1716.33 seconds for germany */
    //multiBoxGeo results = tree.SearchTree(test.createBoundingBox(5.9204, 50.9949, 5.92056, 50.995), "../germany_new");
    /*multiBoxGeo results = tree.SearchTree(test.createBoundingBox(7.73243, 45.2063, 7.73252, 45.2071), "../switzerland");
    for(rTreeValue result : results) {
        std::cout << result.first.min_corner().get<0>() << " " << result.first.min_corner().get<1>() << "," << result.first.max_corner().get<0>()
                  << " " << result.first.max_corner().get<1>() << "," << result.second << std::endl;
    }
    std::cout << "Found " << results.size() << " results:" << std::endl;*/
    //CheckForDuplicateIds(boxes);

    //test.TestConverter(tree, "../osm-germany-100k.tsv");
    /*multiBoxGeo results = tree.LoadEntries("../conversion_test_old");
    for(rTreeValue result : results) {
        std::cout << result.first.min_corner().get<0>() << " " << result.first.min_corner().get<1>() << "," << result.first.max_corner().get<0>()
                  << " " << result.first.max_corner().get<1>() << "," << result.second << std::endl;
    }
    std::cout << "Found " << results.size() << " results:" << std::endl;*/

    /*std::ifstream infile("../switzerland_raw/switzerland_raw_100k.txt");
    //std::filesystem::create_directory("../switzerland_raw");
    std::ofstream convertOfs = std::ofstream("../switzerland_raw/converted_data_100k", std::ios::binary);

    std::cout << "Loading" << std::endl;

    long long id = 1;
    for (std::string line; std::getline(infile, line);)
    {
        std::optional<boxGeo> boundingBox = Rtree::ConvertWordToRtreeEntry(line);
        if (boundingBox) {
            Rtree::SaveEntry(boundingBox.value(), id, convertOfs);
        }
        id++;
    }

    convertOfs.close();

    std::cout << "Loaded entries" << std::endl;
    std::cout << "Building Rtree" << std::endl;

    Rtree rtree = Rtree();
    multiBoxGeo entries = rtree.LoadEntries("../switzerland_raw/test");
    rtree.BuildTree(entries, 16, "../switzerland_raw/rtree_build");

    std::cout << "Finished building the Rtree with " << entries.size() << " entries" << std::endl;*/

    //std::cout << tree.time << std::endl;


    /* External Sorting */

    externalSort("../switzerland_raw/test", "../switzerland_raw/sorted_data", 5000000, 1);

    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    std::cout << "Searched in " << duration.count() / 1000000.0 << " seconds" << std::endl;
}