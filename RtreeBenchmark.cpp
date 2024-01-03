//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include <iostream>
#include <vector>
#include <random>
#include "./Rtree/Rtree.h"
#include "Rtree/RtreeFileReader.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pointGeo = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
using polygonGeo = bg::model::polygon<pointGeo>;
using boxGeo = bg::model::box<pointGeo>;
using rTreeValue = std::pair<boxGeo, long long>;
using multiBox = std::vector<rTreeValue>;

void logMessageWithStream(const std::string& message, std::ofstream& output) {
    std::cout << message << std::endl;
    output << message << std::endl;
}

template<size_t M>
class BoostRtree {
private:
    bgi::rtree<rTreeValue , bgi::rstar<M>> rtree_;

public:
    void createEntries(const multiBox& entries) {
        std::cout << "Creating the entries in the R-tree..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();

        rtree_ = bgi::rtree<rTreeValue, bgi::rstar<M>>(entries.begin(), entries.end());

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << "Created the R-tree with " << entries.size() << " entries in " << duration.count() / 1000000.0 << " seconds" << std::endl;
    }

    multiBox SearchInTree(boxGeo query) {
        //std::cout << "Searching in the R-tree..." << std::endl;
        //auto startTime = std::chrono::high_resolution_clock::now();
        multiBox results;
        std::for_each(rtree_.qbegin(bgi::intersects(query)), rtree_.qend(), [&](rTreeValue const& result) {
            results.push_back(result);
        });

        /*auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
        std::cout << results.size() << " elements were found in " << duration.count() / 1000000.0 << " seconds" << std::endl;*/
        return results;
    }
};

template<size_t M>
std::pair<double, double> searchInBoostRtree(const std::string& file, std::vector<BasicGeometry::BoundingBox>& exampleQueries) {
    BoostRtree<M> tree;
    multiBoxGeo rectangles = FileReaderWithoutIndex::LoadEntries(file);
    multiBox entries;
    auto startTime = std::chrono::high_resolution_clock::now();
    for (RTreeValue entry : rectangles) {
        rTreeValue newEntry = std::make_pair(entry.box, entry.id);
        entries.push_back(newEntry);
    }
    tree.createEntries(entries);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);

    double totalTimeBuild = duration.count();
    double totalTimeSearch = 0;
    uint64_t totalSize = 0;

    for (size_t i = 0; i < 10; i++) {
        for (BasicGeometry::BoundingBox query : exampleQueries) {
            startTime = std::chrono::high_resolution_clock::now();
            multiBox results = tree.SearchInTree(query);
            stopTime = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
            totalTimeSearch += duration.count();
            totalSize += results.size();
        }
    }

    totalTimeSearch /= (double) totalSize / 1000.0;
    return std::make_pair(totalTimeBuild, totalTimeSearch);
}

double searchRtree(Rtree& rtree, std::vector<BasicGeometry::BoundingBox>& exampleQueries) {
    double totalTimeSearch = 0;
    uint64_t totalSize = 0;

    for (size_t i = 0; i < 10; i++) {
        for (BasicGeometry::BoundingBox query : exampleQueries) {
            auto startTime = std::chrono::high_resolution_clock::now();
            multiBoxGeo results = rtree.SearchTree(query);
            auto stopTime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
            totalTimeSearch += (double) duration.count();
            totalSize += results.size();
        }
    }

    totalTimeSearch /= (double) totalSize / 1000.0;
    return totalTimeSearch;
}

void createUniformSample(size_t n, double maxRelativeSize, const std::string& suffix) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distributionX(5.956114,10.492341);
    std::uniform_real_distribution<double> distributionY(45.817943,47.808457);
    double maxArea = maxRelativeSize * ((10.492341 - 5.956114) * (47.808457 - 45.817943));

    // pick one side random but restrict its size
    std::uniform_real_distribution<double> lengthXAbs(0.1 * maxArea,0.9 * maxArea);
    // this will be the percentage of the possible maximum length for Y
    std::uniform_real_distribution<double> lengthYPerc(0.0, 1.0);

    std::ofstream sampleOfs = std::ofstream("../Benchmarks/uniform_sample" + suffix + ".tmp", std::ios::binary);

    for (size_t i = 0; i < n; i++) {
        double x = distributionX(generator);
        double y = distributionY(generator);
        double lengthX = lengthXAbs(generator);
        double xHalf = lengthX / 2;
        double lengthY = lengthYPerc(generator) * (maxArea / lengthX);
        double yHalf = lengthY / 2;

        BasicGeometry::BoundingBox box = BasicGeometry::CreateBoundingBox(x - xHalf, y - yHalf, x + xHalf, y + yHalf);
        FileReaderWithoutIndex::SaveEntry(box, i, sampleOfs);
    }

    sampleOfs.close();
}

void createNormalDistSample(size_t n, double maxRelativeSize, const std::string& suffix) {
    std::default_random_engine generator;
    std::normal_distribution<double> distributionX((10.492341 + 5.956114) / 2.0,(10.492341 - 5.956114) / 4.0);
    std::normal_distribution<double> distributionY((47.808457 + 45.817943) / 2.0,(47.808457 - 45.817943) / 4.0);
    double maxArea = maxRelativeSize * ((10.492341 - 5.956114) * (47.808457 - 45.817943));

    // pick one side random but restrict its size
    std::normal_distribution<double> lengthXAbs((0.9 * maxArea + 0.1 * maxArea) / 2.0, (0.9 * maxArea - 0.1 * maxArea) / 4.0);
    // this will be the percentage of the possible maximum length for Y
    std::uniform_real_distribution<double> lengthYPerc(0.0, 1.0);

    std::ofstream sampleOfs = std::ofstream("../Benchmarks/normal_sample" + suffix + ".tmp", std::ios::binary);

    for (size_t i = 0; i < n; i++) {
        double x = distributionX(generator);
        double y = distributionY(generator);
        double lengthX = lengthXAbs(generator);
        double xHalf = lengthX / 2;
        double lengthY = lengthYPerc(generator) * (maxArea / lengthX);
        double yHalf = lengthY / 2;

        BasicGeometry::BoundingBox box = BasicGeometry::CreateBoundingBox(x - xHalf, y - yHalf, x + xHalf, y + yHalf);
        FileReaderWithoutIndex::SaveEntry(box, i, sampleOfs);
    }

    sampleOfs.close();
}

std::vector<BasicGeometry::BoundingBox> exampleQueries {
    BasicGeometry::CreateBoundingBox(9.88657, 47.38431, 9.88671, 47.6088),  // random
    BasicGeometry::CreateBoundingBox(5.956142, 46.132035, 5.956142, 46.132035),  // most western point
    BasicGeometry::CreateBoundingBox(7.365485, 46.904082, 7.536615, 46.990009),  // Bern
    BasicGeometry::CreateBoundingBox(9.234133, 46.346653, 9.515004, 46.539109),  // random
    BasicGeometry::CreateBoundingBox(7.969524, 45.912223, 8.644526, 47.847424),  // stripe vertical
    BasicGeometry::CreateBoundingBox(6.087310, 46.547426, 10.519028, 46.904035), // stripe horizontal
    BasicGeometry::CreateBoundingBox(9.467977, 47.245444, 9.467977, 47.245444), // random point
    BasicGeometry::CreateBoundingBox(7.57550, 47.54393, 7.58432, 47.54927),
    BasicGeometry::CreateBoundingBox(9.07570, 47.27046, 9.31603, 47.42956),
    BasicGeometry::CreateBoundingBox(9.44316, 46.78160, 10.74229, 47.55278),
    BasicGeometry::CreateBoundingBox(7.88828, 46.57483, 8.08329, 46.63804),
    BasicGeometry::CreateBoundingBox(6.87496, 46.37226, 9.467835, 47.490016),
    BasicGeometry::CreateBoundingBox(9.692774, 46.272107, 10.470066, 47.061925)
};

void Benchmark() {
    const std::string onDiskBase = "../switzerland_raw/converted_data_100k";
    const std::string sourceFile = onDiskBase + ".boundingbox.tmp";
    size_t sampleSize = 100000;

    // preparation
    createUniformSample(sampleSize, 0.001, ".001");
    createUniformSample(sampleSize, 0.0001, ".0001");
    createNormalDistSample(sampleSize, 0.001, ".001");
    createNormalDistSample(sampleSize, 0.0001, ".0001");

    createUniformSample(100000, 0.0001, ".100k");
    createUniformSample(200000, 0.0001, ".200k");
    createUniformSample(300000, 0.0001, ".300k");
    createUniformSample(400000, 0.0001, ".400k");
    createUniformSample(500000, 0.0001, ".500k");
    createUniformSample(600000, 0.0001, ".600k");
    createUniformSample(700000, 0.0001, ".700k");
    createUniformSample(800000, 0.0001, ".800k");
    createUniformSample(900000, 0.0001, ".900k");
    createUniformSample(1000000, 0.0001, ".1M");

    // building own R-tree
    std::ofstream output_building = std::ofstream("../Benchmarks/benchmark_building.txt");
    auto logMessage = [&output_building](const std::string& message) { logMessageWithStream(message, output_building); };
    logMessage("Building my own R-tree with different M");

    Rtree rtree2 = Rtree(99999999999);
    auto startTime = std::chrono::high_resolution_clock::now();
    uint64_t numElements = rtree2.BuildTree(onDiskBase, ".boundingbox", 2, "../Benchmarks/rtree_build_2");
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage(std::to_string(numElements) + " elements were found:");
    logMessage("M = 2: Built in " + std::to_string(duration.count()));

    Rtree rtree4 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree4.BuildTree(onDiskBase, ".boundingbox", 4, "../Benchmarks/rtree_build_4");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 4: Built in " + std::to_string(duration.count()));

    Rtree rtree8 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree8.BuildTree(onDiskBase, ".boundingbox", 8, "../Benchmarks/rtree_build_8");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 8: Built in " + std::to_string(duration.count()));

    Rtree rtree16 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree16.BuildTree(onDiskBase, ".boundingbox", 16, "../Benchmarks/rtree_build_16");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 16: Built in " + std::to_string(duration.count()));

    Rtree rtree32 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree32.BuildTree(onDiskBase, ".boundingbox", 32, "../Benchmarks/rtree_build_32");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 32: Built in " + std::to_string(duration.count()));

    Rtree rtree64Disk = Rtree(500000000);
    startTime = std::chrono::high_resolution_clock::now();
    rtree64Disk.BuildTree(onDiskBase, ".boundingbox", 64, "../Benchmarks/rtree_build_64");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 64 on disk: Built in " + std::to_string(duration.count()));

    Rtree rtree64 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree64.BuildTree(onDiskBase, ".boundingbox", 64, "../Benchmarks/rtree_build_64");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 64: Built in " + std::to_string(duration.count()));

    Rtree rtree128 = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree128.BuildTree(onDiskBase, ".boundingbox", 128, "../Benchmarks/rtree_build_128");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessage("M = 128: Built in " + std::to_string(duration.count()));

    output_building.close();

    // R-tree search
    std::ofstream output_search = std::ofstream("../Benchmarks/benchmark_search.txt");
    auto logMessageSearch = [&output_search](const std::string& message) { logMessageWithStream(message, output_search); };
    logMessageSearch("Search my own R-tree with different M");

    rtree2.SetupForSearch("../Benchmarks/rtree_build_2");
    double rtreeSearch2 = searchRtree(rtree2, exampleQueries);
    logMessageSearch("M = 2: Searched in " + std::to_string(rtreeSearch2));

    rtree4.SetupForSearch("../Benchmarks/rtree_build_4");
    double rtreeSearch4 = searchRtree(rtree4, exampleQueries);
    logMessageSearch("M = 4: Searched in " + std::to_string(rtreeSearch4));

    rtree8.SetupForSearch("../Benchmarks/rtree_build_8");
    double rtreeSearch8 = searchRtree(rtree8, exampleQueries);
    logMessageSearch("M = 8: Searched in " + std::to_string(rtreeSearch8));

    rtree16.SetupForSearch("../Benchmarks/rtree_build_16");
    double rtreeSearch16 = searchRtree(rtree16, exampleQueries);
    logMessageSearch("M = 16: Searched in " + std::to_string(rtreeSearch16));

    rtree32.SetupForSearch("../Benchmarks/rtree_build_32");
    double rtreeSearch32 = searchRtree(rtree32, exampleQueries);
    logMessageSearch("M = 32: Searched in " + std::to_string(rtreeSearch32));

    rtree64.SetupForSearch("../Benchmarks/rtree_build_64");
    double rtreeSearch64 = searchRtree(rtree64, exampleQueries);
    logMessageSearch("M = 64: Searched in " + std::to_string(rtreeSearch64));

    rtree128.SetupForSearch("../Benchmarks/rtree_build_128");
    double rtreeSearch128 = searchRtree(rtree128, exampleQueries);
    logMessageSearch("M = 128: Searched in " + std::to_string(rtreeSearch128));

    output_search.close();

    // build of different sizes
    std::ofstream output_sizes = std::ofstream("../Benchmarks/benchmark_sizes.txt");
    auto logMessageSizes = [&output_sizes](const std::string& message) { logMessageWithStream(message, output_sizes); };
    logMessageSizes("Building the R-tree on different dataset sizes");

    Rtree rtree100k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree100k.BuildTree("../Benchmarks/uniform_sample", ".100k", 64, "../Benchmarks/rtree_build_100k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("100k Built in " + std::to_string(duration.count()));

    Rtree rtree200k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree200k.BuildTree("../Benchmarks/uniform_sample", ".200k", 64, "../Benchmarks/rtree_build_200k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("200k Built in " + std::to_string(duration.count()));

    Rtree rtree300k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree300k.BuildTree("../Benchmarks/uniform_sample", ".300k", 64, "../Benchmarks/rtree_build_300k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("300k Built in " + std::to_string(duration.count()));

    Rtree rtree400k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree400k.BuildTree("../Benchmarks/uniform_sample", ".400k", 64, "../Benchmarks/rtree_build_400k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("400k Built in " + std::to_string(duration.count()));

    Rtree rtree500k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree500k.BuildTree("../Benchmarks/uniform_sample", ".500k", 64, "../Benchmarks/rtree_build_500k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("500k Built in " + std::to_string(duration.count()));

    Rtree rtree600k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree600k.BuildTree("../Benchmarks/uniform_sample", ".600k", 64, "../Benchmarks/rtree_build_600k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("600k Built in " + std::to_string(duration.count()));

    Rtree rtree700k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree700k.BuildTree("../Benchmarks/uniform_sample", ".700k", 64, "../Benchmarks/rtree_build_700k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("700k Built in " + std::to_string(duration.count()));

    Rtree rtree800k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree800k.BuildTree("../Benchmarks/uniform_sample", ".800k", 64, "../Benchmarks/rtree_build_800k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("800k Built in " + std::to_string(duration.count()));

    Rtree rtree900k = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree900k.BuildTree("../Benchmarks/uniform_sample", ".900k", 64, "../Benchmarks/rtree_build_900k");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("900k Built in " + std::to_string(duration.count()));

    Rtree rtree1M = Rtree(99999999999);
    startTime = std::chrono::high_resolution_clock::now();
    rtree1M.BuildTree("../Benchmarks/uniform_sample", ".1M", 64, "../Benchmarks/rtree_build_1M");
    stopTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    logMessageSizes("1M Built in " + std::to_string(duration.count()));

    output_sizes.close();

    // sample data boost vs mine
    /*std::ofstream output_sample = std::ofstream("../Benchmarks/benchmark_sample.txt");
    auto logMessageSample = [&output_sample](const std::string& message) { logMessageWithStream(message, output_sample); };
    logMessageSample("Comparing boost and mine on different generated sample datasets");

    Rtree rtreeUn001 = Rtree(99999999999);
    rtreeUn001.BuildTree("../Benchmarks/uniform_sample", ".001", 64, "../Benchmarks/rtree_build_Un001");
    rtreeUn001.SetupForSearch("../Benchmarks/rtree_build_Un001");
    double rtreeUn001Res = searchRtree(rtreeUn001, exampleQueries);
    logMessageSample("Uniform 0.001 with my R-tree: " + std::to_string(rtreeUn001Res));

    Rtree rtreeUn0001 = Rtree(99999999999);
    rtreeUn0001.BuildTree("../Benchmarks/uniform_sample", ".0001", 64, "../Benchmarks/rtree_build_Un0001");
    rtreeUn0001.SetupForSearch("../Benchmarks/rtree_build_Un0001");
    double rtreeUn0001Res = searchRtree(rtreeUn0001, exampleQueries);
    logMessageSample("Uniform 0.0001 with my R-tree: " + std::to_string(rtreeUn0001Res));

    Rtree rtreeNorm001 = Rtree(99999999999);
    rtreeNorm001.BuildTree("../Benchmarks/normal_sample", ".001", 64, "../Benchmarks/rtree_build_Norm001");
    rtreeNorm001.SetupForSearch("../Benchmarks/rtree_build_Norm001");
    double rtreeNorm001Res = searchRtree(rtreeNorm001, exampleQueries);
    logMessageSample("Normal 0.001 with my R-tree: " + std::to_string(rtreeNorm001Res));

    Rtree rtreeNorm0001 = Rtree(99999999999);
    rtreeNorm0001.BuildTree("../Benchmarks/normal_sample", ".0001", 64, "../Benchmarks/rtree_build_Norm0001");
    rtreeNorm0001.SetupForSearch("../Benchmarks/rtree_build_Norm0001");
    double rtreeNorm0001Res = searchRtree(rtreeNorm0001, exampleQueries);
    logMessageSample("Normal 0.0001 with my R-tree: " + std::to_string(rtreeNorm0001Res));

    std::pair<double, double> boostUn001 = searchInBoostRtree<64>("../Benchmarks/uniform_sample.001.tmp", exampleQueries);
    logMessageSample("Uniform 0.001 with boost R-tree: " + std::to_string(boostUn001.second));

    std::pair<double, double> boostUn0001 = searchInBoostRtree<64>("../Benchmarks/uniform_sample.0001.tmp", exampleQueries);
    logMessageSample("Uniform 0.0001 with boost R-tree: " + std::to_string(boostUn0001.second));

    std::pair<double, double> boostNorm001 = searchInBoostRtree<64>("../Benchmarks/normal_sample.001.tmp", exampleQueries);
    logMessageSample("Normal 0.001 with boost R-tree: " + std::to_string(boostNorm001.second));

    std::pair<double, double> boostNorm0001 = searchInBoostRtree<64>("../Benchmarks/normal_sample.0001.tmp", exampleQueries);
    logMessageSample("Normal 0.0001 with boost R-tree: " + std::to_string(boostNorm0001.second));

    output_sample.close();*/

    // boost R-tree
    std::ofstream output_boost = std::ofstream("../Benchmarks/benchmark_boost.txt");
    auto logMessageBoost = [&output_boost](const std::string& message) { logMessageWithStream(message, output_boost); };
    logMessageBoost("Building the boost R-tree with M=16");

    std::pair<double, double> boost64 = searchInBoostRtree<64>(sourceFile, exampleQueries);
    logMessageBoost("M = 64: Built in " + std::to_string(boost64.first) + " and searched in: " + std::to_string(boost64.second));

    output_boost.close();
}
