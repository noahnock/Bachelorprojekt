//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#ifndef BACHELORPROJEKT_RTREE_H
#define BACHELORPROJEKT_RTREE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <optional>
#include <boost/geometry.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

namespace bg = boost::geometry;

using pointGeo = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
using boxGeo = bg::model::box<pointGeo>;
using rTreeValue = std::pair<boxGeo, long long>;
using multiBoxGeo = std::vector<rTreeValue>;
using rTreeValueWithOrderIndex = std::pair<rTreeValue, std::pair<long long, long long>>;
using multiBoxWithOrderIndex = std::vector<rTreeValueWithOrderIndex>;

using bg::make;

struct SplitResult {
    double bestCost = -1;
    size_t bestDim = 0;
    long long bestIndex = 0;
    rTreeValueWithOrderIndex bestLastElement;
    rTreeValueWithOrderIndex bestElement;
    rTreeValueWithOrderIndex bestMinElement;
    rTreeValueWithOrderIndex bestMaxElement;
};

struct SplitBuffersRam {
    std::shared_ptr<multiBoxWithOrderIndex> s0Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s0Dim1 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1Dim1 = std::make_shared<multiBoxWithOrderIndex>();

    std::shared_ptr<multiBoxWithOrderIndex> s0SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s0SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();
};

struct SplitBuffersDisk {
    SplitBuffersRam splitBuffersRam;
    std::optional<std::ofstream> split0Dim0File;
    std::optional<std::ofstream> split0Dim1File;
    std::optional<std::ofstream> split1Dim0File;
    std::optional<std::ofstream> split1Dim1File;
};

class Node {
protected:
    friend class boost::serialization::access;
    long long id;
    boxGeo boundingBox{};
    bool isLastInnerNode = false;
    multiBoxGeo children;

    template<class Archive>
    void serialize(Archive & a, [[maybe_unused]]const unsigned int version) {
        a & id;
        a & isLastInnerNode;
        a & boundingBox;
        a & children;
    }

    explicit Node(long long id);

public:
    Node();
    Node(long long id, boxGeo boundingBox);
    Node(long long id, boxGeo boundingBox, multiBoxGeo &children, bool isLastInnerNode);
    Node(long long id, double minX, double minY, double maxX, double maxY, bool isLastInnerNode);
    [[nodiscard]] long long GetId() const;
    [[nodiscard]] boxGeo GetBoundingBox() const;
    void AddChild(Node& child);
    void SetIsLastInnerNode(bool isLastInnerNode);
    [[nodiscard]] bool GetIsLastInnerNode() const;
    multiBoxGeo GetChildren();
};

BOOST_CLASS_VERSION(Node, 1)

class Rtree {
private:
    static long long SaveNode(Node &node, bool isLastInnerNode, std::ofstream& nodesOfs);
    static Node LoadNode(long long id, std::ifstream& lookupIfs, std::ifstream& nodesIfs);
    uintmax_t maxBuildingRamUsage;
public:
    void BuildTree(const std::string& onDiskBase, size_t M, const std::string& folder) const;
    static multiBoxGeo SearchTree(boxGeo query, const std::string& folder);
    static std::optional<boxGeo> ConvertWordToRtreeEntry(const std::string& wkt);
    static void SaveEntry(boxGeo boundingBox, uint64_t index, std::ofstream& convertOfs);
    static void SaveEntryWithOrderIndex(rTreeValueWithOrderIndex treeValue, std::ofstream& convertOfs);
    static multiBoxGeo LoadEntries(const std::string& file);
    static multiBoxWithOrderIndex LoadEntriesWithOrderIndex(const std::string& file);
    static boxGeo createBoundingBox(double pointOneX, double pointOneY, double pointTwoX, double pointTwoY);
    explicit Rtree(uintmax_t maxBuildingRamUsage);
};

class OrderedBoxes {
private:
    bool workInRam;
    long long size;
    boxGeo boundingBox;
    std::shared_ptr<multiBoxWithOrderIndex> rectanglesD0InRam;
    std::shared_ptr<multiBoxWithOrderIndex> rectanglesD1InRam;
    std::string rectanglesD0OnDisk;
    std::string rectanglesD1OnDisk;
    std::shared_ptr<multiBoxWithOrderIndex> rectanglesD0Small;
    std::shared_ptr<multiBoxWithOrderIndex> rectanglesD1Small;
    std::pair<OrderedBoxes, OrderedBoxes> SplitAtBestInRam(size_t S, size_t M);
    std::pair<OrderedBoxes, OrderedBoxes> SplitAtBestOnDisk(const std::string& filePath, size_t S, size_t M, long long maxBuildingRamUsage);
    SplitResult GetBestSplit();
    std::pair<boxGeo, boxGeo> PerformSplit(SplitResult splitResult, SplitBuffersRam& splitBuffersRam, size_t M, size_t S);
    std::pair<boxGeo, boxGeo> PerformSplit(SplitResult splitResult, SplitBuffersDisk& splitBuffers, size_t M, size_t S, long long maxBuildingRamUsage);
public:
    [[nodiscard]] bool WorkInRam() const;
    void CreateOrderedBoxesInRam(const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesD1, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD1, boxGeo box); // workInRam = true
    void CreateOrderedBoxesOnDisk(const std::string& rectanglesD0, const std::string& rectanglesD1, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD1, long long size, boxGeo box); // workInRam = false
    boxGeo GetBoundingBox();
    [[nodiscard]] long long GetSize() const;
    std::pair<OrderedBoxes, OrderedBoxes> SplitAtBest(const std::string& filePath, size_t S, size_t M, long long maxBuildingRamUsage);
    std::shared_ptr<multiBoxWithOrderIndex> GetRectanglesInRam();
    std::string GetRectanglesOnDisk();
};

class ConstructionNode: public Node {
private:
    OrderedBoxes orderedBoxes;

public:
    ConstructionNode(long long id, OrderedBoxes orderedBoxes);
    OrderedBoxes GetOrderedBoxes();
    void AddChildrenToItem();
};

class FileReader {
private:
    std::string filePath;
    std::ifstream file;
    long long fileLength;
public:
    explicit FileReader(const std::string& filePath);
    std::optional<rTreeValueWithOrderIndex> GetNextElement();
    void Close();
};

namespace boost::serialization {
    template<class Archive>
    void save(Archive & a, const boxGeo & b, [[maybe_unused]]unsigned int version)
    {
        a << b.min_corner().get<0>();
        a << b.min_corner().get<1>();
        a << b.max_corner().get<0>();
        a << b.max_corner().get<1>();
    }
    template<class Archive>
    void load(Archive & a, boxGeo & b, [[maybe_unused]]unsigned int version)
    {
        double minX = 0;
        a >> minX;
        double minY = 0;
        a >> minY;
        double maxX = 0;
        a >> maxX;
        double maxY = 0;
        a >> maxY;
        b = make<boxGeo>(make<pointGeo>(minX, minY), make<pointGeo>(maxX, maxY));
    }
}
BOOST_SERIALIZATION_SPLIT_FREE(boxGeo);

#endif //BACHELORPROJEKT_RTREE_H
