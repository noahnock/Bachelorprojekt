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

using bg::make;

class Node {
protected:
    friend class boost::serialization::access;
    long long id;
    boxGeo boundingBox;
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
    long long GetId() const;
    boxGeo GetBoundingBox() const;
    void AddChild(Node& child);
    void SetIsLastInnerNode(bool isLastInnerNode);
    bool GetIsLastInnerNode();
    multiBoxGeo GetChildren();
};

BOOST_CLASS_VERSION(Node, 1)

class Rtree {
private:
    long long SaveNode(Node &node, bool isLastInnerNode, std::ofstream& nodesOfs);
    Node LoadNode(long long id, std::ifstream& lookupIfs, std::ifstream& nodesIfs);
    uintmax_t maxBuildingRamUsage;
public:
    void BuildTree(const std::string& onDiskBase, size_t M, const std::string& folder);
    multiBoxGeo SearchTree(boxGeo query, const std::string& folder);
    static std::optional<boxGeo> ConvertWordToRtreeEntry(const std::string& wkt);
    static void SaveEntry(boxGeo boundingBox, uint64_t index, std::ofstream& convertOfs);
    static multiBoxGeo LoadEntries(const std::string& file);
    static boxGeo createBoundingBox(double pointOneX, double pointOneY, double pointTwoX, double pointTwoY);
    explicit Rtree(uintmax_t maxBuildingRamUsage);
};

class OrderedBoxes {
public:
    bool workInRam;
    long long size;
    multiBoxGeo rectanglesD0InRam;
    multiBoxGeo rectanglesD1InRam;
    std::string rectanglesD0OnDisk;
    std::string rectanglesD1OnDisk;
public:
    [[nodiscard]] bool WorkInRam() const;
    void CreateOrderedBoxesInRam(multiBoxGeo& rectanglesD0, multiBoxGeo& rectanglesD1); // workInRam = true
    void CreateOrderedBoxesOnDisk(const std::string& rectanglesD0, const std::string& rectanglesD1, long long size); // workInRam = false
    boxGeo GetBoundingBox();
    long long GetSize() const;
    rTreeValue GetElementAt(size_t dim, long long index);
    std::pair<OrderedBoxes, OrderedBoxes> SplitAt(const std::string& filePath, size_t dim, long long index, boxGeo boundingBox, long long maxBuildingRamUsage);
};

class ConstructionNode: public Node {
private:
    OrderedBoxes orderedBoxes;

public:
    ConstructionNode(long long id, OrderedBoxes orderedBoxes);
    OrderedBoxes GetOrderedBoxes();
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
