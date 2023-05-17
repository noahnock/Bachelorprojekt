//
// Created by Noah on 19.04.23.
//

#ifndef BACHELORPROJEKT_RTREE_H
#define BACHELORPROJEKT_RTREE_H

#include <boost/geometry.hpp>
#include <vector>
#include <iostream>
#include <fstream>
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
    void serialize(Archive & a, const unsigned int version) {
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
    std::ofstream nodesOfs;
    std::ifstream nodesIfs;
    std::ifstream lookupIfs;
    std::ofstream convertOfs;
    long long SaveNode(Node &node, bool isLastInnerNode);
    Node LoadNode(long long id);
public:
    void BuildTree(multiBoxGeo& inputRectangles, size_t M, const std::string& folder);
    multiBoxGeo SearchTree(boxGeo query, const std::string& folder);
    void OpenConversion(const std::string& folder);
    void CloseConversion();
    template <typename data_type>
    void ConvertWordToRtreeEntry(const data_type* data, size_t elementSize, uint64_t index);
};

class ConstructionNode: public Node {
private:
    std::vector<multiBoxGeo> orderedBoxes;

public:
    ConstructionNode(long long id, std::vector<multiBoxGeo> orderedBoxes);
    std::vector<multiBoxGeo> GetOrderedBoxes();
};

namespace boost::serialization {
    template<class Archive>
    void save(Archive & a, const boxGeo & b, unsigned int version)
    {
        a << b.min_corner().get<0>();
        a << b.min_corner().get<1>();
        a << b.max_corner().get<0>();
        a << b.max_corner().get<1>();
    }
    template<class Archive>
    void load(Archive & a, boxGeo & b, unsigned int version)
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
