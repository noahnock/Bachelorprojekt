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
using rTreeValue = std::pair<boxGeo, unsigned int>;
using multiBoxGeo = std::vector<rTreeValue>;

using bg::make;

class Node {
protected:
    friend class boost::serialization::access;
    unsigned int id;
    boxGeo boundingBox;
    bool isLastInnerNode = false;
    std::vector<std::pair<unsigned int, boxGeo>> children;

    template<class Archive>
    void serialize(Archive & a, const unsigned int version) {
        a & id;
        a & isLastInnerNode;
        a & boundingBox;
        a & children;
    }

    explicit Node(unsigned int id);

public:
    Node();
    Node(unsigned int id, boxGeo boundingBox);
    Node(unsigned int id, boxGeo boundingBox, std::vector<std::pair<unsigned int, boxGeo>>& children, bool isLastInnerNode);
    Node(unsigned int id, double minX, double minY, double maxX, double maxY, bool isLastInnerNode);
    unsigned int GetId() const;
    boxGeo GetBoundingBox() const;
    void AddChild(Node& child);
    void SetIsLastInnerNode(bool isLastInnerNode);
};

BOOST_CLASS_VERSION(Node, 1)

class Rtree {
public:
    void BuildTree(multiBoxGeo& inputRectangles, size_t M);
};

class ConstructionNode: public Node {
private:
    std::vector<multiBoxGeo> orderedBoxes;

public:
    ConstructionNode(unsigned int id, std::vector<multiBoxGeo> orderedBoxes);
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

void SaveNode(Node &node, bool isLastInnerNode, const std::string& fileName);
Node loadNode(const std::string& fileName);

#endif //BACHELORPROJEKT_RTREE_H
