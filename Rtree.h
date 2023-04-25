//
// Created by Noah on 19.04.23.
//

#ifndef BACHELORPROJEKT_RTREE_H
#define BACHELORPROJEKT_RTREE_H

#include <boost/geometry.hpp>
#include <vector>
#include <iostream>

namespace bg = boost::geometry;

using pointGeo = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
using boxGeo = bg::model::box<pointGeo>;
using rTreeValue = std::pair<boxGeo, unsigned int>;
using multiBoxGeo = std::vector<rTreeValue>;

using bg::make;

class Node {
protected:
    unsigned int id;
    boxGeo boundingBox{};
    std::vector<Node> children;
    explicit Node(unsigned int id);

public:
    Node(unsigned int id, boxGeo boundingBox);
    unsigned int GetId() const;
    void AddChild(Node& child);
    void Save();
    void SaveAsLastInnerNode();
};

class Rtree {
public:
    void BuildTree(multiBoxGeo& inputRectangles, size_t M);
};

class ConstructionNode: public Node {
private:
    std::vector<multiBoxGeo> orderedBoxes;

public:
    ConstructionNode(unsigned int id, std::vector<multiBoxGeo>& orderedBoxes);
    std::vector<multiBoxGeo> GetOrderedBoxes();
};


#endif //BACHELORPROJEKT_RTREE_H
