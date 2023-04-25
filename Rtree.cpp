//
// Created by Noah on 19.04.23.
//

#include "Rtree.h"

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
    } else {
        // order by centerY
        auto sortRuleLambda = [](rTreeValue b1, rTreeValue b2) -> bool {
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

static boxGeo createBoundingBox(double pointOneX, double pointOneY, double pointTwoX, double pointTwoY) {
    return make<boxGeo>(make<pointGeo>(pointOneX, pointOneY), make<pointGeo>(pointTwoX, pointTwoY));
}

static std::vector<std::vector<multiBoxGeo>> TGSRecursive(const std::vector<multiBoxGeo>& orderedInputRectangles, size_t M, size_t S) {  // https://dl.acm.org/doi/pdf/10.1145/288692.288723
    /*
     * inputRectangles needs to be pre-sorted with centerOrdering for both d0 and d1
     */

    size_t n = orderedInputRectangles[0].size();

    if (n <= S) {
        // stop condition
        return std::vector<std::vector<multiBoxGeo>> { orderedInputRectangles };
    }

    double bestCost = -1;
    size_t bestDim = 0;
    unsigned long bestI = 1;
    boxGeo bestB0 = createBoundingBox(0, 0, 0, 0);

    for (size_t dim = 0; dim <= 1; dim++) {
        for (size_t i = 1; i < M; i++) {
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
    std::vector<std::vector<multiBoxGeo>> result0 = TGSRecursive(s0, S, M);
    std::vector<std::vector<multiBoxGeo>> result1 = TGSRecursive(s1, S, M);

    std::vector<std::vector<multiBoxGeo>> result;
    result.insert(result.begin(), result0.begin(), result0.end());
    result.insert(result.end(), result1.begin(), result1.end());

    return result;
}

void Rtree::BuildTree(multiBoxGeo& inputRectangles, size_t M) {
    // sort the rectangles
    std::vector<multiBoxGeo> orderedInputRectangles;

    multiBoxGeo RectanglesD0(inputRectangles);
    multiBoxGeo RectanglesD1(inputRectangles);

    centerOrdering(RectanglesD0, 0);
    centerOrdering(RectanglesD1, 1);

    orderedInputRectangles.push_back(RectanglesD0);
    orderedInputRectangles.push_back(RectanglesD1);

    // build the tree in a depth first approach
    std::stack<ConstructionNode> layerStack;

    unsigned int newId = 1; // start from 1, because 0 is the root item
    ConstructionNode rootItem = ConstructionNode(0, orderedInputRectangles);
    layerStack.push(rootItem);

    while (!layerStack.empty()) {
        ConstructionNode currentItem = layerStack.top();
        layerStack.pop();

        if (currentItem.GetOrderedBoxes()[0].size() <= M) {
            // reached a leaf
            for(rTreeValue box : currentItem.GetOrderedBoxes()[0]) {
                Node leafNode = Node(box.second, box.first);
                currentItem.AddChild(leafNode);
            }
            currentItem.SaveAsLastInnerNode();
        } else {
            std::vector<std::vector<multiBoxGeo>> tgsResult = TGSRecursive(currentItem.GetOrderedBoxes(), M, std::ceil(currentItem.GetOrderedBoxes()[0].size() / M));
            for (std::vector<multiBoxGeo>& currentOrderedRectangles : tgsResult) {
                ConstructionNode newItem = ConstructionNode(newId, currentOrderedRectangles);
                layerStack.push(newItem);

                currentItem.AddChild(newItem);

                newId++;
            }

            currentItem.Save();
        }
    }
}

ConstructionNode::ConstructionNode(unsigned int id, std::vector<multiBoxGeo>& orderedBoxes)
: Node{id}
{
    this->orderedBoxes = orderedBoxes;

    // calculate the boundingBoxes
    double globalMinX = -1;
    double globalMinY = -1;
    double globalMaxX = -1;
    double globalMaxY = -1;

    for (rTreeValue box : orderedBoxes[0]) {
        double minX = box.first.min_corner().get<0>();
        double minY = box.first.min_corner().get<1>();
        double maxX = box.first.max_corner().get<0>();
        double maxY = box.first.max_corner().get<1>();

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

    this->boundingBox = createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
}

unsigned int Node::GetId() const {
    return this->id;
}

std::vector<multiBoxGeo> ConstructionNode::GetOrderedBoxes() {
    return this->orderedBoxes;
}

Node::Node(unsigned int id, boxGeo boundingBox) {
    this->id = id;
    this->boundingBox = boundingBox;
}

Node::Node(unsigned int id) {
    this->id = id;
}

void Node::AddChild(Node& child) {
    this->children.push_back(child);
}

void Node::Save() {
    std::cout << "Save" << std::endl;
}

void Node::SaveAsLastInnerNode() {
    std::cout << "SaveAsLastInnerNode" << std::endl;
}