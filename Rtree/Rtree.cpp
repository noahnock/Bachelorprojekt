//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include "Rtree.h"

bool intersects(const boxGeo &b1, const boxGeo &b2) {
    bool notIntersecting = b1.min_corner().get<0>() > b2.max_corner().get<0>() ||
                           b2.min_corner().get<0>() > b1.max_corner().get<0>() ||
                           b1.min_corner().get<1>() > b2.max_corner().get<1>() ||
                           b2.min_corner().get<1>() > b1.max_corner().get<1>();

    return !notIntersecting;
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

static void centerOrdering(multiBoxWithOrderIndex& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        auto sortRuleLambda = [] (rTreeValueWithOrderIndex b1, rTreeValueWithOrderIndex b2) -> bool
        {
            double center1 = (b1.first.first.min_corner().get<0>() + b1.first.first.max_corner().get<0>()) / 2;
            double center2 = (b2.first.first.min_corner().get<0>() + b2.first.first.max_corner().get<0>()) / 2;

            if (b1.second.first == b2.second.first)
                return center1 < center2;
            return b1.second.first < b2.second.first;
        };

        std::sort(boxes.begin(), boxes.end(), sortRuleLambda);
    } else {
        // order by centerY
        auto sortRuleLambda = [](rTreeValueWithOrderIndex b1, rTreeValueWithOrderIndex b2) -> bool {
            double center1 = (b1.first.first.min_corner().get<1>() + b1.first.first.max_corner().get<1>()) / 2;
            double center2 = (b2.first.first.min_corner().get<1>() + b2.first.first.max_corner().get<1>()) / 2;

            if (b1.second.second == b2.second.second)
                return center1 < center2;
            return b1.second.second < b2.second.second;
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

static bool boxInMultiBox(rTreeValueWithOrderIndex& item, std::shared_ptr<multiBoxWithOrderIndex> boxes) {
    for (rTreeValueWithOrderIndex currentItem : *boxes) {
        boxGeo currentBox = currentItem.first.first;
        boxGeo box = item.first.first;
        if (currentItem.second == item.second && box.min_corner().get<0>() == currentBox.min_corner().get<0>()
            && box.min_corner().get<1>() == currentBox.min_corner().get<1>()
            && box.max_corner().get<0>() == currentBox.max_corner().get<0>()
            && box.max_corner().get<1>() == currentBox.max_corner().get<1>()) {
            return true;
        }
    }
    return false;
}

boxGeo Rtree::createBoundingBox(double pointOneX, double pointOneY, double pointTwoX, double pointTwoY) {
    return make<boxGeo>(make<pointGeo>(pointOneX, pointOneY), make<pointGeo>(pointTwoX, pointTwoY));
}

static std::vector<OrderedBoxes> TGSRecursive(const std::string& filePath, OrderedBoxes orderedInputRectangles, size_t M, size_t S, long long maxBuildingRamUsage) {  // https://dl.acm.org/doi/pdf/10.1145/288692.288723
    /*
     * inputRectangles needs to be pre-sorted with centerOrdering for both d0 and d1
     */

    unsigned long long n = orderedInputRectangles.GetSize();

    if (n <= S || n <= M) {
        // stop condition
        return std::vector<OrderedBoxes> { orderedInputRectangles };
    }
    // split the rectangles at the best split
    std::pair<OrderedBoxes, OrderedBoxes> split = orderedInputRectangles.SplitAtBest(filePath, S, M, maxBuildingRamUsage);

    if (split.first.GetSize() % S != 0) {
        int fdgfg = 0;
    }

    // recursion
    std::vector<OrderedBoxes> result0 = TGSRecursive(filePath + ".0", split.first, M, S, maxBuildingRamUsage);
    std::vector<OrderedBoxes> result1 = TGSRecursive(filePath + ".1", split.second, M, S, maxBuildingRamUsage);

    std::vector<OrderedBoxes> result;
    result.insert(result.begin(), result0.begin(), result0.end());
    result.insert(result.end(), result1.begin(), result1.end());

    return result;
}

void Rtree::BuildTree(const std::string& onDiskBase, size_t M, const std::string& folder) {
    const std::string file = onDiskBase + ".boundingbox.tmp";

    // prepare the files
    std::filesystem::create_directory(folder);
    std::ofstream nodesOfs = std::ofstream(folder + "/nodes.bin", std::ios::binary);
    std::map<long long, long long> lookup;

    // sort the rectangles
    OrderedBoxes orderedInputRectangles;
    long long fileLines = std::ceil(std::filesystem::file_size(file) / (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long)));
    if ((std::filesystem::file_size(file) + fileLines * 2 * sizeof(long long)) * 4 < this->maxBuildingRamUsage) {
        std::cout << "Building in ram" << std::endl;
        // do it in ram
        multiBoxGeo RectanglesD0 = LoadEntries(file);
        centerOrdering(RectanglesD0, 0);

        double globalMinX = -1;
        double globalMinY = -1;
        double globalMaxX = -1;
        double globalMaxY = -1;

        size_t currentS = std::ceil(((float) RectanglesD0.size()) / ((float) M));
        //size_t iLimit = std::ceil(RectanglesD0.size() / currentS);

        std::shared_ptr<multiBoxWithOrderIndex> R0Small = std::make_shared<multiBoxWithOrderIndex>();
        std::shared_ptr<multiBoxWithOrderIndex> R1Small = std::make_shared<multiBoxWithOrderIndex>();

        std::shared_ptr<multiBoxWithOrderIndex> RectanglesD1WithOrder = std::make_shared<multiBoxWithOrderIndex>();
        for (long long i = 0; i < RectanglesD0.size(); i++) {
            rTreeValue element = RectanglesD0[i];
            rTreeValueWithOrderIndex entry = std::make_pair(element, std::make_pair(i, 0));
            RectanglesD1WithOrder->push_back(entry);

            if (globalMinX == -1 || element.first.min_corner().get<0>() < globalMinX) {
                globalMinX = element.first.min_corner().get<0>();
            }
            if (globalMinY == -1 || element.first.min_corner().get<1>() < globalMinY) {
                globalMinY = element.first.min_corner().get<1>();
            }
            if (element.first.max_corner().get<0>() > globalMaxX) {
                globalMaxX = element.first.max_corner().get<0>();
            }
            if (element.first.max_corner().get<1>() > globalMaxY) {
                globalMaxY = element.first.max_corner().get<1>();
            }
        }

        centerOrdering(*RectanglesD1WithOrder, 1);
        R1Small->push_back((*RectanglesD1WithOrder)[0]);
        rTreeValueWithOrderIndex maxElementDim1 = (*RectanglesD1WithOrder)[RectanglesD1WithOrder->size() - 1];
        maxElementDim1.second.second = RectanglesD1WithOrder->size() - 1;
        R1Small->push_back(maxElementDim1);
        for (long long i = 0; i < RectanglesD1WithOrder->size(); i++) {
            (*RectanglesD1WithOrder)[i].second.second = i;

            if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1/* && (i + 1) / currentS < iLimit*/)
                || (i % currentS == 0 && i / currentS >= 1/* && i / currentS < iLimit*/)) {
                // index i * S - 1 or i * S
                R1Small->push_back((*RectanglesD1WithOrder)[i]);
            }
        }

        std::shared_ptr<multiBoxWithOrderIndex> RectanglesD0WithOrder = std::make_shared<multiBoxWithOrderIndex>(*RectanglesD1WithOrder);
        centerOrdering(*RectanglesD0WithOrder, 0);

        R0Small->push_back((*RectanglesD0WithOrder)[0]);
        rTreeValueWithOrderIndex maxElementDim0 = (*RectanglesD0WithOrder)[RectanglesD0WithOrder->size() - 1];
        maxElementDim0.second.second = RectanglesD0WithOrder->size() - 1;
        R0Small->push_back(maxElementDim0);
        for (long long i = 0; i < RectanglesD0WithOrder->size(); i++) {
            //(*RectanglesD0WithOrder)[i].second.second = i;

            if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1/* && (i + 1) / currentS < iLimit*/)
                || (i % currentS == 0 && i / currentS >= 1/* && i / currentS < iLimit*/)) {
                // index i * S - 1 or i * S
                R0Small->push_back((*RectanglesD0WithOrder)[i]);
            }
        }

        orderedInputRectangles.CreateOrderedBoxesInRam(RectanglesD0WithOrder, RectanglesD1WithOrder, R0Small, R1Small,
                                                       createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY));
        std::cout << "Finished initial sorting" << std::endl;
    } else {
        std::cout << "Building on disk" << std::endl;
        // do it on disk

        // TODO sort on disk instead of ram
        // START
        multiBoxGeo RectanglesD0 = LoadEntries(file);

        centerOrdering(RectanglesD0, 0);

        long long xSize = 0;
        double globalMinX = -1;
        double globalMinY = -1;
        double globalMaxX = -1;
        double globalMaxY = -1;

        std::ofstream r0File = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
        for (rTreeValue element : RectanglesD0) {
            rTreeValueWithOrderIndex entry = std::make_pair(element, std::make_pair(xSize, 0));
            Rtree::SaveEntryWithOrderIndex(entry, r0File);
            xSize++;

            if (globalMinX == -1 || element.first.min_corner().get<0>() < globalMinX) {
                globalMinX = element.first.min_corner().get<0>();
            }
            if (globalMinY == -1 || element.first.min_corner().get<1>() < globalMinY) {
                globalMinY = element.first.min_corner().get<1>();
            }
            if (element.first.max_corner().get<0>() > globalMaxX) {
                globalMaxX = element.first.max_corner().get<0>();
            }
            if (element.first.max_corner().get<1>() > globalMaxY) {
                globalMaxY = element.first.max_corner().get<1>();
            }
        }
        r0File.close();
        RectanglesD0.clear();

        multiBoxWithOrderIndex RectanglesD1 = LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d0.tmp");
        centerOrdering(RectanglesD1, 1);

        size_t currentS = std::ceil(((float) xSize) / ((float) M));
        size_t iLimit = std::ceil(xSize / currentS);

        long long ySize = 0;
        std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
        std::shared_ptr<multiBoxWithOrderIndex> r1Small = std::make_shared<multiBoxWithOrderIndex>();
        r1Small->push_back(RectanglesD1[0]);
        rTreeValueWithOrderIndex maxElementDim1 = RectanglesD1[RectanglesD1.size() - 1];
        maxElementDim1.second.second = RectanglesD1.size() - 1;
        r1Small->push_back(maxElementDim1);
        for (rTreeValueWithOrderIndex element : RectanglesD1) {
            element.second.second = ySize;
            Rtree::SaveEntryWithOrderIndex(element, r1File);

            if (((ySize + 1) % currentS == 0 && (ySize + 1) / currentS >= 1 && (ySize + 1) / currentS < iLimit)
            || (ySize % currentS == 0 && ySize / currentS >= 1 && ySize / currentS < iLimit)) {
                // index i * S - 1 or i * S
                r1Small->push_back(element);
            }

            ySize++;
        }
        r1File.close();
        RectanglesD1.clear();

        multiBoxWithOrderIndex RectanglesD0Second = LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d1.tmp");
        centerOrdering(RectanglesD0Second, 0);

        long long currentX = 0;
        std::ofstream r0FileSecond = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
        std::shared_ptr<multiBoxWithOrderIndex> r0Small = std::make_shared<multiBoxWithOrderIndex>();
        r0Small->push_back(RectanglesD0Second[0]);
        r0Small->push_back(RectanglesD0Second[RectanglesD0Second.size() - 1]);
        for (rTreeValueWithOrderIndex element : RectanglesD0Second) {
            Rtree::SaveEntryWithOrderIndex(element, r0FileSecond);

            if (((currentX + 1) % currentS == 0 && (currentX + 1) / currentS >= 1 && (currentX + 1) / currentS < iLimit)
                || (currentX % currentS == 0 && currentX / currentS >= 1 && currentX / currentS < iLimit)) {
                // index i * S - 1 or i * S
                r0Small->push_back(element);
            }

            currentX++;
        }
        r0FileSecond.close();
        RectanglesD0Second.clear();
        // END

        orderedInputRectangles.CreateOrderedBoxesOnDisk(onDiskBase + ".boundingbox.d0", onDiskBase + ".boundingbox.d1", r0Small, r1Small, xSize,
                                                        createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY));
        std::cout << "Finished initial sorting" << std::endl;
    }

    // build the tree in a depth first approach
    std::stack<ConstructionNode> layerStack;

    long long newId = 1; // start from 1, because 0 is the root item
    ConstructionNode rootItem = ConstructionNode(0, orderedInputRectangles);
    layerStack.push(rootItem);
    size_t layer = 0;

    while (!layerStack.empty()) {
        ConstructionNode currentItem = layerStack.top();
        layerStack.pop();

        if (currentItem.GetOrderedBoxes().GetSize() <= M) {
            // reached a leaf
            currentItem.AddChildrenToItem();
            long long nodePtr = SaveNode(currentItem, true, nodesOfs);
            lookup[currentItem.GetId()] = nodePtr;
        } else {
            std::vector<OrderedBoxes> tgsResult = TGSRecursive(onDiskBase + ".boundingbox." + std::to_string(layer), currentItem.GetOrderedBoxes(), M, std::ceil(((float) currentItem.GetOrderedBoxes().GetSize()) / ((float) M)), this->maxBuildingRamUsage);
            if (tgsResult.size() > M) {
                std::cout << "Error while building the Rtree: Number of children (" << tgsResult.size() << ") is not equal to M (" << M << ")" << std::endl;
            }
            for (OrderedBoxes& currentOrderedRectangles : tgsResult) {
                ConstructionNode newItem = ConstructionNode(newId, currentOrderedRectangles);
                layerStack.push(newItem);

                currentItem.AddChild(newItem);

                newId++;
            }

            long long nodePtr = SaveNode(currentItem, false, nodesOfs);
            lookup[currentItem.GetId()] = nodePtr;
        }
        layer++;
    }
    nodesOfs.close();

    std::ofstream lookupOfs(folder + "/lookup.bin", std::ios::binary);
    for (unsigned int i = 0; i < newId; i++) {
        long long nodePtr = lookup[i];
        lookupOfs.write(reinterpret_cast<const char *>(&nodePtr), sizeof(long long));
    }
    lookupOfs.close();
}

multiBoxGeo Rtree::SearchTree(boxGeo query, const std::string &folder) {
    std::ifstream lookupIfs = std::ifstream(folder + "/lookup.bin", std::ios::binary);
    std::ifstream nodesIfs = std::ifstream(folder + "/nodes.bin", std::ios::binary);

    Node rootNode = LoadNode(0, lookupIfs, nodesIfs);
    multiBoxGeo results;
    std::stack<Node> nodes;
    nodes.push(rootNode);

    while(!nodes.empty()) {
        Node currentNode = nodes.top();
        nodes.pop();

        for (rTreeValue child : currentNode.GetChildren()) {
            if (intersects(query, child.first)) {
                if (currentNode.GetIsLastInnerNode()) {
                    results.push_back(child);
                } else {
                    Node newNode = LoadNode(child.second, lookupIfs, nodesIfs);
                    nodes.push(newNode);
                }
            }
        }
    }

    lookupIfs.close();
    nodesIfs.close();
    return results;
}

ConstructionNode::ConstructionNode(long long id, OrderedBoxes orderedBoxes)
        : Node{id}
{
    this->orderedBoxes = orderedBoxes;

    // calculate the boundingBoxes
    this->boundingBox = orderedBoxes.GetBoundingBox();
}

void ConstructionNode::AddChildrenToItem() {
    if (this->GetOrderedBoxes().WorkInRam()) {
        for(rTreeValueWithOrderIndex box : *this->GetOrderedBoxes().GetRectanglesD0InRam()) {
            Node leafNode = Node(box.first.second, box.first.first);
            this->AddChild(leafNode);
        }
    } else {
        std::ifstream loadEntriesIfs = std::ifstream(this->GetOrderedBoxes().GetRectanglesD0OnDisk(), std::ios::binary);
        loadEntriesIfs.seekg (0, std::ifstream::end);
        long long fileLength = loadEntriesIfs.tellg();
        loadEntriesIfs.seekg (0, std::ifstream::beg);

        double minX;
        double minY;
        double maxX;
        double maxY;
        uint64_t id;
        long long orderX;
        long long orderY;

        while (loadEntriesIfs.tellg() < fileLength) {
            loadEntriesIfs.read(reinterpret_cast<char*>(&minX), sizeof(double));
            loadEntriesIfs.read(reinterpret_cast<char*>(&minY), sizeof(double));
            loadEntriesIfs.read(reinterpret_cast<char*>(&maxX), sizeof(double));
            loadEntriesIfs.read(reinterpret_cast<char*>(&maxY), sizeof(double));
            loadEntriesIfs.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
            loadEntriesIfs.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
            loadEntriesIfs.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

            boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
            rTreeValue boxWithId = std::make_pair(box, id);
            Node leafNode = Node(boxWithId.second, boxWithId.first);
            this->AddChild(leafNode);
        }

        loadEntriesIfs.close();
    }
}

long long Node::GetId() const {
    return this->id;
}

OrderedBoxes ConstructionNode::GetOrderedBoxes() {
    return this->orderedBoxes;
}

Node::Node(long long id, boxGeo boundingBox) {
    this->id = id;
    this->boundingBox = boundingBox;
}

Node::Node(long long id) {
    this->id = id;
}

Node::Node() {}

Node::Node(long long id, boxGeo boundingBox, multiBoxGeo &children, bool isLastInnerNode) {
    this->id = id;
    this->boundingBox = boundingBox;
    this->children = children;
    this->isLastInnerNode = isLastInnerNode;
}

Node::Node(long long id, double minX, double minY, double maxX, double maxY, bool isLastInnerNode) {
    this->id = id;
    this->boundingBox = Rtree::createBoundingBox(minX, minY, maxX, maxY);
    this->isLastInnerNode = isLastInnerNode;
}

void Node::AddChild(Node& child) {
    boxGeo box = child.GetBoundingBox();
    unsigned long long entryId = child.GetId();
    rTreeValue entry = std::make_pair(box, entryId);
    this->children.push_back(entry);
}

boxGeo Node::GetBoundingBox() const {
    return this->boundingBox;
}

void Node::SetIsLastInnerNode(bool _isLastInnerNode) {
    this->isLastInnerNode = _isLastInnerNode;
}

bool Node::GetIsLastInnerNode() {
    return this->isLastInnerNode;
}

multiBoxGeo Node::GetChildren() {
    return this->children;
}

long long Rtree::SaveNode(Node &node, bool isLastInnerNode, std::ofstream& nodesOfs) {
    node.SetIsLastInnerNode(isLastInnerNode);

    long long pos = static_cast<long long>(nodesOfs.tellp());
    boost::archive::binary_oarchive archive(nodesOfs);
    archive << node;
    nodesOfs.write(" ", 1);

    return pos;
}

Node Rtree::LoadNode(long long id, std::ifstream& lookupIfs, std::ifstream& nodesIfs) {
    Node newNode;

    long long offset = id * (long long)sizeof(long long);
    lookupIfs.seekg(offset, std::ios::beg);

    long long nodePtr;
    lookupIfs.read(reinterpret_cast<char*>(&nodePtr), sizeof(long long));

    nodesIfs.seekg(nodePtr);
    boost::archive::binary_iarchive ia(nodesIfs);
    ia >> newNode;

    return newNode;
}

std::optional<boxGeo> GetBoundingBoxFromWKT(const std::string& wkt) {
    bool lookingForX = true;
    bool readingDouble = false;
    std::string currentDouble;

    double minX = -1;
    double maxX = -1;
    double minY = -1;
    double maxY = -1;

    for (char c : wkt) {
        if (isdigit(c)) {
            readingDouble = true;
            currentDouble += c;
        } else if (c == '.') {
            readingDouble = true;
            currentDouble += '.';
        } else if (c == ' ') {
            if (readingDouble && lookingForX) {
                // x is completely read in
                readingDouble = false;
                lookingForX = false;
                double x;
                try {
                    x = std::stod(currentDouble);
                } catch(...) {
                    return { };
                }
                currentDouble = "";
                if (x < minX || minX == -1) {
                    minX = x;
                }

                if (x > maxX) {
                    maxX = x;
                }
            }
        } else {
            if (readingDouble && !lookingForX) {
                // y is completely read in
                readingDouble = false;
                lookingForX = true;
                double y;
                try {
                    y = std::stod(currentDouble);
                } catch(...) {
                    return { };
                }
                currentDouble = "";
                if (y < minY || minY == -1) {
                    minY = y;
                }

                if (y > maxY) {
                    maxY = y;
                }
            }
        }
    }

    return { Rtree::createBoundingBox(minX, minY, maxX, maxY) };
}

std::optional<boxGeo> Rtree::ConvertWordToRtreeEntry(const std::string& wkt) {
    std::optional<boxGeo> boundingBox;

    /* Get the bounding box(es) of either a multipolygon, polygon or a linestring */
    std::size_t posWKTStart = wkt.find("MULTIPOLYGON(((") + 14;
    std::size_t posWKTEnd = wkt.find(")))", posWKTStart);
    if (posWKTStart != std::string::npos && posWKTEnd != std::string::npos) {
        std::string newWkt = wkt.substr(posWKTStart, posWKTEnd - posWKTStart + 1);
        boundingBox = GetBoundingBoxFromWKT(newWkt);
    } else {
        posWKTStart = wkt.find("POLYGON((") + 8;
        posWKTEnd = wkt.find("))", posWKTStart);
        if (posWKTStart != std::string::npos && posWKTEnd != std::string::npos) {
            std::string newWkt = wkt.substr(posWKTStart, posWKTEnd - posWKTStart + 1);
            boundingBox = GetBoundingBoxFromWKT(newWkt);
        } else {
            posWKTStart = wkt.find("LINESTRING(") + 10;
            posWKTEnd = wkt.find(')', posWKTStart);
            if (posWKTStart != std::string::npos && posWKTEnd != std::string::npos) {
                std::string newWkt = wkt.substr(posWKTStart, posWKTEnd - posWKTStart + 1);
                boundingBox = GetBoundingBoxFromWKT(newWkt);
            } else {
                return { };
            }
        }
    }

    return boundingBox;
}

void Rtree::SaveEntry(boxGeo boundingBox, uint64_t index, std::ofstream& convertOfs) {
    /* write the boundingbox value to file */
    double minX = boundingBox.min_corner().get<0>();
    double minY = boundingBox.min_corner().get<1>();
    double maxX = boundingBox.max_corner().get<0>();
    double maxY = boundingBox.max_corner().get<1>();

    convertOfs.write(reinterpret_cast<const char *>(&minX), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&minY), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&index), sizeof(uint64_t));
}

void Rtree::SaveEntryWithOrderIndex(rTreeValueWithOrderIndex treeValue, std::ofstream& convertOfs) {
    /* write the boundingbox value to file */
    double minX = treeValue.first.first.min_corner().get<0>();
    double minY = treeValue.first.first.min_corner().get<1>();
    double maxX = treeValue.first.first.max_corner().get<0>();
    double maxY = treeValue.first.first.max_corner().get<1>();

    convertOfs.write(reinterpret_cast<const char *>(&minX), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&minY), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
    convertOfs.write(reinterpret_cast<const char *>(&treeValue.first.second), sizeof(uint64_t));
    convertOfs.write(reinterpret_cast<const char *>(&treeValue.second.first), sizeof(long long));
    convertOfs.write(reinterpret_cast<const char *>(&treeValue.second.second), sizeof(long long));
}

multiBoxGeo Rtree::LoadEntries(const std::string& file) {
    multiBoxGeo boxes;

    std::ifstream loadEntriesIfs = std::ifstream(file, std::ios::binary);
    loadEntriesIfs.seekg (0, std::ifstream::end);
    long long fileLength = loadEntriesIfs.tellg();
    loadEntriesIfs.seekg (0, std::ifstream::beg);

    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t id;

    while (loadEntriesIfs.tellg() < fileLength) {
        loadEntriesIfs.read(reinterpret_cast<char*>(&minX), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&minY), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

        boxGeo box = createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);
        boxes.push_back(boxWithId);
    }

    loadEntriesIfs.close();
    return boxes;
}

multiBoxWithOrderIndex Rtree::LoadEntriesWithOrderIndex(const std::string& file) {
    multiBoxWithOrderIndex boxes;

    std::ifstream loadEntriesIfs = std::ifstream(file, std::ios::binary);
    loadEntriesIfs.seekg (0, std::ifstream::end);
    long long fileLength = loadEntriesIfs.tellg();
    loadEntriesIfs.seekg (0, std::ifstream::beg);

    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t id;
    long long orderX;
    long long orderY;

    while (loadEntriesIfs.tellg() < fileLength) {
        loadEntriesIfs.read(reinterpret_cast<char*>(&minX), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&minY), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        loadEntriesIfs.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        loadEntriesIfs.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        loadEntriesIfs.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);
        boxes.emplace_back(boxWithId, std::make_pair(orderX, orderY));
    }

    loadEntriesIfs.close();
    return boxes;
}

Rtree::Rtree(uintmax_t maxBuildingRamUsage) {
    this->maxBuildingRamUsage = maxBuildingRamUsage;
}

bool OrderedBoxes::WorkInRam() const{
    return this->workInRam;
}

void OrderedBoxes::CreateOrderedBoxesInRam(const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesD1, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD1, boxGeo box) {
    this->workInRam = true;
    this->rectanglesD0InRam = rectanglesD0;
    this->rectanglesD1InRam = rectanglesD1;
    this->rectanglesD0Small = rectanglesSmallD0;
    this->rectanglesD1Small = rectanglesSmallD1;
    this->size = (*rectanglesD0).size();
    this->boundingBox = box;
}

void OrderedBoxes::CreateOrderedBoxesOnDisk(const std::string& rectanglesD0, const std::string& rectanglesD1, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD0, const std::shared_ptr<multiBoxWithOrderIndex>& rectanglesSmallD1, long long size, boxGeo box) {
    this->workInRam = false;
    this->rectanglesD0OnDisk = rectanglesD0 + ".tmp";
    this->rectanglesD1OnDisk = rectanglesD1 + ".tmp";
    this->rectanglesD0Small = rectanglesSmallD0;
    this->rectanglesD1Small = rectanglesSmallD1;
    this->size = size;
    this->boundingBox = box;
}

boxGeo OrderedBoxes::GetBoundingBox() {
    return this->boundingBox;
}

long long OrderedBoxes::GetSize() const {
    return this->size;
}

std::shared_ptr<multiBoxWithOrderIndex> OrderedBoxes::GetRectanglesD0InRam() {
    return this->rectanglesD0InRam;
}

std::shared_ptr<multiBoxWithOrderIndex> OrderedBoxes::GetRectanglesD1InRam() {
    return this->rectanglesD1InRam;
}

std::string OrderedBoxes::GetRectanglesD0OnDisk() {
    return this->rectanglesD0OnDisk;
}

std::string OrderedBoxes::GetRectanglesD1OnDisk() {
    return this->rectanglesD1OnDisk;
}

SplitResult OrderedBoxes::GetBestSplit() {
    struct SplitResult splitResult;

    rTreeValueWithOrderIndex minElement;
    rTreeValueWithOrderIndex maxElement;
    rTreeValueWithOrderIndex currentLastElement;
    rTreeValueWithOrderIndex currentElement;

    bool currentlyAtSTimesI = false;

    // dim 0

    for (long long i = 0; i < (*this->rectanglesD0Small).size(); i++) {
        currentElement = (*this->rectanglesD0Small)[i];

        if (i == 0) {
            // this is the min element
            minElement = currentElement;
            continue;
        }

        if (i == 1) {
            // this is the max element
            maxElement = currentElement;
            continue;
        }

        double minXB0 = 0;
        double maxXB0 = 1;
        double minXB1 = 0;
        double maxXB1 = 1;

        if (!currentlyAtSTimesI) {
            currentLastElement = currentElement;
            currentlyAtSTimesI = true;
            continue;
        }

        if (currentlyAtSTimesI && currentElement.first.second != maxElement.first.second) {

            minXB0 = (minElement.first.first.min_corner().get<0>() + minElement.first.first.max_corner().get<0>()) / 2;
            maxXB0 = (currentLastElement.first.first.min_corner().get<0>() + currentLastElement.first.first.max_corner().get<0>()) / 2;

            minXB1 = (currentElement.first.first.min_corner().get<0>() + currentElement.first.first.max_corner().get<0>()) / 2;
            maxXB1 = (maxElement.first.first.min_corner().get<0>() + maxElement.first.first.max_corner().get<0>()) / 2;

            currentlyAtSTimesI = false;
        } else {
            break;
        }

        boxGeo b0 = Rtree::createBoundingBox(minXB0, 0, maxXB0, 1);
        boxGeo b1 = Rtree::createBoundingBox(minXB1, 0, maxXB1, 1);


        double cost = costFunctionTGS(b0, b1);

        if (splitResult.bestCost == -1 || cost < splitResult.bestCost) {
            splitResult.bestCost = cost;
            splitResult.bestDim = 0;
            splitResult.bestB0 = b0;
            splitResult.bestLastElement = currentLastElement;
            splitResult.bestElement = currentElement;
            splitResult.bestMinElement = minElement;
            splitResult.bestMaxElement = maxElement;
            splitResult.bestIndex = i;
        }
    }

    // dim 1

    currentlyAtSTimesI = false;

    for (long long i = 0; i < (*this->rectanglesD1Small).size(); i++) {
        currentElement = (*this->rectanglesD1Small)[i];

        if (i == 0) {
            // this is the min element
            minElement = currentElement;
            continue;
        }

        if (i == 1) {
            // this is the max element
            maxElement = currentElement;
            continue;
        }

        double minYB0 = 0;
        double maxYB0 = 1;
        double minYB1 = 0;
        double maxYB1 = 1;

        if (!currentlyAtSTimesI) {
            currentLastElement = currentElement;
            currentlyAtSTimesI = true;
            continue;
        }

        if (currentlyAtSTimesI && currentElement.first.second != maxElement.first.second) {
            minYB0 = (minElement.first.first.min_corner().get<1>() + minElement.first.first.max_corner().get<1>()) / 2;
            maxYB0 = (currentLastElement.first.first.min_corner().get<1>() + currentLastElement.first.first.max_corner().get<1>()) / 2;

            minYB1 = (currentElement.first.first.min_corner().get<1>() + currentElement.first.first.max_corner().get<1>()) / 2;
            maxYB1 = (maxElement.first.first.min_corner().get<1>() + maxElement.first.first.max_corner().get<1>()) / 2;

            currentlyAtSTimesI = false;
        } else {
            break;
        }

        boxGeo b0 = Rtree::createBoundingBox(0, minYB0, 1, maxYB0);
        boxGeo b1 = Rtree::createBoundingBox(0, minYB1, 1, maxYB1);


        double cost = costFunctionTGS(b0, b1);

        if (splitResult.bestCost == -1 || cost < splitResult.bestCost) {
            splitResult.bestCost = cost;
            splitResult.bestDim = 1;
            splitResult.bestB0 = b0;
            splitResult.bestLastElement = currentLastElement;
            splitResult.bestElement = currentElement;
            splitResult.bestMinElement = minElement;
            splitResult.bestMaxElement = maxElement;
            splitResult.bestIndex = i;
        }
    }

    return splitResult;
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBest(const std::string& filePath, size_t S, size_t M, long long maxBuildingRamUsage) {
    if (this->workInRam) {
        return this->SplitAtBestInRam(S, M);
    } else {
        return this->SplitAtBestOnDisk(filePath, S, M, maxBuildingRamUsage);
    }
}

bool CheckIfListsContainSameElements(std::shared_ptr<multiBoxWithOrderIndex> list1, std::shared_ptr<multiBoxWithOrderIndex> list2) {
    if (list1->size() != list2->size()) return false;

    std::map m = std::map<long long, std::pair<long long, long long>>();

    for (rTreeValueWithOrderIndex element : *list1) {
        m[element.first.second] = element.second;

        if (element.first.second == 5199) {
            int i = 0;
        }
    }

    if (list1->size() != m.size()) return false;

    for (rTreeValueWithOrderIndex element : *list2) {
        if (element.first.second == 5199) {
            int i = 0;
        }
        if (m.find(element.first.second) == m.end()) return false;
        long long id = element.first.second;
        std::pair<long long, long long> pair1 = m[element.first.second];
        std::pair<long long, long long> pair2 = element.second;
        if (m[element.first.second] != element.second) return false;
        m.erase(element.first.second);
    }

    return true;
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBestInRam(size_t S, size_t M) {
    struct SplitResult splitResult = this->GetBestSplit();

    OrderedBoxes split0;
    OrderedBoxes split1;

    std::shared_ptr<multiBoxWithOrderIndex> s0Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s0Dim1 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1Dim1 = std::make_shared<multiBoxWithOrderIndex>();

    std::shared_ptr<multiBoxWithOrderIndex> s0SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s0SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> s1SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();

    long long sizeLeft = std::ceil((splitResult.bestIndex - 2) / 2.0) * S;
    long long sizeRight = this->size - sizeLeft;
    size_t SSplit0 = sizeLeft <= S ? std::ceil(sizeLeft / (double) M) : S;
    size_t SSplit1 = sizeRight <= S ? std::ceil(sizeRight / (double) M) : S;
    /*SSplit0 = 6251;
    SSplit1 = 6251;*/

    /*if (sizeLeft % (int)std::ceil(this->size / (float) M) != 0 || sizeRight % (int)std::ceil(this->size / (float) M) != 0) {
        int fdg = 0;
    }*/
    //long long iLimit = std::ceil(this->size / S);

    double globalMinXS0 = -1;
    double globalMinYS0 = -1;
    double globalMaxXS0 = -1;
    double globalMaxYS0 = -1;

    double globalMinXS1 = -1;
    double globalMinYS1 = -1;
    double globalMaxXS1 = -1;
    double globalMaxYS1 = -1;

    // placeholder
    if (splitResult.bestDim == 0) {
        s0SmallDim0->push_back((*this->rectanglesD0InRam)[0]);
        s0SmallDim0->push_back(splitResult.bestLastElement);
        s1SmallDim0->push_back(splitResult.bestElement);
        s1SmallDim0->push_back((*this->rectanglesD0InRam)[this->rectanglesD0InRam->size() - 1]);
        s0SmallDim1->push_back(rTreeValueWithOrderIndex());
        s0SmallDim1->push_back(rTreeValueWithOrderIndex());
        s1SmallDim1->push_back(rTreeValueWithOrderIndex());
        s1SmallDim1->push_back(rTreeValueWithOrderIndex());
    } else {
        s0SmallDim1->push_back((*this->rectanglesD1InRam)[0]);
        s0SmallDim1->push_back(splitResult.bestLastElement);
        s1SmallDim1->push_back(splitResult.bestElement);
        s1SmallDim1->push_back((*this->rectanglesD1InRam)[this->rectanglesD1InRam->size() - 1]);
        s0SmallDim0->push_back(rTreeValueWithOrderIndex());
        s0SmallDim0->push_back(rTreeValueWithOrderIndex());
        s1SmallDim0->push_back(rTreeValueWithOrderIndex());
        s1SmallDim0->push_back(rTreeValueWithOrderIndex());
    }

    long long currentXSplit0 = 0;
    long long currentXSplit1 = 0;
    for (rTreeValueWithOrderIndex element : *this->rectanglesD0InRam) {
        if ((splitResult.bestDim == 0 && element.second.first < splitResult.bestElement.second.first)
            || (splitResult.bestDim == 1 && element.second.second < splitResult.bestElement.second.second)) {
            // split 0

            s0Dim0->push_back(element);
            if (((currentXSplit0 + 1) % SSplit0 == 0 && (currentXSplit0 + 1) / SSplit0 >= 1 && (currentXSplit0 + 1) / SSplit0 < M)
                || (currentXSplit0 % SSplit0 == 0 && currentXSplit0 / SSplit0 >= 1 && currentXSplit0 / SSplit0 < M)) {
                // index i * S - 1 or i * S
                s0SmallDim0->push_back(element);
            }

            if (globalMinXS0 == -1 || element.first.first.min_corner().get<0>() < globalMinXS0) {
                globalMinXS0 = element.first.first.min_corner().get<0>();
            }
            if (globalMinYS0 == -1 || element.first.first.min_corner().get<1>() < globalMinYS0) {
                globalMinYS0 = element.first.first.min_corner().get<1>();
            }
            if (element.first.first.max_corner().get<0>() > globalMaxXS0) {
                globalMaxXS0 = element.first.first.max_corner().get<0>();
            }
            if (element.first.first.max_corner().get<1>() > globalMaxYS0) {
                globalMaxYS0 = element.first.first.max_corner().get<1>();
            }

            currentXSplit0++;
        } else {
            // split 1

            s1Dim0->push_back(element);
            if (((currentXSplit1 + 1) % SSplit1 == 0 && (currentXSplit1 + 1) / SSplit1 >= 1 && (currentXSplit1 + 1) / SSplit1 < M)
                || (currentXSplit1 % SSplit1 == 0 && currentXSplit1 / SSplit1 >= 1 && currentXSplit1 / SSplit1 < M)) {
                // index i * S - 1 or i * S
                s1SmallDim0->push_back(element);
            }

            if (globalMinXS1 == -1 || element.first.first.min_corner().get<0>() < globalMinXS1) {
                globalMinXS1 = element.first.first.min_corner().get<0>();
            }
            if (globalMinYS1 == -1 || element.first.first.min_corner().get<1>() < globalMinYS1) {
                globalMinYS1 = element.first.first.min_corner().get<1>();
            }
            if (element.first.first.max_corner().get<0>() > globalMaxXS1) {
                globalMaxXS1 = element.first.first.max_corner().get<0>();
            }
            if (element.first.first.max_corner().get<1>() > globalMaxYS1) {
                globalMaxYS1 = element.first.first.max_corner().get<1>();
            }

            currentXSplit1++;
        }
    }

    long long currentYSplit0 = 0;
    long long currentYSplit1 = 0;
    for (rTreeValueWithOrderIndex element : *this->rectanglesD1InRam) {
        if ((splitResult.bestDim == 0 && element.second.first < splitResult.bestElement.second.first)
            || (splitResult.bestDim == 1 && element.second.second < splitResult.bestElement.second.second)) {
            // split 0

            s0Dim1->push_back(element);
            if (((currentYSplit0 + 1) % SSplit0 == 0 && (currentYSplit0 + 1) / SSplit0 >= 1 && (currentYSplit0 + 1) / SSplit0 < M)
                || (currentYSplit0 % SSplit0 == 0 && currentYSplit0 / SSplit0 >= 1 && currentYSplit0 / SSplit0 < M)) {
                // index i * S - 1 or i * S
                s0SmallDim1->push_back(element);
            }

            currentYSplit0++;
        } else {
            // split 1

            s1Dim1->push_back(element);
            if (((currentYSplit1 + 1) % SSplit1 == 0 && (currentYSplit1 + 1) / SSplit1 >= 1 && (currentYSplit1 + 1) / SSplit1 < M)
                || (currentYSplit1 % SSplit1 == 0 && currentYSplit1 / SSplit1 >= 1 && currentYSplit1 / SSplit1 < M)) {
                // index i * S - 1 or i * S
                s1SmallDim1->push_back(element);
            }

            currentYSplit1++;
        }
    }

    if (s0Dim0->size() != s0Dim1->size() || s1Dim0->size() != s1Dim1->size() || s0Dim0->size() == 0 || s1Dim0->size() == 0)
    {
        long long a = s0Dim0->size();
        long long b = s0Dim1->size();
        long long c = s1Dim0->size();
        long long d = s1Dim1->size();
        long long e = this->rectanglesD0InRam->size();
        long long f = this->rectanglesD1InRam->size();
        bool test = true;
    }
    // replace the placeholder
    if (splitResult.bestDim == 0) {
        (*s0SmallDim1)[0] = (*s0Dim1)[0];
        (*s0SmallDim1)[1] = (*s0Dim1)[s0Dim1->size() - 1];
        (*s1SmallDim1)[0] = (*s1Dim1)[0];
        (*s1SmallDim1)[1] = (*s1Dim1)[s1Dim1->size() - 1];
    } else {
        (*s0SmallDim0)[0] = (*s0Dim0)[0];
        (*s0SmallDim0)[1] = (*s0Dim0)[s0Dim0->size() - 1];
        (*s1SmallDim0)[0] = (*s1Dim0)[0];
        (*s1SmallDim0)[1] = (*s1Dim0)[s1Dim0->size() - 1];
    }

    if (s0SmallDim0->size() <= 2 || s0SmallDim1->size() <= 2 || s1SmallDim0->size() <= 2 || s1SmallDim1->size() <= 2) {
        long long a = s0SmallDim0->size();
        long long b = s0SmallDim1->size();
        long long c = s1SmallDim0->size();
        long long d = s1SmallDim1->size();
        long long e = this->rectanglesD0InRam->size();
        long long f = this->rectanglesD1InRam->size();
        bool test = true;
    }

    split0.CreateOrderedBoxesInRam(s0Dim0, s0Dim1, s0SmallDim0, s0SmallDim1, Rtree::createBoundingBox(globalMinXS0, globalMinYS0, globalMaxXS0, globalMaxYS0));
    split1.CreateOrderedBoxesInRam(s1Dim0, s1Dim1, s1SmallDim0, s1SmallDim1, Rtree::createBoundingBox(globalMinXS1, globalMinYS1, globalMaxXS1, globalMaxYS1));

    (*this->rectanglesD0InRam).clear();
    (*this->rectanglesD1InRam).clear();
    (*this->rectanglesD0Small).clear();
    (*this->rectanglesD1Small).clear();
    (*this->rectanglesD0InRam).shrink_to_fit();
    (*this->rectanglesD1InRam).shrink_to_fit();
    (*this->rectanglesD0Small).shrink_to_fit();
    (*this->rectanglesD1Small).shrink_to_fit();

    return std::make_pair(split0, split1);
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBestOnDisk(const std::string& filePath, size_t S, size_t M, long long maxBuildingRamUsage) {
    OrderedBoxes split0;
    OrderedBoxes split1;

    struct SplitResult splitResult = this->GetBestSplit();

    // perfrom the split
    long long sizeLeft = std::ceil((splitResult.bestIndex - 2) / 2.0) * S;
    long long sizeRight = this->size - sizeLeft;
    size_t SSplit0 = sizeLeft <= S ? std::ceil(sizeLeft / (double) M) : S;
    size_t SSplit1 = sizeRight <= S ? std::ceil(sizeRight / (double) M) : S;
    long long split0ByteSize = sizeLeft * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    long long split1ByteSize = sizeRight * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    bool split0InRam = split0ByteSize * 4 < maxBuildingRamUsage;
    bool split1InRam = split1ByteSize * 4 < maxBuildingRamUsage;

    std::shared_ptr<multiBoxWithOrderIndex> split0Dim0InRam = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> split0Dim1InRam = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> split1Dim0InRam = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> split1Dim1InRam = std::make_shared<multiBoxWithOrderIndex>();
    std::ofstream split0Dim0File;
    std::ofstream split0Dim1File;
    std::shared_ptr<multiBoxWithOrderIndex> split0Dim0Small = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> split0Dim1Small = std::make_shared<multiBoxWithOrderIndex>();
    std::ofstream split1Dim0File;
    std::ofstream split1Dim1File;
    std::shared_ptr<multiBoxWithOrderIndex> split1Dim0Small = std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> split1Dim1Small = std::make_shared<multiBoxWithOrderIndex>();

    if (!split0InRam) {
        split0Dim0File = std::ofstream(filePath + ".0.dim0.tmp", std::ios::binary);
        split0Dim1File = std::ofstream(filePath + ".0.dim1.tmp", std::ios::binary);
    }

    if (splitResult.bestDim == 0) {
        split0Dim0Small->push_back(splitResult.bestMinElement);
        split0Dim0Small->push_back(splitResult.bestLastElement);
        split1Dim0Small->push_back(splitResult.bestElement);
        split1Dim0Small->push_back(splitResult.bestMaxElement);

        // placeholder
        split0Dim1Small->push_back(rTreeValueWithOrderIndex());
        split0Dim1Small->push_back(rTreeValueWithOrderIndex());
        split1Dim1Small->push_back(rTreeValueWithOrderIndex());
        split1Dim1Small->push_back(rTreeValueWithOrderIndex());
    } else {
        split0Dim1Small->push_back(splitResult.bestMinElement);
        split0Dim1Small->push_back(splitResult.bestLastElement);
        split1Dim1Small->push_back(splitResult.bestElement);
        split1Dim1Small->push_back(splitResult.bestMaxElement);

        // placeholder
        split0Dim0Small->push_back(rTreeValueWithOrderIndex());
        split0Dim0Small->push_back(rTreeValueWithOrderIndex());
        split1Dim0Small->push_back(rTreeValueWithOrderIndex());
        split1Dim0Small->push_back(rTreeValueWithOrderIndex());
    }

    if (!split1InRam) {
        split1Dim0File = std::ofstream(filePath + ".1.dim0.tmp", std::ios::binary);
        split1Dim1File = std::ofstream(filePath + ".1.dim1.tmp", std::ios::binary);
    }

    rTreeValueWithOrderIndex minSplit0OtherDim;
    rTreeValueWithOrderIndex maxSplit0OtherDim;
    rTreeValueWithOrderIndex minSplit1OtherDim;
    rTreeValueWithOrderIndex maxSplit1OtherDim;

    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t id;
    long long orderX;
    long long orderY;

    double globalMinXS0 = -1;
    double globalMinYS0 = -1;
    double globalMaxXS0 = -1;
    double globalMaxYS0 = -1;

    double globalMinXS1 = -1;
    double globalMinYS1 = -1;
    double globalMaxXS1 = -1;
    double globalMaxYS1 = -1;

    std::ifstream dim0WholeFile = std::ifstream(this->rectanglesD0OnDisk, std::ios::binary);
    dim0WholeFile.seekg (0, std::ifstream::end);
    long long fileLength = dim0WholeFile.tellg();
    dim0WholeFile.seekg (0, std::ifstream::beg);

    long long currentXSplit0 = 0;
    long long currentXSplit1 = 0;
    while (dim0WholeFile.tellg() < fileLength) {
        dim0WholeFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
        dim0WholeFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
        dim0WholeFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        dim0WholeFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        dim0WholeFile.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        dim0WholeFile.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        dim0WholeFile.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);
        rTreeValueWithOrderIndex currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

        if ((splitResult.bestDim == 0 && currentElement.second.first < splitResult.bestElement.second.first)
        || (splitResult.bestDim == 1 && currentElement.second.second < splitResult.bestElement.second.second)) {
            // it's in split 0
            if (((currentXSplit0 + 1) % SSplit0 == 0 && (currentXSplit0 + 1) / SSplit0 >= 1 && (currentXSplit0 + 1) / SSplit0 < M)
                || (currentXSplit0 % SSplit0 == 0 && currentXSplit0 / SSplit0 >= 1 && currentXSplit0 / SSplit0 < M)) {
                // index i * S - 1 or i * S
                split0Dim0Small->push_back(currentElement);
            }

            if (split0InRam) {
                split0Dim0InRam->push_back(currentElement);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim0File);
            }

            if (globalMinXS0 == -1 || currentElement.first.first.min_corner().get<0>() < globalMinXS0) {
                globalMinXS0 = currentElement.first.first.min_corner().get<0>();
            }
            if (globalMinYS0 == -1 || currentElement.first.first.min_corner().get<1>() < globalMinYS0) {
                globalMinYS0 = currentElement.first.first.min_corner().get<1>();
            }
            if (currentElement.first.first.max_corner().get<0>() > globalMaxXS0) {
                globalMaxXS0 = currentElement.first.first.max_corner().get<0>();
            }
            if (currentElement.first.first.max_corner().get<1>() > globalMaxYS0) {
                globalMaxYS0 = currentElement.first.first.max_corner().get<1>();
            }

            if (splitResult.bestDim == 1) {
                if (currentXSplit0 == 0) {
                    minSplit0OtherDim = currentElement;
                    maxSplit0OtherDim = currentElement;
                }
                if (currentElement.second.first > maxSplit0OtherDim.second.first) {
                    maxSplit0OtherDim = currentElement;
                }
            }

            currentXSplit0++;
        } else {
            // it's in split 1
            if (((currentXSplit1 + 1) % SSplit1 == 0 && (currentXSplit1 + 1) / SSplit1 >= 1 && (currentXSplit1 + 1) / SSplit1 < M)
                || (currentXSplit1 % SSplit1 == 0 && currentXSplit1 / SSplit1 >= 1 && currentXSplit1 / SSplit1 < M)) {
                // index i * S - 1 or i * S
                split1Dim0Small->push_back(currentElement);
            }

            if (split1InRam) {
                split1Dim0InRam->push_back(currentElement);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim0File);
            }

            if (globalMinXS1 == -1 || currentElement.first.first.min_corner().get<0>() < globalMinXS1) {
                globalMinXS1 = currentElement.first.first.min_corner().get<0>();
            }
            if (globalMinYS1 == -1 || currentElement.first.first.min_corner().get<1>() < globalMinYS1) {
                globalMinYS1 = currentElement.first.first.min_corner().get<1>();
            }
            if (currentElement.first.first.max_corner().get<0>() > globalMaxXS1) {
                globalMaxXS1 = currentElement.first.first.max_corner().get<0>();
            }
            if (currentElement.first.first.max_corner().get<1>() > globalMaxYS1) {
                globalMaxYS1 = currentElement.first.first.max_corner().get<1>();
            }

            if (splitResult.bestDim == 1) {
                if (currentXSplit1 == 0) {
                    minSplit1OtherDim = currentElement;
                    maxSplit1OtherDim = currentElement;
                }
                if (currentElement.second.first > maxSplit1OtherDim.second.first) {
                    maxSplit1OtherDim = currentElement;
                }
            }

            currentXSplit1++;
        }
    }
    dim0WholeFile.close();

    std::ifstream dim1WholeFile = std::ifstream(this->rectanglesD1OnDisk, std::ios::binary);
    dim1WholeFile.seekg (0, std::ifstream::end);
    fileLength = dim1WholeFile.tellg();
    dim1WholeFile.seekg (0, std::ifstream::beg);

    long long currentYSplit0 = 0;
    long long currentYSplit1 = 0;
    while (dim1WholeFile.tellg() < fileLength) {
        dim1WholeFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
        dim1WholeFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
        dim1WholeFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        dim1WholeFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        dim1WholeFile.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        dim1WholeFile.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        dim1WholeFile.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);
        rTreeValueWithOrderIndex currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

        if ((splitResult.bestDim == 1 && currentElement.second.second < splitResult.bestElement.second.second)
        || (splitResult.bestDim == 0 && currentElement.second.first < splitResult.bestElement.second.first)) {
            // it's in split 0
            if (((currentYSplit0 + 1) % SSplit0 == 0 && (currentYSplit0 + 1) / SSplit0 >= 1 && (currentYSplit0 + 1) / SSplit0 < M)
                || (currentYSplit0 % SSplit0 == 0 && currentYSplit0 / SSplit0 >= 1 && currentYSplit0 / SSplit0 < M)) {
                // index i * S - 1 or i * S
                split0Dim1Small->push_back(currentElement);
            }

            if (split0InRam) {
                split0Dim1InRam->push_back(currentElement);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim1File);
            }

            if (splitResult.bestDim == 0) {
                if (currentYSplit0 == 0) {
                    minSplit0OtherDim = currentElement;
                    maxSplit0OtherDim = currentElement;
                }
                if (currentElement.second.second > maxSplit0OtherDim.second.second) {
                    maxSplit0OtherDim = currentElement;
                }
            }

            currentYSplit0++;
        } else {
            // it's in split 1
            if (((currentYSplit1 + 1) % SSplit1 == 0 && (currentYSplit1 + 1) / SSplit1 >= 1 && (currentYSplit1 + 1) / SSplit1 < M)
                || (currentYSplit1 % SSplit1 == 0 && currentYSplit1 / SSplit1 >= 1 && currentYSplit1 / SSplit1 < M)) {
                // index i * S - 1 or i * S
                split1Dim1Small->push_back(currentElement);
            }

            if (split1InRam) {
                split1Dim1InRam->push_back(currentElement);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim1File);
            }

            if (splitResult.bestDim == 0) {
                if (currentYSplit1 == 0) {
                    minSplit1OtherDim = currentElement;
                    maxSplit1OtherDim = currentElement;
                }
                if (currentElement.second.second > maxSplit1OtherDim.second.second) {
                    maxSplit1OtherDim = currentElement;
                }
            }

            currentYSplit1++;
        }
    }
    dim1WholeFile.close();

    // replace the placeholder with actual values
    if (splitResult.bestDim == 0) {
        (*split0Dim1Small)[0] = minSplit0OtherDim;
        (*split0Dim1Small)[1] = maxSplit0OtherDim;

        (*split1Dim1Small)[0] = minSplit1OtherDim;
        (*split1Dim1Small)[1] = maxSplit1OtherDim;
    } else {
        (*split0Dim0Small)[0] = minSplit0OtherDim;
        (*split0Dim0Small)[1] = maxSplit0OtherDim;

        (*split1Dim0Small)[0] = minSplit1OtherDim;
        (*split1Dim0Small)[1] = maxSplit1OtherDim;
    }

    boxGeo boxSplit0 = Rtree::createBoundingBox(globalMinXS0, globalMinYS0, globalMaxXS0, globalMaxYS0);
    boxGeo boxSplit1 = Rtree::createBoundingBox(globalMinXS1, globalMinYS1, globalMaxXS1, globalMaxYS1);

    if (!split0InRam) {
        split0Dim0File.close();
        split0Dim1File.close();

        split0.CreateOrderedBoxesOnDisk(filePath + ".0.dim0", filePath + ".0.dim1", split0Dim0Small, split0Dim1Small, sizeLeft, boxSplit0);
    } else {
        split0.CreateOrderedBoxesInRam(split0Dim0InRam, split0Dim1InRam, split0Dim0Small, split0Dim1Small, boxSplit0);
    }

    if (!split1InRam) {
        split1Dim0File.close();
        split1Dim1File.close();

        split1.CreateOrderedBoxesOnDisk(filePath + ".1.dim0", filePath + ".1.dim1", split1Dim0Small, split1Dim1Small, sizeRight, boxSplit1);
    } else {
        split1.CreateOrderedBoxesInRam(split1Dim0InRam, split1Dim1InRam, split1Dim0Small, split1Dim1Small, boxSplit1);
    }

    std::remove(this->rectanglesD0OnDisk.c_str());
    std::remove(this->rectanglesD1OnDisk.c_str());

    return std::make_pair(split0, split1);
}
