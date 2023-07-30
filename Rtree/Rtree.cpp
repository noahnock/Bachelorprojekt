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
            return center1 < center2;
        };

        std::sort(boxes.begin(), boxes.end(), sortRuleLambda);
    } else {
        // order by centerY
        auto sortRuleLambda = [](rTreeValueWithOrderIndex b1, rTreeValueWithOrderIndex b2) -> bool {
            double center1 = (b1.first.first.min_corner().get<1>() + b1.first.first.max_corner().get<1>()) / 2;
            double center2 = (b2.first.first.min_corner().get<1>() + b2.first.first.max_corner().get<1>()) / 2;
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

static bool boxInMultiBox(rTreeValue& item, multiBoxGeo& boxes) {
    for (rTreeValue currentItem : boxes) {
        boxGeo currentBox = currentItem.first;
        boxGeo box = item.first;
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
    //if (orderedInputRectangles[0].size() != orderedInputRectangles[1].size()) {
    //    throw std::length_error("The sizes of the x and y sortings do not match");
    //}

    unsigned long long n = orderedInputRectangles.GetSize();

    if (n <= S || n <= M) {
        // stop condition
        return std::vector<OrderedBoxes> { orderedInputRectangles };
    }

    // split the rectangles at the best split
    std::pair<OrderedBoxes, OrderedBoxes> split = orderedInputRectangles.SplitAtBest(filePath, S, maxBuildingRamUsage);

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
    if (std::filesystem::file_size(file) * 4 < this->maxBuildingRamUsage) {
        // do it in ram
        multiBoxGeo inputRectangles = LoadEntries(file);

        multiBoxGeo RectanglesD0(inputRectangles);
        multiBoxGeo RectanglesD1(inputRectangles);

        centerOrdering(RectanglesD0, 0);
        centerOrdering(RectanglesD1, 1);

        orderedInputRectangles.CreateOrderedBoxesInRam(RectanglesD0, RectanglesD1);
    } else {
        // do it on disk

        // TODO sort on disk instead of ram
        // START
        multiBoxGeo RectanglesD0 = LoadEntries(file);

        centerOrdering(RectanglesD0, 0);

        long long xSize = 0;

        std::ofstream r0File = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
        for (rTreeValue element : RectanglesD0) {
            rTreeValueWithOrderIndex entry = std::make_pair(element, std::make_pair(xSize, 0));
            Rtree::SaveEntryWithOrderIndex(entry, r0File);
            xSize++;
        }
        r0File.close();
        RectanglesD0.clear();

        multiBoxWithOrderIndex RectanglesD1 = LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d0.tmp");
        centerOrdering(RectanglesD1, 1);

        size_t currentS = std::ceil(((float) xSize) / ((float) M));
        size_t iLimit = std::ceil(xSize / currentS);

        long long ySize = 0;
        std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
        std::ofstream r1FileSmall = std::ofstream(onDiskBase + ".boundingbox.d1.small.tmp", std::ios::binary);
        Rtree::SaveEntryWithOrderIndex(RectanglesD1[0], r1FileSmall);
        rTreeValueWithOrderIndex maxElementDim1 = RectanglesD1[RectanglesD1.size() - 1];
        maxElementDim1.second.second = RectanglesD1.size() - 1;
        Rtree::SaveEntryWithOrderIndex(maxElementDim1, r1FileSmall);
        for (rTreeValueWithOrderIndex element : RectanglesD1) {
            element.second.second = ySize;
            Rtree::SaveEntryWithOrderIndex(element, r1File);

            if (((ySize + 1) % currentS == 0 && (ySize + 1) / currentS >= 1 && (ySize + 1) / currentS < iLimit)
            || (ySize % currentS == 0 && ySize / currentS >= 1 && ySize / currentS < iLimit)) {
                // index i * S - 1 or i * S
                Rtree::SaveEntryWithOrderIndex(element, r1FileSmall);
            }

            ySize++;
        }
        r1File.close();
        r1FileSmall.close();
        RectanglesD1.clear();

        multiBoxWithOrderIndex RectanglesD0Second = LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d1.tmp");
        centerOrdering(RectanglesD0Second, 0);

        long long currentX = 0;
        std::ofstream r0FileSecond = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
        std::ofstream r0FileSmall = std::ofstream(onDiskBase + ".boundingbox.d0.small.tmp", std::ios::binary);
        Rtree::SaveEntryWithOrderIndex(RectanglesD0Second[0], r0FileSmall);
        Rtree::SaveEntryWithOrderIndex(RectanglesD0Second[RectanglesD0Second.size() - 1], r0FileSmall);
        for (rTreeValueWithOrderIndex element : RectanglesD0Second) {
            Rtree::SaveEntryWithOrderIndex(element, r0FileSecond);

            if (((currentX + 1) % currentS == 0 && (currentX + 1) / currentS >= 1 && (currentX + 1) / currentS < iLimit)
                || (currentX % currentS == 0 && currentX / currentS >= 1 && currentX / currentS < iLimit)) {
                // index i * S - 1 or i * S
                Rtree::SaveEntryWithOrderIndex(element, r0FileSmall);
            }

            currentX++;
        }
        r0FileSecond.close();
        r0FileSmall.close();
        RectanglesD0Second.clear();
        // END

        orderedInputRectangles.CreateOrderedBoxesOnDisk(onDiskBase + ".boundingbox.d0", onDiskBase + ".boundingbox.d1", xSize);
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
            for(size_t i = 0; i < currentItem.GetOrderedBoxes().GetSize(); i++) {
                rTreeValue box = currentItem.GetOrderedBoxes().GetElementAt(0, i);
                Node leafNode = Node(box.second, box.first);
                currentItem.AddChild(leafNode);
            }
            long long nodePtr = SaveNode(currentItem, true, nodesOfs);
            lookup[currentItem.GetId()] = nodePtr;
        } else {
            std::vector<OrderedBoxes> tgsResult = TGSRecursive(onDiskBase + ".boundingbox." + std::to_string(layer), currentItem.GetOrderedBoxes(), M, std::ceil(((float) currentItem.GetOrderedBoxes().GetSize()) / ((float) M)), this->maxBuildingRamUsage);
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

void OrderedBoxes::CreateOrderedBoxesInRam(multiBoxGeo& rectanglesD0, multiBoxGeo& rectanglesD1) {
    this->workInRam = true;
    this->rectanglesD0InRam = rectanglesD0;
    this->rectanglesD1InRam = rectanglesD1;
    this->size = rectanglesD0.size();
}

void OrderedBoxes::CreateOrderedBoxesOnDisk(const std::string& rectanglesD0, const std::string& rectanglesD1, long long size) {
    this->workInRam = false;
    this->rectanglesD0OnDisk = rectanglesD0 + ".tmp";
    this->rectanglesD1OnDisk = rectanglesD1 + ".tmp";
    this->rectanglesD0SmallOnDisk = rectanglesD0 + ".small.tmp";
    this->rectanglesD1SmallOnDisk = rectanglesD1 + ".small.tmp";
    this->size = size;
}

boxGeo OrderedBoxes::GetBoundingBox() {
    double globalMinX = -1;
    double globalMinY = -1;
    double globalMaxX = -1;
    double globalMaxY = -1;

    if (this->workInRam) {
        for (rTreeValue box : rectanglesD0InRam) {
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
    } else {
        std::ifstream inFile = std::ifstream(this->rectanglesD0OnDisk, std::ios::binary);
        inFile.seekg (0, std::ifstream::end);
        long long fileLength = inFile.tellg();
        inFile.seekg (0, std::ifstream::beg);

        double minX;
        double minY;
        double maxX;
        double maxY;

        while (inFile.tellg() < fileLength) {
            inFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
            inFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
            inFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
            inFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));

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

        inFile.close();
    }

    return Rtree::createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
}

long long OrderedBoxes::GetSize() const {
    return this->size;
}

rTreeValue OrderedBoxes::GetElementAt(size_t dim, long long index) {
    if (this->workInRam) {
        if (dim == 0) {
            return this->rectanglesD0InRam[index];
        }
        return this->rectanglesD1InRam[index];
    }

    std::ifstream inFile;
    if (dim == 0) {
        inFile = std::ifstream(this->rectanglesD0OnDisk, std::ios::binary);
    } else {
        inFile = std::ifstream(this->rectanglesD1OnDisk, std::ios::binary);
    }
    long long fileIndex = index * (sizeof(double) * 4 + sizeof(uint64_t));
    inFile.seekg (fileIndex, std::ifstream::beg);

    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t id;

    inFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

    inFile.close();

    boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
    return std::make_pair(box, id);
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBest(const std::string& filePath, size_t S, long long maxBuildingRamUsage) {
    if (this->workInRam) {
        return this->SplitAtBestInRam(S);
    } else {
        return this->SplitAtBestOnDisk(filePath, S, maxBuildingRamUsage);
    }
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBestInRam(size_t S) {
    double bestCost = -1;
    size_t bestDim = 0;
    unsigned long bestI = 1;
    boxGeo bestB0 = Rtree::createBoundingBox(0, 0, 0, 0);

    for (size_t dim = 0; dim <= 1; dim++) {
        for (size_t i = 1; i < std::ceil(((float) this->size) / ((float) S)); i++) {
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
                minXB0 = (this->rectanglesD0InRam[0].first.min_corner().get<0>() + this->rectanglesD0InRam[0].first.max_corner().get<0>()) / 2;
                maxXB0 = (this->rectanglesD0InRam[i * S - 1].first.min_corner().get<0>() + this->rectanglesD0InRam[i * S - 1].first.max_corner().get<0>()) / 2;

                minXB1 = (this->rectanglesD0InRam[i * S].first.min_corner().get<0>() + this->rectanglesD0InRam[i * S].first.max_corner().get<0>()) / 2;
                maxXB1 = (this->rectanglesD0InRam[this->size - 1].first.min_corner().get<0>() + this->rectanglesD0InRam[this->size - 1].first.max_corner().get<0>()) / 2;
            } else {
                minYB0 = (this->rectanglesD1InRam[0].first.min_corner().get<1>() + this->rectanglesD1InRam[0].first.max_corner().get<1>()) / 2;
                maxYB0 = (this->rectanglesD1InRam[i * S - 1].first.min_corner().get<1>() + this->rectanglesD1InRam[i * S - 1].first.max_corner().get<1>()) / 2;

                minYB1 = (this->rectanglesD1InRam[i * S].first.min_corner().get<1>() + this->rectanglesD1InRam[i * S].first.max_corner().get<1>()) / 2;
                maxYB1 = (this->rectanglesD1InRam[this->size - 1].first.min_corner().get<1>() + this->rectanglesD1InRam[this->size - 1].first.max_corner().get<1>()) / 2;
            }

            boxGeo b0 = Rtree::createBoundingBox(minXB0, minYB0, maxXB0, maxYB0);
            boxGeo b1 = Rtree::createBoundingBox(minXB1, minYB1, maxXB1, maxYB1);


            double cost = costFunctionTGS(b0, b1);

            if (bestCost == -1 || cost < bestCost) {
                bestCost = cost;
                bestDim = dim;
                bestI = i;
                bestB0 = b0;
            }
        }
    }

    OrderedBoxes split0;
    OrderedBoxes split1;

    multiBoxGeo s0Dim0;
    multiBoxGeo s0Dim1;
    multiBoxGeo s1Dim0;
    multiBoxGeo s1Dim1;

    if (bestDim == 0) {
        s0Dim0 = multiBoxGeo(this->rectanglesD0InRam.begin(), this->rectanglesD0InRam.begin() + bestI * S);
        s1Dim0 = multiBoxGeo(this->rectanglesD0InRam.begin() + bestI * S, this->rectanglesD0InRam.end());
    } else {
        s0Dim1 = multiBoxGeo(this->rectanglesD1InRam.begin(), this->rectanglesD1InRam.begin() + bestI * S);
        s1Dim1 = multiBoxGeo(this->rectanglesD1InRam.begin() + bestI * S, this->rectanglesD1InRam.end());
    }

    for (rTreeValue box : bestDim == 0 ? this->rectanglesD1InRam : this->rectanglesD0InRam) {
        pointGeo boxCenter;
        // check if it is exactly on the border. In this case it is not certain that the box belongs to b0
        if (bestDim == 0) {
            boxCenter = make<pointGeo>((box.first.min_corner().get<0>() + box.first.max_corner().get<0>()) / 2, 0.5);
            if (boxCenter.get<0>() == bestB0.max_corner().get<0>()) {
                if (boxInMultiBox(box, s0Dim0)) {
                    s0Dim1.push_back(box);
                } else {
                    s1Dim1.push_back(box);
                }
                continue;
            }
        } else {
            boxCenter = make<pointGeo>(0.5, (box.first.min_corner().get<1>() + box.first.max_corner().get<1>()) / 2);
            if (boxCenter.get<1>() == bestB0.max_corner().get<1>()) {
                if (boxInMultiBox(box, s0Dim1)) {
                    s0Dim0.push_back(box);
                } else {
                    s1Dim0.push_back(box);
                }
                continue;
            }
        }

        if (pointWithinBox(boxCenter, bestB0)) {
            if (bestDim == 0) {
                s0Dim1.push_back(box);
            } else {
                s0Dim0.push_back(box);
            }
        } else {
            if (bestDim == 0) {
                s1Dim1.push_back(box);
            } else {
                s1Dim0.push_back(box);
            }
        }
    }

    split0.CreateOrderedBoxesInRam(s0Dim0, s0Dim1);
    split1.CreateOrderedBoxesInRam(s1Dim0, s1Dim1);

    this->rectanglesD0InRam.clear();
    this->rectanglesD1InRam.clear();

    return std::make_pair(split0, split1);
}

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBestOnDisk(const std::string& filePath, size_t S, long long maxBuildingRamUsage) const {
    OrderedBoxes split0;
    OrderedBoxes split1;

    double bestCost = -1;
    size_t bestDim = 0;
    boxGeo bestB0 = Rtree::createBoundingBox(0, 0, 0, 0);
    long long bestIndex = 0;
    rTreeValueWithOrderIndex bestLastElement;
    rTreeValueWithOrderIndex bestElement;
    rTreeValueWithOrderIndex bestMinElement;
    rTreeValueWithOrderIndex bestMaxElement;

    double minX;
    double minY;
    double maxX;
    double maxY;
    uint64_t id;
    long long orderX;
    long long orderY;

    long long i = 0;

    rTreeValueWithOrderIndex minElement;
    rTreeValueWithOrderIndex maxElement;
    rTreeValueWithOrderIndex currentLastElement;
    rTreeValueWithOrderIndex currentElement;

    bool currentlyAtSTimesI = false;

    // dim 0

    std::ifstream dim0File = std::ifstream(this->rectanglesD0SmallOnDisk, std::ios::binary);
    dim0File.seekg (0, std::ifstream::end);
    long long fileLength = dim0File.tellg();
    dim0File.seekg (0, std::ifstream::beg);

    while (dim0File.tellg() < fileLength) {
        dim0File.read(reinterpret_cast<char*>(&minX), sizeof(double));
        dim0File.read(reinterpret_cast<char*>(&minY), sizeof(double));
        dim0File.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        dim0File.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        dim0File.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        dim0File.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        dim0File.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);

        if (i == 0) {
            // this is the min element
            minElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            i++;
            continue;
        }

        if (i == 1) {
            // this is the max element
            maxElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            i++;
            continue;
        }

        double minXB0 = 0;
        double maxXB0 = 1;
        double minXB1 = 0;
        double maxXB1 = 1;

        if (!currentlyAtSTimesI) {
            currentLastElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            currentlyAtSTimesI = true;
            i++;
            continue;
        }

        if (currentlyAtSTimesI && id != maxElement.first.second) {

            currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

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

        if (bestCost == -1 || cost < bestCost) {
            bestCost = cost;
            bestDim = 0;
            bestB0 = b0;
            bestLastElement = currentLastElement;
            bestElement = currentElement;
            bestMinElement = minElement;
            bestMaxElement = maxElement;
            bestIndex = i;
        }

        i++;
    }
    dim0File.close();

    // dim 1

    i = 0;
    currentlyAtSTimesI = false;

    std::ifstream dim1File = std::ifstream(this->rectanglesD1SmallOnDisk, std::ios::binary);
    dim1File.seekg (0, std::ifstream::end);
    fileLength = dim1File.tellg();
    dim1File.seekg (0, std::ifstream::beg);

    while (dim1File.tellg() < fileLength) {
        dim1File.read(reinterpret_cast<char*>(&minX), sizeof(double));
        dim1File.read(reinterpret_cast<char*>(&minY), sizeof(double));
        dim1File.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        dim1File.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        dim1File.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        dim1File.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        dim1File.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);

        if (i == 0) {
            // this is the min element
            minElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            i++;
            continue;
        }

        if (i == 1) {
            // this is the max element
            maxElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            i++;
            continue;
        }

        double minYB0 = 0;
        double maxYB0 = 1;
        double minYB1 = 0;
        double maxYB1 = 1;

        if (!currentlyAtSTimesI) {
            currentLastElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            currentlyAtSTimesI = true;
            i++;
            continue;
        }

        if (currentlyAtSTimesI && id != maxElement.first.second) {

            currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

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

        if (bestCost == -1 || cost < bestCost) {
            bestCost = cost;
            bestDim = 1;
            bestB0 = b0;
            bestLastElement = currentLastElement;
            bestElement = currentElement;
            bestMinElement = minElement;
            bestMaxElement = maxElement;
            bestIndex = i;
        }

        i++;
    }
    dim1File.close();

    // perfrom the split
    long long sizeLeft = std::ceil((bestIndex - 2) / 2) * S;
    long long sizeRight = this->size - sizeLeft;
    long long split0ByteSize = sizeLeft * (4 * sizeof(double) + sizeof(uint64_t));
    long long split1ByteSize = sizeRight * (4 * sizeof(double) + sizeof(uint64_t));
    bool split0InRam = split0ByteSize * 4 < maxBuildingRamUsage;
    bool split1InRam = split1ByteSize * 4 < maxBuildingRamUsage;

    multiBoxGeo split0Dim0InRam;
    multiBoxGeo split0Dim1InRam;
    multiBoxGeo split1Dim0InRam;
    multiBoxGeo split1Dim1InRam;
    std::ofstream split0Dim0File;
    std::ofstream split0Dim1File;
    std::ofstream split0Dim0FileSmall;
    std::ofstream split0Dim1FileSmall;
    std::ofstream split1Dim0File;
    std::ofstream split1Dim1File;
    std::ofstream split1Dim0FileSmall;
    std::ofstream split1Dim1FileSmall;

    if (!split0InRam) {
        split0Dim0File = std::ofstream(filePath + ".0.dim0.tmp", std::ios::binary);
        split0Dim1File = std::ofstream(filePath + ".0.dim1.tmp", std::ios::binary);
        split0Dim0FileSmall = std::ofstream(filePath + ".0.dim0.small.tmp", std::ios::binary);
        split0Dim1FileSmall = std::ofstream(filePath + ".0.dim1.small.tmp", std::ios::binary);

        if (bestDim == 0) {
            Rtree::SaveEntryWithOrderIndex(bestMinElement, split0Dim0FileSmall);
            Rtree::SaveEntryWithOrderIndex(bestLastElement, split0Dim0FileSmall);
        } else {
            Rtree::SaveEntryWithOrderIndex(bestMinElement, split0Dim1FileSmall);
            Rtree::SaveEntryWithOrderIndex(bestLastElement, split0Dim1FileSmall);
        }
    }

    if (!split1InRam) {
        split1Dim0File = std::ofstream(filePath + ".1.dim0.tmp", std::ios::binary);
        split1Dim1File = std::ofstream(filePath + ".1.dim1.tmp", std::ios::binary);
        split1Dim0FileSmall = std::ofstream(filePath + ".1.dim0.small.tmp", std::ios::binary);
        split1Dim1FileSmall = std::ofstream(filePath + ".1.dim1.small.tmp", std::ios::binary);

        if (bestDim == 0) {
            Rtree::SaveEntryWithOrderIndex(bestElement, split1Dim0FileSmall);
            Rtree::SaveEntryWithOrderIndex(bestMaxElement, split1Dim0FileSmall);
        } else {
            Rtree::SaveEntryWithOrderIndex(bestElement, split1Dim1FileSmall);
            Rtree::SaveEntryWithOrderIndex(bestMaxElement, split1Dim1FileSmall);
        }
    }

    rTreeValueWithOrderIndex minSplit0OtherDim;
    rTreeValueWithOrderIndex maxSplit0OtherDim;
    rTreeValueWithOrderIndex minSplit1OtherDim;
    rTreeValueWithOrderIndex maxSplit1OtherDim;

    std::ifstream dim0WholeFile = std::ifstream(this->rectanglesD0OnDisk, std::ios::binary);
    dim0WholeFile.seekg (0, std::ifstream::end);
    fileLength = dim0WholeFile.tellg();
    dim0WholeFile.seekg (0, std::ifstream::beg);

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
        currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

        long long currentXSplit0 = 0;
        long long currentXSplit1 = 0;
        if ((bestDim == 0 && currentElement.second.first < bestElement.second.first)
        || (bestDim == 1 && currentElement.second.second < bestElement.second.second)) {
            // it's in split 0
            if (split0InRam) {
                split0Dim0InRam.push_back(currentElement.first);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim0File);

                if (bestDim == 0) {
                    if (((currentXSplit0 + 1) % S == 0 && (currentXSplit0 + 1) / S >= 1) || (currentXSplit0 % S == 0 && currentXSplit0 / S >= 1)) {
                        // index i * S - 1 or i * S
                        Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim0FileSmall);
                    }
                } else {
                    if (currentXSplit0 == 0) {
                        minSplit0OtherDim = currentElement;
                        maxSplit0OtherDim = currentElement;
                    }
                    if (currentElement.second.first > maxSplit0OtherDim.second.first) {
                        maxSplit0OtherDim = currentElement;
                    }
                }
                currentXSplit0++;
            }
        } else {
            // it's in split 1
            if (split1InRam) {
                split1Dim0InRam.push_back(currentElement.first);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim0File);

                if (bestDim == 0) {
                    if (((currentXSplit1 + 1) % S == 0 && (currentXSplit1 + 1) / S >= 1) || (currentXSplit1 % S == 0 && currentXSplit1 / S >= 1)) {
                        // index i * S - 1 or i * S
                        Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim0FileSmall);
                    }
                } else {
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
    }
    dim0WholeFile.close();

    std::ifstream dim1WholeFile = std::ifstream(this->rectanglesD1OnDisk, std::ios::binary);
    dim1WholeFile.seekg (0, std::ifstream::end);
    fileLength = dim1WholeFile.tellg();
    dim1WholeFile.seekg (0, std::ifstream::beg);

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
        currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

        long long currentYSplit0 = 0;
        long long currentYSplit1 = 0;
        if ((bestDim == 0 && currentElement.second.first < bestElement.second.first)
            || (bestDim == 1 && currentElement.second.second < bestElement.second.second)) {
            // it's in split 0
            if (split0InRam) {
                split0Dim1InRam.push_back(currentElement.first);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim1File);

                if (bestDim == 1) {
                    if (((currentYSplit0 + 1) % S == 0 && (currentYSplit0 + 1) / S >= 1) || (currentYSplit0 % S == 0 && currentYSplit0 / S >= 1)) {
                        // index i * S - 1 or i * S
                        Rtree::SaveEntryWithOrderIndex(currentElement, split0Dim1FileSmall);
                    }
                } else {
                    if (currentYSplit0 == 0) {
                        minSplit0OtherDim = currentElement;
                        maxSplit0OtherDim = currentElement;
                    }
                    if (currentElement.second.second > maxSplit0OtherDim.second.second) {
                        maxSplit0OtherDim = currentElement;
                    }
                }
                currentYSplit0++;
            }
        } else {
            // it's in split 1
            if (split1InRam) {
                split1Dim1InRam.push_back(currentElement.first);
            } else {
                Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim1File);

                if (bestDim == 1) {
                    if (((currentYSplit1 + 1) % S == 0 && (currentYSplit1 + 1) / S >= 1) || (currentYSplit1 % S == 0 && currentYSplit1 / S >= 1)) {
                        // index i * S - 1 or i * S
                        Rtree::SaveEntryWithOrderIndex(currentElement, split1Dim1FileSmall);
                    }
                } else {
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
    }
    dim1WholeFile.close();

    if (!split0InRam) {
        split0Dim0File.close();
        split0Dim1File.close();

        std::ifstream split0WholeFile;
        if (bestDim == 0) {
            // load the other dim
            split0WholeFile = std::ifstream(filePath + ".0.dim1.tmp", std::ios::binary);
        } else {
            split0WholeFile = std::ifstream(filePath + ".0.dim0.tmp", std::ios::binary);
        }

        // create the missing small file
        Rtree::SaveEntryWithOrderIndex(minSplit0OtherDim, bestDim == 0 ? split0Dim1FileSmall : split0Dim0FileSmall);
        Rtree::SaveEntryWithOrderIndex(maxSplit0OtherDim, bestDim == 0 ? split0Dim1FileSmall : split0Dim0FileSmall);

        split0WholeFile.seekg (0, std::ifstream::end);
        fileLength = split0WholeFile.tellg();
        split0WholeFile.seekg (0, std::ifstream::beg);

        long long currentIndexSplit0 = 0;
        while (split0WholeFile.tellg() < fileLength) {
            split0WholeFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
            split0WholeFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
            split0WholeFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
            split0WholeFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));
            split0WholeFile.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
            split0WholeFile.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
            split0WholeFile.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

            boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
            rTreeValue boxWithId = std::make_pair(box, id);
            currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            if (((currentIndexSplit0 + 1) % S == 0 && (currentIndexSplit0 + 1) / S >= 1) || (currentIndexSplit0 % S == 0 && currentIndexSplit0 / S >= 1)) {
                // index i * S - 1 or i * S
                Rtree::SaveEntryWithOrderIndex(currentElement, bestDim == 0 ? split0Dim1FileSmall : split0Dim0FileSmall);
            }
            currentIndexSplit0++;
        }
        split0WholeFile.close();
        split0Dim0FileSmall.close();
        split0Dim1FileSmall.close();

        split0.CreateOrderedBoxesOnDisk(filePath + ".0.dim0", filePath + ".0.dim1", sizeLeft);
    } else {
        split0.CreateOrderedBoxesInRam(split0Dim0InRam, split0Dim1InRam);
    }

    if (!split1InRam) {
        split1Dim0File.close();
        split1Dim1File.close();

        std::ifstream split1WholeFile;
        if (bestDim == 0) {
            // load the other dim
            split1WholeFile = std::ifstream(filePath + ".1.dim1.tmp", std::ios::binary);
        } else {
            split1WholeFile = std::ifstream(filePath + ".1.dim0.tmp", std::ios::binary);
        }

        // create the missing small file
        Rtree::SaveEntryWithOrderIndex(minSplit1OtherDim, bestDim == 0 ? split1Dim1FileSmall : split1Dim0FileSmall);
        Rtree::SaveEntryWithOrderIndex(maxSplit1OtherDim, bestDim == 0 ? split1Dim1FileSmall : split1Dim0FileSmall);

        split1WholeFile.seekg (0, std::ifstream::end);
        fileLength = split1WholeFile.tellg();
        split1WholeFile.seekg (0, std::ifstream::beg);

        long long currentIndexSplit1 = 0;
        while (split1WholeFile.tellg() < fileLength) {
            split1WholeFile.read(reinterpret_cast<char*>(&minX), sizeof(double));
            split1WholeFile.read(reinterpret_cast<char*>(&minY), sizeof(double));
            split1WholeFile.read(reinterpret_cast<char*>(&maxX), sizeof(double));
            split1WholeFile.read(reinterpret_cast<char*>(&maxY), sizeof(double));
            split1WholeFile.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
            split1WholeFile.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
            split1WholeFile.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

            boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
            rTreeValue boxWithId = std::make_pair(box, id);
            currentElement = std::make_pair(boxWithId, std::make_pair(orderX, orderY));
            if (((currentIndexSplit1 + 1) % S == 0 && (currentIndexSplit1 + 1) / S >= 1) || (currentIndexSplit1 % S == 0 && currentIndexSplit1 / S >= 1)) {
                // index i * S - 1 or i * S
                Rtree::SaveEntryWithOrderIndex(currentElement, bestDim == 0 ? split1Dim1FileSmall : split1Dim0FileSmall);
            }
            currentIndexSplit1++;
        }
        split1WholeFile.close();
        split1Dim0FileSmall.close();
        split1Dim1FileSmall.close();

        split1.CreateOrderedBoxesOnDisk(filePath + ".1.dim0", filePath + ".1.dim1", sizeRight);
    } else {
        split1.CreateOrderedBoxesInRam(split1Dim0InRam, split1Dim1InRam);
    }

    std::remove(this->rectanglesD0OnDisk.c_str());
    std::remove(this->rectanglesD0SmallOnDisk.c_str());
    std::remove(this->rectanglesD1OnDisk.c_str());
    std::remove(this->rectanglesD1SmallOnDisk.c_str());

    return std::make_pair(split0, split1);
}
