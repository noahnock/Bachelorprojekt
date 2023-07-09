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

    if (n <= S) {
        // stop condition
        return std::vector<OrderedBoxes> { orderedInputRectangles };
    }

    double bestCost = -1;
    size_t bestDim = 0;
    unsigned long bestI = 1;
    boxGeo bestB0 = Rtree::createBoundingBox(0, 0, 0, 0);

    centerOrdering(orderedInputRectangles.rectanglesD0InRam, 0);
    centerOrdering(orderedInputRectangles.rectanglesD1InRam, 1);

    for (size_t dim = 0; dim <= 1; dim++) {
        for (size_t i = 1; i < std::ceil(((float) n) / ((float) S)); i++) {
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
                minXB0 = (orderedInputRectangles.GetElementAt(dim, 0).first.min_corner().get<0>() + orderedInputRectangles.GetElementAt(dim, 0).first.max_corner().get<0>()) / 2;
                maxXB0 = (orderedInputRectangles.GetElementAt(dim, i * S - 1).first.min_corner().get<0>() + orderedInputRectangles.GetElementAt(dim, i * S - 1).first.max_corner().get<0>()) / 2;

                minXB1 = (orderedInputRectangles.GetElementAt(dim, i * S).first.min_corner().get<0>() + orderedInputRectangles.GetElementAt(dim, i * S).first.max_corner().get<0>()) / 2;
                maxXB1 = (orderedInputRectangles.GetElementAt(dim, orderedInputRectangles.GetSize() - 1).first.min_corner().get<0>() + orderedInputRectangles.GetElementAt(dim, orderedInputRectangles.GetSize() - 1).first.max_corner().get<0>()) / 2;
            } else {
                minYB0 = (orderedInputRectangles.GetElementAt(dim, 0).first.min_corner().get<1>() + orderedInputRectangles.GetElementAt(dim, 0).first.max_corner().get<1>()) / 2;
                maxYB0 = (orderedInputRectangles.GetElementAt(dim, i * S - 1).first.min_corner().get<1>() + orderedInputRectangles.GetElementAt(dim, i * S - 1).first.max_corner().get<1>()) / 2;

                minYB1 = (orderedInputRectangles.GetElementAt(dim, i * S).first.min_corner().get<1>() + orderedInputRectangles.GetElementAt(dim, i * S).first.max_corner().get<1>()) / 2;
                maxYB1 = (orderedInputRectangles.GetElementAt(dim, orderedInputRectangles.GetSize() - 1).first.min_corner().get<1>() + orderedInputRectangles.GetElementAt(dim, orderedInputRectangles.GetSize() - 1).first.max_corner().get<1>()) / 2;
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

    // split the rectangles
    std::pair<OrderedBoxes, OrderedBoxes> split = orderedInputRectangles.SplitAt(filePath, bestDim, bestI * S, bestB0, maxBuildingRamUsage);

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
        multiBoxGeo inputRectangles = LoadEntries(file);

        multiBoxGeo RectanglesD0(inputRectangles);
        multiBoxGeo RectanglesD1(inputRectangles);

        centerOrdering(RectanglesD0, 0);
        centerOrdering(RectanglesD1, 1);

        long long size = 0;

        std::ofstream r0File = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
        for (rTreeValue element : RectanglesD0) {
            Rtree::SaveEntry(element.first, element.second, r0File);
            size++;
        }
        r0File.close();
        std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
        for (rTreeValue element : RectanglesD1) {
            Rtree::SaveEntry(element.first, element.second, r1File);
        }
        r1File.close();
        RectanglesD0.clear();
        RectanglesD1.clear();
        // END

        orderedInputRectangles.CreateOrderedBoxesOnDisk(onDiskBase + ".boundingbox.d0.tmp", onDiskBase + ".boundingbox.d1.tmp", size);
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
        long long test = loadEntriesIfs.tellg();
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
    this->rectanglesD0OnDisk = rectanglesD0;
    this->rectanglesD1OnDisk = rectanglesD1;
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

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAt(const std::string& filePath, size_t dim, long long index, boxGeo boundingBox, long long maxBuildingRamUsage) {
    OrderedBoxes split0;
    OrderedBoxes split1;

    if (this->workInRam) {
        multiBoxGeo s0BestDim;
        multiBoxGeo s1BestDim;

        if (dim == 0) {
            s0BestDim = multiBoxGeo(this->rectanglesD0InRam.begin(), this->rectanglesD0InRam.begin() + index);
            s1BestDim = multiBoxGeo(this->rectanglesD0InRam.begin() + index, this->rectanglesD0InRam.end());
        } else {
            s0BestDim = multiBoxGeo(this->rectanglesD1InRam.begin(), this->rectanglesD1InRam.begin() + index);
            s1BestDim = multiBoxGeo(this->rectanglesD1InRam.begin() + index, this->rectanglesD1InRam.end());
        }

        multiBoxGeo s0OtherDim;
        multiBoxGeo s1OtherDim;

        for (rTreeValue box : dim == 0 ? this->rectanglesD1InRam : this->rectanglesD0InRam) {
            pointGeo boxCenter;
            // check if it is exactly on the border. In this case it is not certain that the box belongs to b0
            if (dim == 0) {
                boxCenter = make<pointGeo>((box.first.min_corner().get<0>() + box.first.max_corner().get<0>()) / 2, 0.5);
                if (boxCenter.get<0>() == boundingBox.max_corner().get<0>()) {
                    if (boxInMultiBox(box, s0BestDim)) {
                        s0OtherDim.push_back(box);
                    } else {
                        s1OtherDim.push_back(box);
                    }
                    continue;
                }
            } else {
                boxCenter = make<pointGeo>(0.5, (box.first.min_corner().get<1>() + box.first.max_corner().get<1>()) / 2);
                if (boxCenter.get<1>() == boundingBox.max_corner().get<1>()) {
                    if (boxInMultiBox(box, s0BestDim)) {
                        s0OtherDim.push_back(box);
                    } else {
                        s1OtherDim.push_back(box);
                    }
                    continue;
                }
            }

            if (pointWithinBox(boxCenter, boundingBox)) {
                s0OtherDim.push_back(box);
            } else {
                s1OtherDim.push_back(box);
            }
        }

        split0.CreateOrderedBoxesInRam(s0BestDim, s0OtherDim);
        split1.CreateOrderedBoxesInRam(s1BestDim, s1OtherDim);

        this->rectanglesD0InRam.clear();
        this->rectanglesD1InRam.clear();
    } else {
        // on disk
        std::ifstream fileBestDim;
        std::ofstream bestDimSplit0;
        std::ofstream bestDimSplit1;
        std::ofstream otherDimSplit0;
        std::ofstream otherDimSplit1;

        multiBoxGeo bestDimSplit0OnRam;
        multiBoxGeo bestDimSplit1OnRam;
        multiBoxGeo otherDimSplit0OnRam;
        multiBoxGeo otherDimSplit1OnRam;

        long long split0ByteSize = index * (sizeof(double) * 4 + sizeof(uint64_t));
        long long split1ByteSize = this->size * (sizeof(double) * 4 + sizeof(uint64_t)) - split0ByteSize;
        bool split0InRam = split0ByteSize * 4 < maxBuildingRamUsage;
        bool split1InRam = split1ByteSize * 4 < maxBuildingRamUsage;

        if (dim == 0) {
            fileBestDim = std::ifstream(this->rectanglesD0OnDisk, std::ios::binary);
            if (!split0InRam) {
                bestDimSplit0 = std::ofstream(filePath + ".0.dim0.tmp", std::ios::binary);
                otherDimSplit0 = std::ofstream(filePath + ".0.dim1.tmp", std::ios::binary);
            }
            if (!split1InRam) {
                bestDimSplit1 = std::ofstream(filePath + ".1.dim0.tmp", std::ios::binary);
                otherDimSplit1 = std::ofstream(filePath + ".1.dim1.tmp", std::ios::binary);
            }
        } else {
            fileBestDim = std::ifstream(this->rectanglesD1OnDisk, std::ios::binary);
            if (!split0InRam) {
                bestDimSplit0 = std::ofstream(filePath + ".0.dim1.tmp", std::ios::binary);
                otherDimSplit0 = std::ofstream(filePath + ".0.dim0.tmp", std::ios::binary);
            }
            if (!split1InRam) {
                bestDimSplit1 = std::ofstream(filePath + ".1.dim1.tmp", std::ios::binary);
                otherDimSplit1 = std::ofstream(filePath + ".1.dim0.tmp", std::ios::binary);
            }
        }
        fileBestDim.seekg (0, std::ifstream::end);
        long long fileLength = fileBestDim.tellg();
        fileBestDim.seekg (0, std::ifstream::beg);

        double minX;
        double minY;
        double maxX;
        double maxY;
        uint64_t id;

        long long currentElementIndex = 0;

        while (fileBestDim.tellg() < fileLength) {
            fileBestDim.read(reinterpret_cast<char*>(&minX), sizeof(double));
            fileBestDim.read(reinterpret_cast<char*>(&minY), sizeof(double));
            fileBestDim.read(reinterpret_cast<char*>(&maxX), sizeof(double));
            fileBestDim.read(reinterpret_cast<char*>(&maxY), sizeof(double));
            fileBestDim.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

            if (currentElementIndex < index) {
                if (!split0InRam) {
                    bestDimSplit0.write(reinterpret_cast<const char *>(&minX), sizeof(double));
                    bestDimSplit0.write(reinterpret_cast<const char *>(&minY), sizeof(double));
                    bestDimSplit0.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
                    bestDimSplit0.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
                    bestDimSplit0.write(reinterpret_cast<const char *>(&id), sizeof(uint64_t));
                } else {
                    boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
                    rTreeValue boxWithId = std::make_pair(box, id);
                    bestDimSplit0OnRam.push_back(boxWithId);
                }
            } else {
                if (!split1InRam) {
                    bestDimSplit1.write(reinterpret_cast<const char *>(&minX), sizeof(double));
                    bestDimSplit1.write(reinterpret_cast<const char *>(&minY), sizeof(double));
                    bestDimSplit1.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
                    bestDimSplit1.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
                    bestDimSplit1.write(reinterpret_cast<const char *>(&id), sizeof(uint64_t));
                } else {
                    boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
                    rTreeValue boxWithId = std::make_pair(box, id);
                    bestDimSplit1OnRam.push_back(boxWithId);
                }
            }
            currentElementIndex++;
        }

        fileBestDim.close();
        if (!split0InRam) {
            bestDimSplit0.close();
        }
        if (!split1InRam) {
            bestDimSplit1.close();
        }

        //TODO sort on disk for other dimension
        // START
        multiBoxGeo rectSplit0;
        multiBoxGeo rectSplit1;
        if (!split0InRam) {
            rectSplit0 = Rtree::LoadEntries(dim == 0 ? filePath + ".0.dim0.tmp" : filePath + ".0.dim1.tmp");
        } else {
            rectSplit0 = multiBoxGeo(bestDimSplit0OnRam);
        }
        if (!split1InRam) {
            rectSplit1 = Rtree::LoadEntries(dim == 0 ? filePath + ".1.dim0.tmp" : filePath + ".1.dim1.tmp");
        } else {
            rectSplit1 = multiBoxGeo(bestDimSplit1OnRam);
        }

        centerOrdering(rectSplit0, 1 - dim);
        centerOrdering(rectSplit1, 1 - dim);

        for (rTreeValue element : rectSplit0) {
            if (!split0InRam) {
                Rtree::SaveEntry(element.first, element.second, otherDimSplit0);
            } else {
                otherDimSplit0OnRam.push_back(element);
            }
        }

        for (rTreeValue element : rectSplit1) {
            if (!split1InRam) {
                Rtree::SaveEntry(element.first, element.second, otherDimSplit1);
            } else {
                otherDimSplit1OnRam.push_back(element);
            }
        }

        if (!split0InRam) {
            otherDimSplit0.close();
        }
        if (!split1InRam) {
            otherDimSplit1.close();
        }
        rectSplit0.clear();
        rectSplit1.clear();
        // END

        std::remove(this->rectanglesD0OnDisk.c_str());
        std::remove(this->rectanglesD1OnDisk.c_str());
        if (!split0InRam) {
            split0.CreateOrderedBoxesOnDisk(filePath + ".0.dim0.tmp", filePath + ".0.dim1.tmp", index);
        } else {
            split0.CreateOrderedBoxesInRam(bestDimSplit0OnRam, otherDimSplit0OnRam);
        }

        if (!split1InRam) {
            split1.CreateOrderedBoxesOnDisk(filePath + ".1.dim0.tmp", filePath + ".1.dim1.tmp", this->size - index);
        } else {
            split1.CreateOrderedBoxesInRam(bestDimSplit1OnRam, otherDimSplit1OnRam);
        }
    }

    return std::make_pair(split0, split1);
}