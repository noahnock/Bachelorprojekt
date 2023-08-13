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

    // recursion
    std::vector<OrderedBoxes> result0 = TGSRecursive(filePath + ".0", split.first, M, S, maxBuildingRamUsage);
    std::vector<OrderedBoxes> result1 = TGSRecursive(filePath + ".1", split.second, M, S, maxBuildingRamUsage);

    std::vector<OrderedBoxes> result;
    result.insert(result.begin(), result0.begin(), result0.end());
    result.insert(result.end(), result1.begin(), result1.end());

    return result;
}

void Rtree::BuildTree(const std::string& onDiskBase, size_t M, const std::string& folder) const {
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

            if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
                || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
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
            if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
                || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
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

            if (((ySize + 1) % currentS == 0 && (ySize + 1) / currentS >= 1 && (ySize + 1) / currentS < M)
            || (ySize % currentS == 0 && ySize / currentS >= 1 && ySize / currentS < M)) {
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

            if (((currentX + 1) % currentS == 0 && (currentX + 1) / currentS >= 1 && (currentX + 1) / currentS < M)
                || (currentX % currentS == 0 && currentX / currentS >= 1 && currentX / currentS < M)) {
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
        for(rTreeValueWithOrderIndex box : *this->GetOrderedBoxes().GetRectanglesInRam()) {
            Node leafNode = Node(box.first.second, box.first.first);
            this->AddChild(leafNode);
        }
    } else {
        FileReader fileReader = FileReader(this->GetOrderedBoxes().GetRectanglesOnDisk());

        std::optional<rTreeValueWithOrderIndex> element = fileReader.GetNextElement();
        while(element) {
            Node leafNode = Node(element.value().first.second, element.value().first.first);
            this->AddChild(leafNode);
            element = fileReader.GetNextElement();
        }

        fileReader.Close();
    }
}

long long Node::GetId() const {
    return this->id;
}

OrderedBoxes ConstructionNode::GetOrderedBoxes() {
    return this->orderedBoxes;
}

Node::Node(long long id, boxGeo boundingbox) {
    this->id = id;
    this->boundingBox = boundingbox;
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

bool Node::GetIsLastInnerNode() const {
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
    FileReader fileReader = FileReader(file);

    std::optional<rTreeValueWithOrderIndex> element = fileReader.GetNextElement();
    while (element) {
        boxes.push_back(element.value());
        element = fileReader.GetNextElement();
    }

    fileReader.Close();

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

std::shared_ptr<multiBoxWithOrderIndex> OrderedBoxes::GetRectanglesInRam() {
    return this->rectanglesD0InRam;
}

std::string OrderedBoxes::GetRectanglesOnDisk() {
    return this->rectanglesD0OnDisk;
}

SplitResult OrderedBoxes::GetBestSplit() {
    struct SplitResult splitResult;

    rTreeValueWithOrderIndex minElement;
    rTreeValueWithOrderIndex maxElement;
    rTreeValueWithOrderIndex currentLastElement;
    rTreeValueWithOrderIndex currentElement;

    bool currentlyAtSTimesI = false;

    for (size_t dim = 0; dim < 2; dim++) {
        for (long long i = 0; i < this->rectanglesD0Small->size(); i++) {
            currentElement = dim == 0 ? (*this->rectanglesD0Small)[i] : (*this->rectanglesD1Small)[i];

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

            if (!currentlyAtSTimesI) {
                currentLastElement = currentElement;
                currentlyAtSTimesI = true;
                continue;
            }

            double minXB0 = 0;
            double maxXB0 = 1;
            double minXB1 = 0;
            double maxXB1 = 1;
            double minYB0 = 0;
            double maxYB0 = 1;
            double minYB1 = 0;
            double maxYB1 = 1;

            if (currentlyAtSTimesI && currentElement.first.second != maxElement.first.second) {

                if (dim == 0) {
                    minXB0 = (minElement.first.first.min_corner().get<0>() + minElement.first.first.max_corner().get<0>()) / 2;
                    maxXB0 = (currentLastElement.first.first.min_corner().get<0>() + currentLastElement.first.first.max_corner().get<0>()) / 2;

                    minXB1 = (currentElement.first.first.min_corner().get<0>() + currentElement.first.first.max_corner().get<0>()) / 2;
                    maxXB1 = (maxElement.first.first.min_corner().get<0>() + maxElement.first.first.max_corner().get<0>()) / 2;
                } else {
                    minYB0 = (minElement.first.first.min_corner().get<1>() + minElement.first.first.max_corner().get<1>()) / 2;
                    maxYB0 = (currentLastElement.first.first.min_corner().get<1>() + currentLastElement.first.first.max_corner().get<1>()) / 2;

                    minYB1 = (currentElement.first.first.min_corner().get<1>() + currentElement.first.first.max_corner().get<1>()) / 2;
                    maxYB1 = (maxElement.first.first.min_corner().get<1>() + maxElement.first.first.max_corner().get<1>()) / 2;
                }

                currentlyAtSTimesI = false;
            } else {
                break;
            }

            boxGeo b0 = Rtree::createBoundingBox(minXB0, minYB0, maxXB0, maxYB0);
            boxGeo b1 = Rtree::createBoundingBox(minXB1, minYB1, maxXB1, maxYB1);


            double cost = costFunctionTGS(b0, b1);

            if (splitResult.bestCost == -1 || cost < splitResult.bestCost) {
                splitResult.bestCost = cost;
                splitResult.bestDim = dim;
                splitResult.bestLastElement = currentLastElement;
                splitResult.bestElement = currentElement;
                splitResult.bestMinElement = minElement;
                splitResult.bestMaxElement = maxElement;
                splitResult.bestIndex = i;
            }
        }
        currentlyAtSTimesI = false;
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

std::pair<OrderedBoxes, OrderedBoxes> OrderedBoxes::SplitAtBestInRam(size_t S, size_t M) {
    struct SplitResult splitResult = this->GetBestSplit();

    OrderedBoxes split0;
    OrderedBoxes split1;

    struct SplitBuffersRam splitBuffers;

    splitBuffers.s0Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s0Dim1 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s1Dim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s1Dim1 = std::make_shared<multiBoxWithOrderIndex>();

    splitBuffers.s0SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s0SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s1SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffers.s1SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();

    std::pair<boxGeo, boxGeo> boundingBoxes = PerformSplit(splitResult, splitBuffers, M, S);

    split0.CreateOrderedBoxesInRam(splitBuffers.s0Dim0, splitBuffers.s0Dim1, splitBuffers.s0SmallDim0, splitBuffers.s0SmallDim1, boundingBoxes.first);
    split1.CreateOrderedBoxesInRam(splitBuffers.s1Dim0, splitBuffers.s1Dim1, splitBuffers.s1SmallDim0, splitBuffers.s1SmallDim1, boundingBoxes.second);

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

    struct SplitBuffersDisk splitBuffers;
    struct SplitBuffersRam splitBuffersRam;

    // perfrom the split
    long long sizeLeft = std::ceil((splitResult.bestIndex - 2) / 2.0) * S;
    long long sizeRight = this->size - sizeLeft;
    long long split0ByteSize = sizeLeft * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    long long split1ByteSize = sizeRight * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    bool split0InRam = split0ByteSize * 4 < maxBuildingRamUsage;
    bool split1InRam = split1ByteSize * 4 < maxBuildingRamUsage;

    splitBuffersRam.s0SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffersRam.s0SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffersRam.s1SmallDim0 = std::make_shared<multiBoxWithOrderIndex>();
    splitBuffersRam.s1SmallDim1 = std::make_shared<multiBoxWithOrderIndex>();

    if (!split0InRam) {
        splitBuffers.split0Dim0File = { std::ofstream(filePath + ".0.dim0.tmp", std::ios::binary) };
        splitBuffers.split0Dim1File = { std::ofstream(filePath + ".0.dim1.tmp", std::ios::binary) };
    } else {
        splitBuffersRam.s0Dim0 = std::make_shared<multiBoxWithOrderIndex>();
        splitBuffersRam.s0Dim1 = std::make_shared<multiBoxWithOrderIndex>();
    }

    if (!split1InRam) {
        splitBuffers.split1Dim0File = { std::ofstream(filePath + ".1.dim0.tmp", std::ios::binary) };
        splitBuffers.split1Dim1File = { std::ofstream(filePath + ".1.dim1.tmp", std::ios::binary) };
    } else {
        splitBuffersRam.s1Dim0 = std::make_shared<multiBoxWithOrderIndex>();
        splitBuffersRam.s1Dim1 = std::make_shared<multiBoxWithOrderIndex>();
    }

    splitBuffers.splitBuffersRam = splitBuffersRam;

    std::pair<boxGeo, boxGeo> boundingBoxes = PerformSplit(splitResult, splitBuffers, M, S, maxBuildingRamUsage);

    if (!split0InRam) {
        splitBuffers.split0Dim0File.value().close();
        splitBuffers.split0Dim1File.value().close();

        split0.CreateOrderedBoxesOnDisk(filePath + ".0.dim0", filePath + ".0.dim1", splitBuffers.splitBuffersRam.s0SmallDim0, splitBuffers.splitBuffersRam.s0SmallDim1, sizeLeft, boundingBoxes.first);
    } else {
        split0.CreateOrderedBoxesInRam(splitBuffers.splitBuffersRam.s0Dim0, splitBuffers.splitBuffersRam.s0Dim1, splitBuffers.splitBuffersRam.s0SmallDim0, splitBuffers.splitBuffersRam.s0SmallDim1, boundingBoxes.first);
    }

    if (!split1InRam) {
        splitBuffers.split1Dim0File.value().close();
        splitBuffers.split1Dim1File.value().close();

        split1.CreateOrderedBoxesOnDisk(filePath + ".1.dim0", filePath + ".1.dim1", splitBuffers.splitBuffersRam.s1SmallDim0, splitBuffers.splitBuffersRam.s1SmallDim1, sizeRight, boundingBoxes.second);
    } else {
        split1.CreateOrderedBoxesInRam(splitBuffers.splitBuffersRam.s1Dim0, splitBuffers.splitBuffersRam.s1Dim1, splitBuffers.splitBuffersRam.s1SmallDim0, splitBuffers.splitBuffersRam.s1SmallDim1, boundingBoxes.second);
    }

    std::remove(this->rectanglesD0OnDisk.c_str());
    std::remove(this->rectanglesD1OnDisk.c_str());

    return std::make_pair(split0, split1);
}

std::pair<boxGeo, boxGeo> OrderedBoxes::PerformSplit(SplitResult splitResult, SplitBuffersRam& splitBuffersRam, size_t M, size_t S) {
    struct SplitBuffersDisk splitBuffersDisk;

    splitBuffersDisk.splitBuffersRam = splitBuffersRam;
    splitBuffersDisk.split0Dim0File = {};
    splitBuffersDisk.split0Dim1File = {};
    splitBuffersDisk.split1Dim0File = {};
    splitBuffersDisk.split1Dim1File = {};

    std::pair<boxGeo, boxGeo> boundingBoxes = PerformSplit(splitResult, splitBuffersDisk, M, S, 0);

    splitBuffersRam = splitBuffersDisk.splitBuffersRam;

    return boundingBoxes;
}

std::pair<boxGeo, boxGeo> OrderedBoxes::PerformSplit(SplitResult splitResult, SplitBuffersDisk& splitBuffers, size_t M, size_t S, long long maxBuildingRamUsage) {
    long long sizeLeft = std::ceil((splitResult.bestIndex - 2) / 2.0) * S;
    long long sizeRight = this->size - sizeLeft;
    size_t SSplit0 = sizeLeft <= S ? std::ceil(sizeLeft / (double) M) : S;
    size_t SSplit1 = sizeRight <= S ? std::ceil(sizeRight / (double) M) : S;
    long long split0ByteSize = sizeLeft * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    long long split1ByteSize = sizeRight * (4 * sizeof(double) + sizeof(uint64_t) + 2 * sizeof(long long));
    bool split0InRam = maxBuildingRamUsage == 0 || split0ByteSize * 4 < maxBuildingRamUsage;
    bool split1InRam = maxBuildingRamUsage == 0 || split1ByteSize * 4 < maxBuildingRamUsage;

    double globalMinXS0 = -1;
    double globalMinYS0 = -1;
    double globalMaxXS0 = -1;
    double globalMaxYS0 = -1;

    double globalMinXS1 = -1;
    double globalMinYS1 = -1;
    double globalMaxXS1 = -1;
    double globalMaxYS1 = -1;

    rTreeValueWithOrderIndex minSplit0OtherDim;
    rTreeValueWithOrderIndex maxSplit0OtherDim;
    rTreeValueWithOrderIndex minSplit1OtherDim;
    rTreeValueWithOrderIndex maxSplit1OtherDim;

    if (splitResult.bestDim == 0) {
        splitBuffers.splitBuffersRam.s0SmallDim0->push_back(splitResult.bestMinElement);
        splitBuffers.splitBuffersRam.s0SmallDim0->push_back(splitResult.bestLastElement);
        splitBuffers.splitBuffersRam.s1SmallDim0->push_back(splitResult.bestElement);
        splitBuffers.splitBuffersRam.s1SmallDim0->push_back(splitResult.bestMaxElement);

        // placeholder
        splitBuffers.splitBuffersRam.s0SmallDim1->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s0SmallDim1->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s1SmallDim1->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s1SmallDim1->push_back(rTreeValueWithOrderIndex());
    } else {
        splitBuffers.splitBuffersRam.s0SmallDim1->push_back(splitResult.bestMinElement);
        splitBuffers.splitBuffersRam.s0SmallDim1->push_back(splitResult.bestLastElement);
        splitBuffers.splitBuffersRam.s1SmallDim1->push_back(splitResult.bestElement);
        splitBuffers.splitBuffersRam.s1SmallDim1->push_back(splitResult.bestMaxElement);

        // placeholder
        splitBuffers.splitBuffersRam.s0SmallDim0->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s0SmallDim0->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s1SmallDim0->push_back(rTreeValueWithOrderIndex());
        splitBuffers.splitBuffersRam.s1SmallDim0->push_back(rTreeValueWithOrderIndex());
    }

    std::optional<rTreeValueWithOrderIndex> elementOpt;
    std::optional<FileReader> fileReaderDim0;
    std::optional<FileReader> fileReaderDim1;
    if (!this->workInRam) {
        fileReaderDim0 = { FileReader(this->rectanglesD0OnDisk) };
        fileReaderDim1 = { FileReader(this->rectanglesD1OnDisk) };
    }
    long long currentXSplit0 = 0;
    long long currentXSplit1 = 0;
    long long currentYSplit0 = 0;
    long long currentYSplit1 = 0;
    for (size_t dim = 0; dim < 2; dim++) {
        long long i = 0;

        if (!this->workInRam) {
            if (dim == 0)
                elementOpt = fileReaderDim0.value().GetNextElement();
            if (dim == 1)
                elementOpt = fileReaderDim1.value().GetNextElement();
        }

        while ((this->workInRam && i < this->size) || (!this->workInRam && elementOpt)) {
            rTreeValueWithOrderIndex element;
            if (this->workInRam) {
                element = dim == 0 ? (*this->rectanglesD0InRam)[i] : (*this->rectanglesD1InRam)[i];
            } else {
                element = elementOpt.value();
            }

            if ((splitResult.bestDim == 0 && element.second.first < splitResult.bestElement.second.first)
                || (splitResult.bestDim == 1 && element.second.second < splitResult.bestElement.second.second)) {
                // split 0

                if (dim == 0) {
                    if (split0InRam || this->workInRam) {
                        splitBuffers.splitBuffersRam.s0Dim0->push_back(element);
                    } else {
                        Rtree::SaveEntryWithOrderIndex(element, splitBuffers.split0Dim0File.value());
                    }
                    if (((currentXSplit0 + 1) % SSplit0 == 0 && (currentXSplit0 + 1) / SSplit0 >= 1 && (currentXSplit0 + 1) / SSplit0 < M)
                        || (currentXSplit0 % SSplit0 == 0 && currentXSplit0 / SSplit0 >= 1 && currentXSplit0 / SSplit0 < M)) {
                        // index i * S - 1 or i * S
                        splitBuffers.splitBuffersRam.s0SmallDim0->push_back(element);
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

                    if (splitResult.bestDim == 1) {
                        if (currentXSplit0 == 0) {
                            minSplit0OtherDim = element;
                            maxSplit0OtherDim = element;
                        }
                        if (element.second.first > maxSplit0OtherDim.second.first) {
                            maxSplit0OtherDim = element;
                        }
                    }

                    currentXSplit0++;
                } else {
                    if (split0InRam || this->workInRam) {
                        splitBuffers.splitBuffersRam.s0Dim1->push_back(element);
                    } else {
                        Rtree::SaveEntryWithOrderIndex(element, splitBuffers.split0Dim1File.value());
                    }

                    if (((currentYSplit0 + 1) % SSplit0 == 0 && (currentYSplit0 + 1) / SSplit0 >= 1 && (currentYSplit0 + 1) / SSplit0 < M)
                        || (currentYSplit0 % SSplit0 == 0 && currentYSplit0 / SSplit0 >= 1 && currentYSplit0 / SSplit0 < M)) {
                        // index i * S - 1 or i * S
                        splitBuffers.splitBuffersRam.s0SmallDim1->push_back(element);
                    }

                    if (splitResult.bestDim == 0) {
                        if (currentYSplit0 == 0) {
                            minSplit0OtherDim = element;
                            maxSplit0OtherDim = element;
                        }
                        if (element.second.first > maxSplit0OtherDim.second.first) {
                            maxSplit0OtherDim = element;
                        }
                    }

                    currentYSplit0++;
                }
            } else {
                // split 1

                if (dim == 0) {
                    if (split1InRam || this->workInRam) {
                        splitBuffers.splitBuffersRam.s1Dim0->push_back(element);
                    } else {
                        Rtree::SaveEntryWithOrderIndex(element, splitBuffers.split1Dim0File.value());
                    }
                    if (((currentXSplit1 + 1) % SSplit1 == 0 && (currentXSplit1 + 1) / SSplit1 >= 1 && (currentXSplit1 + 1) / SSplit1 < M)
                        || (currentXSplit1 % SSplit1 == 0 && currentXSplit1 / SSplit1 >= 1 && currentXSplit1 / SSplit1 < M)) {
                        // index i * S - 1 or i * S
                        splitBuffers.splitBuffersRam.s1SmallDim0->push_back(element);
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

                    if (splitResult.bestDim == 1) {
                        if (currentXSplit1 == 0) {
                            minSplit1OtherDim = element;
                            maxSplit1OtherDim = element;
                        }
                        if (element.second.first > maxSplit1OtherDim.second.first) {
                            maxSplit1OtherDim = element;
                        }
                    }

                    currentXSplit1++;
                } else {
                    if (split1InRam || this->workInRam) {
                        splitBuffers.splitBuffersRam.s1Dim1->push_back(element);
                    } else {
                        Rtree::SaveEntryWithOrderIndex(element, splitBuffers.split1Dim1File.value());
                    }
                    if (((currentYSplit1 + 1) % SSplit1 == 0 && (currentYSplit1 + 1) / SSplit1 >= 1 && (currentYSplit1 + 1) / SSplit1 < M)
                        || (currentYSplit1 % SSplit1 == 0 && currentYSplit1 / SSplit1 >= 1 && currentYSplit1 / SSplit1 < M)) {
                        // index i * S - 1 or i * S
                        splitBuffers.splitBuffersRam.s1SmallDim1->push_back(element);
                    }

                    if (splitResult.bestDim == 0) {
                        if (currentYSplit1 == 0) {
                            minSplit1OtherDim = element;
                            maxSplit1OtherDim = element;
                        }
                        if (element.second.first > maxSplit1OtherDim.second.first) {
                            maxSplit1OtherDim = element;
                        }
                    }

                    currentYSplit1++;
                }
            }
            i++;

            if (!this->workInRam) {
                if (dim == 0)
                    elementOpt = fileReaderDim0.value().GetNextElement();
                if (dim == 1)
                    elementOpt = fileReaderDim1.value().GetNextElement();
            }
        }
    }

    if (!this->workInRam) {
        fileReaderDim0.value().Close();
        fileReaderDim1.value().Close();
    }

    // replace the placeholder
    if (splitResult.bestDim == 0) {
        (*splitBuffers.splitBuffersRam.s0SmallDim1)[0] = minSplit0OtherDim;
        (*splitBuffers.splitBuffersRam.s0SmallDim1)[1] = maxSplit0OtherDim;
        (*splitBuffers.splitBuffersRam.s1SmallDim1)[0] = minSplit1OtherDim;
        (*splitBuffers.splitBuffersRam.s1SmallDim1)[1] = maxSplit1OtherDim;
    } else {
        (*splitBuffers.splitBuffersRam.s0SmallDim0)[0] = minSplit0OtherDim;
        (*splitBuffers.splitBuffersRam.s0SmallDim0)[1] = maxSplit0OtherDim;
        (*splitBuffers.splitBuffersRam.s1SmallDim0)[0] = minSplit1OtherDim;
        (*splitBuffers.splitBuffersRam.s1SmallDim0)[1] = maxSplit1OtherDim;
    }

    boxGeo boxSplit0 = Rtree::createBoundingBox(globalMinXS0, globalMinYS0, globalMaxXS0, globalMaxYS0);
    boxGeo boxSplit1 = Rtree::createBoundingBox(globalMinXS1, globalMinYS1, globalMaxXS1, globalMaxYS1);

    return std::make_pair(boxSplit0, boxSplit1);
}

FileReader::FileReader(const std::string& filePath) {
    this->filePath = filePath;

    this->file = std::ifstream(this->filePath, std::ios::binary);
    this->file.seekg (0, std::ifstream::end);
    this->fileLength = this->file.tellg();
    this->file.seekg (0, std::ifstream::beg);
}

std::optional<rTreeValueWithOrderIndex> FileReader::GetNextElement() {
    if (this->file.tellg() < this->fileLength) {
        double minX;
        double minY;
        double maxX;
        double maxY;
        uint64_t id;
        long long orderX;
        long long orderY;

        this->file.read(reinterpret_cast<char*>(&minX), sizeof(double));
        this->file.read(reinterpret_cast<char*>(&minY), sizeof(double));
        this->file.read(reinterpret_cast<char*>(&maxX), sizeof(double));
        this->file.read(reinterpret_cast<char*>(&maxY), sizeof(double));
        this->file.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));
        this->file.read(reinterpret_cast<char*>(&orderX), sizeof(long long));
        this->file.read(reinterpret_cast<char*>(&orderY), sizeof(long long));

        boxGeo box = Rtree::createBoundingBox(minX, minY, maxX, maxY);
        rTreeValue boxWithId = std::make_pair(box, id);
        rTreeValueWithOrderIndex element = std::make_pair(boxWithId, std::make_pair(orderX, orderY));

        return { element };
    } else {
        return {};
    }
}

void FileReader::Close() {
    this->file.close();
}
