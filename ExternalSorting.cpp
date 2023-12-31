//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include "Rtree/Rtree.h"
#include "Rtree/RtreeFileReader.h"

template <size_t dimension>
struct SortRuleLambda2 {
    // comparison function
    bool operator()(const RTreeValue& b1, const RTreeValue& b2) const {
        double center1 = dimension == 0 ? std::midpoint(b1.MinX(), b1.MaxX())
                                        : std::midpoint(b1.MinY(), b1.MaxY());
        double center2 = dimension == 0 ? std::midpoint(b2.MinX(), b2.MaxX())
                                        : std::midpoint(b2.MinY(), b2.MaxY());
        return center1 < center2;
    }

    // Value that is strictly smaller than any input element.
    static RTreeValue min_value() {
        return {BasicGeometry::CreateBoundingBox(-std::numeric_limits<double>::max(),
                                                 -std::numeric_limits<double>::max(),
                                                 -std::numeric_limits<double>::max(),
                                                 -std::numeric_limits<double>::max()),
                0};
    }

    // Value that is strictly larger than any input element.
    static RTreeValue max_value() {
        return {BasicGeometry::CreateBoundingBox(std::numeric_limits<double>::max(),
                                                 std::numeric_limits<double>::max(),
                                                 std::numeric_limits<double>::max(),
                                                 std::numeric_limits<double>::max()),
                0};
    }
};

template <size_t dimension>
struct SortRuleLambdaWithIndex2 {
    uint64_t RTreeValueWithOrderIndex::*orderSelected =
            dimension == 0 ? &RTreeValueWithOrderIndex::orderX
                           : &RTreeValueWithOrderIndex::orderY;

    // comparison function
    bool operator()(const RTreeValueWithOrderIndex& b1,
                    const RTreeValueWithOrderIndex& b2) const {
        double center1 = dimension == 0 ? std::midpoint(b1.MinX(), b1.MaxX())
                                        : std::midpoint(b1.MinY(), b1.MaxY());
        double center2 = dimension == 0 ? std::midpoint(b2.MinX(), b2.MaxX())
                                        : std::midpoint(b2.MinY(), b2.MaxY());

        if (b1.*orderSelected == b2.*orderSelected) return center1 < center2;
        return b1.*orderSelected < b2.*orderSelected;
    }

    // Value that is strictly smaller than any input element.
    static RTreeValueWithOrderIndex min_value() {
        return {{BasicGeometry::CreateBoundingBox(-std::numeric_limits<double>::max(),
                                                  -std::numeric_limits<double>::max(),
                                                  -std::numeric_limits<double>::max(),
                                                  -std::numeric_limits<double>::max()),
                 0},
                0,
                0};
    }

    // Value that is strictly larger than any input element.
    static RTreeValueWithOrderIndex max_value() {
        return {{BasicGeometry::CreateBoundingBox(std::numeric_limits<double>::max(),
                                                  std::numeric_limits<double>::max(),
                                                  std::numeric_limits<double>::max(),
                                                  std::numeric_limits<double>::max()),
                 0},
                std::numeric_limits<long long>::max(),
                std::numeric_limits<long long>::max()};
    }
};

static void centerOrderingExt(multiBoxGeo& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        SortRuleLambda2<0> comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    } else {
        // order by centerY
        auto sortRuleLambda = [](RTreeValue b1, RTreeValue b2) -> bool {
            double center1 = (b1.box.min_corner().get<1>() + b1.box.max_corner().get<1>()) / 2;
            double center2 = (b2.box.min_corner().get<1>() + b2.box.max_corner().get<1>()) / 2;
            return center1 < center2;
        };

        std::sort(boxes.begin(), boxes.end(), sortRuleLambda);
    }
}

static void centerOrderingExt(multiBoxWithOrderIndex& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        SortRuleLambdaWithIndex2<0> comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    } else {
        // order by centerY
        SortRuleLambdaWithIndex2<1> comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    }
}

OrderedBoxes ExternalSort(const std::string& onDiskBase, size_t M, uintmax_t maxBuildingRamUsage) {
    OrderedBoxes orderedInputRectangles;
    multiBoxGeo RectanglesD0 = FileReaderWithoutIndex::LoadEntries(onDiskBase + ".boundingbox.tmp");

    centerOrderingExt(RectanglesD0, 0);

    uint64_t xSize = 0;
    double globalMinX = -1;
    double globalMinY = -1;
    double globalMaxX = -1;
    double globalMaxY = -1;

    std::ofstream r0File = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
    for (RTreeValue element : RectanglesD0) {
        RTreeValueWithOrderIndex entry = {{element.box, element.id}, xSize, 0};
        FileReader::SaveEntryWithOrderIndex(entry, r0File);
        xSize++;

        if (globalMinX == -1 || element.box.min_corner().get<0>() < globalMinX) {
            globalMinX = element.box.min_corner().get<0>();
        }
        if (globalMinY == -1 || element.box.min_corner().get<1>() < globalMinY) {
            globalMinY = element.box.min_corner().get<1>();
        }
        if (element.box.max_corner().get<0>() > globalMaxX) {
            globalMaxX = element.box.max_corner().get<0>();
        }
        if (element.box.max_corner().get<1>() > globalMaxY) {
            globalMaxY = element.box.max_corner().get<1>();
        }
    }
    r0File.close();
    RectanglesD0.clear();

    multiBoxWithOrderIndex RectanglesD1 = FileReader::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d0.tmp");
    centerOrderingExt(RectanglesD1, 1);

    size_t currentS = std::ceil(((float) xSize) / ((float) M));

    uint64_t ySize = 0;
    std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
    multiBoxWithOrderIndex r1Small = multiBoxWithOrderIndex();
    r1Small.push_back(RectanglesD1[0]);
    RTreeValueWithOrderIndex maxElementDim1 = RectanglesD1[RectanglesD1.size() - 1];
    maxElementDim1.orderY = RectanglesD1.size() - 1;
    r1Small.push_back(maxElementDim1);
    for (RTreeValueWithOrderIndex element : RectanglesD1) {
        element.orderY = ySize;
        FileReader::SaveEntryWithOrderIndex(element, r1File);

        if (((ySize + 1) % currentS == 0 && (ySize + 1) / currentS >= 1 && (ySize + 1) / currentS < M)
            || (ySize % currentS == 0 && ySize / currentS >= 1 && ySize / currentS < M)) {
            // index i * S - 1 or i * S
            r1Small.push_back(element);
        }

        ySize++;
    }
    r1File.close();
    RectanglesD1.clear();

    multiBoxWithOrderIndex RectanglesD0Second = FileReader::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d1.tmp");
    centerOrderingExt(RectanglesD0Second, 0);

    uint64_t currentX = 0;
    std::ofstream r0FileSecond = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
    multiBoxWithOrderIndex r0Small = multiBoxWithOrderIndex();
    r0Small.push_back(RectanglesD0Second[0]);
    r0Small.push_back(RectanglesD0Second[RectanglesD0Second.size() - 1]);
    for (RTreeValueWithOrderIndex element : RectanglesD0Second) {
        FileReader::SaveEntryWithOrderIndex(element, r0FileSecond);

        if (((currentX + 1) % currentS == 0 && (currentX + 1) / currentS >= 1 && (currentX + 1) / currentS < M)
            || (currentX % currentS == 0 && currentX / currentS >= 1 && currentX / currentS < M)) {
            // index i * S - 1 or i * S
            r0Small.push_back(element);
        }

        currentX++;
    }
    r0FileSecond.close();
    RectanglesD0Second.clear();

    BasicGeometry::BoundingBox boundingBox = BasicGeometry::CreateBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
    RectanglesForOrderedBoxes d0WithOrder;
    d0WithOrder.rectangles = onDiskBase + ".boundingbox.d0.tmp";
    d0WithOrder.rectanglesSmall = r0Small;
    RectanglesForOrderedBoxes d1WithOrder;
    d1WithOrder.rectangles = onDiskBase + ".boundingbox.d1.tmp";
    d1WithOrder.rectanglesSmall = r1Small;
    orderedInputRectangles.SetOrderedBoxesToDisk(d0WithOrder, d1WithOrder, xSize, boundingBox);
    return orderedInputRectangles;
}

OrderedBoxes InternalSort(const std::string& onDiskBase, size_t M) {
    OrderedBoxes orderedInputRectangles;
    multiBoxGeo RectanglesD0 = FileReaderWithoutIndex::LoadEntries(onDiskBase + ".boundingbox.tmp");
    centerOrderingExt(RectanglesD0, 0);

    double globalMinX = -1;
    double globalMinY = -1;
    double globalMaxX = -1;
    double globalMaxY = -1;

    size_t currentS = std::ceil(((float) RectanglesD0.size()) / ((float) M));

    multiBoxWithOrderIndex R0Small = multiBoxWithOrderIndex();
    multiBoxWithOrderIndex R1Small = multiBoxWithOrderIndex();

    multiBoxWithOrderIndex RectanglesD1WithOrder = multiBoxWithOrderIndex();
    for (uint64_t i = 0; i < RectanglesD0.size(); i++) {
        RTreeValue element = RectanglesD0[i];
        RTreeValueWithOrderIndex entry = {{element.box, element.id}, i, 0};
        RectanglesD1WithOrder.push_back(entry);

        if (globalMinX == -1 || element.box.min_corner().get<0>() < globalMinX) {
            globalMinX = element.box.min_corner().get<0>();
        }
        if (globalMinY == -1 || element.box.min_corner().get<1>() < globalMinY) {
            globalMinY = element.box.min_corner().get<1>();
        }
        if (element.box.max_corner().get<0>() > globalMaxX) {
            globalMaxX = element.box.max_corner().get<0>();
        }
        if (element.box.max_corner().get<1>() > globalMaxY) {
            globalMaxY = element.box.max_corner().get<1>();
        }
    }

    centerOrderingExt(RectanglesD1WithOrder, 1);

    R1Small.push_back((RectanglesD1WithOrder)[0]);
    RTreeValueWithOrderIndex maxElementDim1 = (RectanglesD1WithOrder)[RectanglesD1WithOrder.size() - 1];
    maxElementDim1.orderY = RectanglesD1WithOrder.size() - 1;
    R1Small.push_back(maxElementDim1);
    for (uint64_t i = 0; i < RectanglesD1WithOrder.size(); i++) {
        (RectanglesD1WithOrder)[i].orderY = i;

        if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
            || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
            // index i * S - 1 or i * S
            R1Small.push_back((RectanglesD1WithOrder)[i]);
        }
    }

    multiBoxWithOrderIndex RectanglesD0WithOrder = multiBoxWithOrderIndex(RectanglesD1WithOrder);
    centerOrderingExt(RectanglesD0WithOrder, 0);

    R0Small.push_back((RectanglesD0WithOrder)[0]);
    RTreeValueWithOrderIndex maxElementDim0 = (RectanglesD0WithOrder)[RectanglesD0WithOrder.size() - 1];
    maxElementDim0.orderY = RectanglesD0WithOrder.size() - 1;
    R0Small.push_back(maxElementDim0);
    for (uint64_t i = 0; i < RectanglesD0WithOrder.size(); i++) {
        if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
            || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
            // index i * S - 1 or i * S
            R0Small.push_back((RectanglesD0WithOrder)[i]);
        }
    }

    BasicGeometry::BoundingBox boundingBox = BasicGeometry::CreateBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
    RectanglesForOrderedBoxes d0WithOrder;
    d0WithOrder.rectangles = RectanglesD0WithOrder;
    d0WithOrder.rectanglesSmall = R0Small;
    RectanglesForOrderedBoxes d1WithOrder;
    d1WithOrder.rectangles = RectanglesD1WithOrder;
    d1WithOrder.rectanglesSmall = R1Small;
    orderedInputRectangles.SetOrderedBoxesToRam(d0WithOrder, d1WithOrder, boundingBox);
    return orderedInputRectangles;
}
