//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include "Rtree/Rtree.h"

static void centerOrderingExt(multiBoxGeo& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        sortRuleLambdaX comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    } else {
        // order by centerY
        auto sortRuleLambda = [](rTreeValue b1, rTreeValue b2) -> bool {
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
        sortRuleLambdaXWithIndex comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    } else {
        // order by centerY
        sortRuleLambdaYWithIndex comp;

        std::sort(boxes.begin(), boxes.end(), comp);
    }
}

OrderedBoxes ExternalSort(const std::string& onDiskBase, size_t M, uintmax_t maxBuildingRamUsage) {
    OrderedBoxes orderedInputRectangles;
    multiBoxGeo RectanglesD0 = Rtree::LoadEntries(onDiskBase + ".boundingbox.tmp");

    centerOrderingExt(RectanglesD0, 0);

    long long xSize = 0;
    double globalMinX = -1;
    double globalMinY = -1;
    double globalMaxX = -1;
    double globalMaxY = -1;

    std::ofstream r0File = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
    for (rTreeValue element : RectanglesD0) {
        rTreeValueWithOrderIndex entry = rTreeValueWithOrderIndex(element.box, element.id, xSize, 0);
        Rtree::SaveEntryWithOrderIndex(entry, r0File);
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

    multiBoxWithOrderIndex RectanglesD1 = Rtree::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d0.tmp");
    centerOrderingExt(RectanglesD1, 1);

    size_t currentS = std::ceil(((float) xSize) / ((float) M));

    long long ySize = 0;
    std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
    std::shared_ptr<multiBoxWithOrderIndex> r1Small = std::make_shared<multiBoxWithOrderIndex>();
    r1Small->push_back(RectanglesD1[0]);
    rTreeValueWithOrderIndex maxElementDim1 = RectanglesD1[RectanglesD1.size() - 1];
    maxElementDim1.orderY = RectanglesD1.size() - 1;
    r1Small->push_back(maxElementDim1);
    for (rTreeValueWithOrderIndex element : RectanglesD1) {
        element.orderY = ySize;
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

    multiBoxWithOrderIndex RectanglesD0Second = Rtree::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d1.tmp");
    centerOrderingExt(RectanglesD0Second, 0);

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

    boxGeo boundingBox = Rtree::createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
    orderedInputRectangles.CreateOrderedBoxesOnDisk(onDiskBase + ".boundingbox.d0", onDiskBase + ".boundingbox.d1", r0Small, r1Small, xSize, boundingBox);
    return orderedInputRectangles;
}

OrderedBoxes InternalSort(const std::string& onDiskBase, size_t M) {
    OrderedBoxes orderedInputRectangles;
    multiBoxGeo RectanglesD0 = Rtree::LoadEntries(onDiskBase + ".boundingbox.tmp");
    centerOrderingExt(RectanglesD0, 0);

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
        rTreeValueWithOrderIndex entry = rTreeValueWithOrderIndex(element.box, element.id, i, 0);
        RectanglesD1WithOrder->push_back(entry);

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

    centerOrderingExt(*RectanglesD1WithOrder, 1);

    R1Small->push_back((*RectanglesD1WithOrder)[0]);
    rTreeValueWithOrderIndex maxElementDim1 = (*RectanglesD1WithOrder)[RectanglesD1WithOrder->size() - 1];
    maxElementDim1.orderY = RectanglesD1WithOrder->size() - 1;
    R1Small->push_back(maxElementDim1);
    for (long long i = 0; i < RectanglesD1WithOrder->size(); i++) {
        (*RectanglesD1WithOrder)[i].orderY = i;

        if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
            || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
            // index i * S - 1 or i * S
            R1Small->push_back((*RectanglesD1WithOrder)[i]);
        }
    }

    std::shared_ptr<multiBoxWithOrderIndex> RectanglesD0WithOrder = std::make_shared<multiBoxWithOrderIndex>(*RectanglesD1WithOrder);
    centerOrderingExt(*RectanglesD0WithOrder, 0);

    R0Small->push_back((*RectanglesD0WithOrder)[0]);
    rTreeValueWithOrderIndex maxElementDim0 = (*RectanglesD0WithOrder)[RectanglesD0WithOrder->size() - 1];
    maxElementDim0.orderY = RectanglesD0WithOrder->size() - 1;
    R0Small->push_back(maxElementDim0);
    for (long long i = 0; i < RectanglesD0WithOrder->size(); i++) {
        if (((i + 1) % currentS == 0 && (i + 1) / currentS >= 1 && (i + 1) / currentS < M)
            || (i % currentS == 0 && i / currentS >= 1 && i / currentS < M)) {
            // index i * S - 1 or i * S
            R0Small->push_back((*RectanglesD0WithOrder)[i]);
        }
    }

    boxGeo boundingBox = Rtree::createBoundingBox(globalMinX, globalMinY, globalMaxX, globalMaxY);
    orderedInputRectangles.CreateOrderedBoxesInRam(RectanglesD0WithOrder, RectanglesD1WithOrder, R0Small, R1Small, boundingBox);
    return orderedInputRectangles;
}
