//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include "Rtree/Rtree.h"

static void centerOrderingExt(multiBoxGeo& boxes, size_t dim) {
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

static void centerOrderingExt(multiBoxWithOrderIndex& boxes, size_t dim) {
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

std::pair<boxGeo, long long> ExternalSort(const std::string& onDiskBase, std::shared_ptr<multiBoxWithOrderIndex>& r0Small, std::shared_ptr<multiBoxWithOrderIndex>& r1Small, size_t M) {
    multiBoxGeo RectanglesD0 = Rtree::LoadEntries(onDiskBase + ".boundingbox.tmp");

    centerOrderingExt(RectanglesD0, 0);

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

    multiBoxWithOrderIndex RectanglesD1 = Rtree::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d0.tmp");
    centerOrderingExt(RectanglesD1, 1);

    size_t currentS = std::ceil(((float) xSize) / ((float) M));

    long long ySize = 0;
    std::ofstream r1File = std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
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

    multiBoxWithOrderIndex RectanglesD0Second = Rtree::LoadEntriesWithOrderIndex(onDiskBase + ".boundingbox.d1.tmp");
    centerOrderingExt(RectanglesD0Second, 0);

    long long currentX = 0;
    std::ofstream r0FileSecond = std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
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
    return std::make_pair(boundingBox, xSize);
}
