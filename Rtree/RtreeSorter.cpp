//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#include "./Rtree.h"
//#include <util/BackgroundStxxlSorter.h>
#include "../ExternalSorting.cpp"

template <size_t dimension>
struct SortRuleLambda {
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
struct SortRuleLambdaWithIndex {
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

static void centerOrdering(multiBoxGeo& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        std::sort(boxes.begin(), boxes.end(), SortRuleLambda<0>{});
    } else {
        // order by centerY
        std::sort(boxes.begin(), boxes.end(), SortRuleLambda<1>{});
    }
}

static void centerOrdering(multiBoxWithOrderIndex& boxes, size_t dim) {
    if (dim == 0) {
        // order by centerX
        std::sort(boxes.begin(), boxes.end(), SortRuleLambdaWithIndex<0>{});
    } else {
        // order by centerY
        std::sort(boxes.begin(), boxes.end(), SortRuleLambdaWithIndex<1>{});
    }
}

/*OrderedBoxes SortInput(const std::filesystem::path& onDiskBase, size_t M,
                       uintmax_t maxBuildingRamUsage, bool workInRam) {
    OrderedBoxes orderedInputRectangles;

    auto maxRamForSorter =
            std::ceil((maxBuildingRamUsage < 9999999999.0 ? maxBuildingRamUsage
                                                          : 9999999999.0) /
                      3.0);
    ad_utility::BackgroundStxxlSorter<RTreeValue, SortRuleLambda<0>>
                                                                  sorterRectsD0Basic =
                                                          ad_utility::BackgroundStxxlSorter<RTreeValue, SortRuleLambda<0>>(
            maxRamForSorter);
    multiBoxGeo rectsD0Basic;

    if (workInRam) {
        rectsD0Basic = Rtree::LoadEntries(onDiskBase + ".boundingbox.tmp");
        centerOrdering(rectsD0Basic, 0);
    } else {
        for (const RTreeValue& rectD0Element :
                FileReaderWithoutIndex(onDiskBase + ".boundingbox.tmp")) {
            sorterRectsD0Basic.push(rectD0Element);
        }
    }

    uint64_t xSize = 0;
    Rtree::BoundingBox boundingBox = Rtree::createBoundingBox(0, 0, 0, 0);

    ad_utility::BackgroundStxxlSorter<RTreeValueWithOrderIndex,
            SortRuleLambdaWithIndex<1>>
                                     sorterRectsD1 =
                    ad_utility::BackgroundStxxlSorter<RTreeValueWithOrderIndex,
            SortRuleLambdaWithIndex<1>>(
                    maxRamForSorter);
    std::shared_ptr<multiBoxWithOrderIndex> RectanglesD1WithOrder =
            std::make_shared<multiBoxWithOrderIndex>();

    if (workInRam) {
        for (RTreeValue element : rectsD0Basic) {
            RTreeValueWithOrderIndex entry = {{element.box, element.id}, xSize, 0};
            RectanglesD1WithOrder->push_back(entry);
            xSize++;

            boundingBox = Rtree::combineBoundingBoxes(boundingBox, element.box);
        }
        centerOrdering(*RectanglesD1WithOrder, 1);
    } else {
        for (RTreeValue element : sorterRectsD0Basic.sortedView()) {
            RTreeValueWithOrderIndex entry = {{element.box, element.id}, xSize, 0};
            sorterRectsD1.push(entry);
            xSize++;

            boundingBox = Rtree::combineBoundingBoxes(boundingBox, element.box);
        }
    }
    sorterRectsD0Basic.clear();

    size_t currentS = std::ceil(((float)xSize) / ((float)M));
    if (xSize <= M * M) {
        // in this case S can just be M
        currentS = M;
    }

    uint64_t ySize = 0;
    std::ofstream r1File =
            std::ofstream(onDiskBase + ".boundingbox.d1.tmp", std::ios::binary);
    ad_utility::BackgroundStxxlSorter<RTreeValueWithOrderIndex,
            SortRuleLambdaWithIndex<0>>
                                     sorterRectsD0 =
                    ad_utility::BackgroundStxxlSorter<RTreeValueWithOrderIndex,
            SortRuleLambdaWithIndex<0>>(
                    maxRamForSorter);
    std::shared_ptr<multiBoxWithOrderIndex> RectanglesD0WithOrder =
            std::make_shared<multiBoxWithOrderIndex>();
    std::shared_ptr<multiBoxWithOrderIndex> r1Small =
            std::make_shared<multiBoxWithOrderIndex>();
    // placeholder
    r1Small->push_back(RTreeValueWithOrderIndex());
    r1Small->push_back(RTreeValueWithOrderIndex());
    RTreeValueWithOrderIndex minD1;
    RTreeValueWithOrderIndex maxD1;

    auto processD1Element = [&ySize, currentS, M, &r1Small, &minD1,
            &maxD1](RTreeValueWithOrderIndex& element) {
        element.orderY = ySize;

        if (isBorderOfSplitCandidate(ySize, currentS, M)) {
            // index i * S - 1 or i * S
            r1Small->push_back(element);
        }

        if (ySize == 0) {
            minD1 = element;
        }

        maxD1 = element;

        ySize++;
    };

    if (workInRam) {
        for (RTreeValueWithOrderIndex element : *RectanglesD1WithOrder) {
            processD1Element(element);

            RectanglesD0WithOrder->push_back(element);
        }
        centerOrdering(*RectanglesD0WithOrder, 0);
    } else {
        for (RTreeValueWithOrderIndex element : sorterRectsD1.sortedView()) {
            processD1Element(element);

            Rtree::SaveEntryWithOrderIndex(element, r1File);
            sorterRectsD0.push(element);
        }
    }

    r1File.close();
    sorterRectsD1.clear();

    // replace the placeholder
    (*r1Small)[0] = minD1;
    (*r1Small)[1] = maxD1;

    uint64_t currentX = 0;
    std::ofstream r0File =
            std::ofstream(onDiskBase + ".boundingbox.d0.tmp", std::ios::binary);
    std::shared_ptr<multiBoxWithOrderIndex> r0Small =
            std::make_shared<multiBoxWithOrderIndex>();
    // placeholder
    r0Small->push_back(RTreeValueWithOrderIndex());
    r0Small->push_back(RTreeValueWithOrderIndex());
    RTreeValueWithOrderIndex minD0;
    RTreeValueWithOrderIndex maxD0;

    auto processD0Element = [&currentX, currentS, M, &r0Small, &minD0,
            &maxD0](RTreeValueWithOrderIndex& element) {
        if (isBorderOfSplitCandidate(currentX, currentS, M)) {
            // index i * S - 1 or i * S
            r0Small->push_back(element);
        }

        if (currentX == 0) {
            minD0 = element;
        }
        maxD0 = element;

        currentX++;
    };

    if (workInRam) {
        for (RTreeValueWithOrderIndex element : *RectanglesD0WithOrder) {
            processD0Element(element);
        }
    } else {
        for (RTreeValueWithOrderIndex element : sorterRectsD0.sortedView()) {
            Rtree::SaveEntryWithOrderIndex(element, r0File);

            processD0Element(element);
        }
    }

    r0File.close();
    sorterRectsD0.clear();

    // replace the placeholder
    (*r0Small)[0] = minD0;
    (*r0Small)[1] = maxD0;

    RectanglesForOrderedBoxes rectsD0;
    RectanglesForOrderedBoxes rectsD1;
    rectsD0.rectanglesSmall = r0Small;
    rectsD1.rectanglesSmall = r1Small;
    if (workInRam) {
        rectsD0.rectanglesInRam = RectanglesD0WithOrder;
        rectsD1.rectanglesInRam = RectanglesD1WithOrder;
        orderedInputRectangles.SetOrderedBoxesToRam(rectsD0, rectsD1, boundingBox);
    } else {
        rectsD0.rectanglesOnDisk = onDiskBase + ".boundingbox.d0";
        rectsD1.rectanglesOnDisk = onDiskBase + ".boundingbox.d1";
        orderedInputRectangles.SetOrderedBoxesToDisk(rectsD0, rectsD1, xSize,
                                                     boundingBox);
    }
    return orderedInputRectangles;
}*/

OrderedBoxes SortInput(const std::filesystem::path& onDiskBase, size_t M, uintmax_t
maxBuildingRamUsage, bool workInRam) { if (workInRam) { return
                InternalSort(onDiskBase, M); } else { return ExternalSort(onDiskBase, M,
                                                                          maxBuildingRamUsage);
    }
}
