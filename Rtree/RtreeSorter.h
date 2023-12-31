//  Copyright 2023, University of Freiburg,
//                  Chair of Algorithms and Data Structures.
//  Author: Noah Nock <noah.v.nock@gmail.com>

#ifndef BACHELORPROJEKT_RTREESORTER_H
#define BACHELORPROJEKT_RTREESORTER_H

#include "./Rtree.h"

OrderedBoxes SortInput(const std::string& onDiskBase,
                       const std::string& fileSuffix, size_t M,
                       uintmax_t maxBuildingRamUsage, bool workInRam);

#endif  // BACHELORPROJEKT_RTREESORTER_H