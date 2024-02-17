#pragma once

enum class BSPConstructType
{
	AABB_MIDDLE_SPLIT, SDM, SAH, ObbMiddel, Gravity_SPLIT, SDM_OBB
};

#define MAX_TREE_DEPTH 10
#define LEAF_SHRESHOLD_NUM 25
#define TestType BSPConstructType::Gravity_SPLIT

using planeRecord = pair<Plane3d, bool>;

using planeRecords = vector<planeRecord>;
using bspFaceRenderingInfo = vector<pair<Plane3d, planeRecords>>;

#include "BSPNode.h"
#include "BSPTree.h"

