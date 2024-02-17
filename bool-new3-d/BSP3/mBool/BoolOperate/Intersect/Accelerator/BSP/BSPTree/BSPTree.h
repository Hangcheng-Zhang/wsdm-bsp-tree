#pragma once


class BSPTree {


public:
	BSPTree(Mesh& p_meshA, Mesh& p_meshB, Tolerance& p_Tolerance, int p_MaxDepth, int p_LeafShresholdNum, BSPConstructType p_Type = BSPConstructType::AABB_MIDDLE_SPLIT);

protected:
	bool Initialize();

protected:
	vector<TriangleBSP> m_MeshA;
	vector<TriangleBSP> m_MeshB;

	BSPConstructType m_Type;
	Tolerance& m_Tolerance;

	BSPTreeNode* m_Header;

public:
	void GetPartitionPlane(bspFaceRenderingInfo& r_Pplanes);


public:
	void GetAllLeafNode(vector<BSPTreeNode*>& r_AllLeafNodes);
};