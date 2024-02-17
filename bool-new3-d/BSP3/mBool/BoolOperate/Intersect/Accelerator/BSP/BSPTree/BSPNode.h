#pragma once


class BSPTreeNode {

public:
	BSPTreeNode(vector<TriangleBSP> p_FacesA, vector<TriangleBSP> p_FacesB, Tolerance& p_Tolerance, int p_Depth, int p_MaxDepth, int p_LeafShresholdNum, bool p_Isroot = false, BSPConstructType p_Type = BSPConstructType::AABB_MIDDLE_SPLIT);
	//BSPTreeNode(vector<TriangleBSP> p_FacesA, vector<TriangleBSP> p_FacesB, Plane3d p_PartitionPlane, Tolerance& p_Tolerance, bool p_Isroot = false, BSPConstructType p_Type = BSPConstructType::AABB);

public:
	Plane3d m_PartitionPlane;

	BSPTreeNode* m_PositiveSon;
	BSPTreeNode* m_NegativeSon;

	bool m_IsRoot;
	bool m_IsLeaf;
	int m_Depth;

	int m_MaxDepth;
	int m_LeafShresholdNum;


	BSPConstructType m_Type;

	Tolerance& m_Tolerance;

	vector<TriangleBSP> m_FacesA;
	vector<TriangleBSP> m_FacesB;
public:
	void GetLeafNodes(vector<BSPTreeNode*>& r_AllLeafNodes);


	bool IsLeaf();
	void GetIntersectPair(IntersectTriangleCheckList& p_IntersectCheckList);

	bool IsValidNode();
	void GetValidFaces(unordered_set<FaceId>& facesA, unordered_set<FaceId>& facesB);


	void GetPartitionPlane(bspFaceRenderingInfo& r_Pplanes, planeRecords p_CurrentRecords);

	void ConstructAABBSons();
	void ConstructSDMSons();
	void ConstructSAHSons();
	void ConstructObbMiddleSons();
	void ConstructGravity_SPLITSons();
	void ConstructSDM_OBBSons();

	void GetNegativePositiveLabelWithSplitTriangle(vector<TriangleBSP>& negativeFaces, vector<TriangleBSP>& positiveFaces, Plane3d partitionPlane);
	void GetNegativePositiveLabelWithOutSplitTriangle(vector<TriangleBSP>& negativeFaces, vector<TriangleBSP>& positiveFaces, Plane3d partitionPlane);


	void SeparationTrianglesFromDifferentMesh(vector<TriangleBSP> facesOrigin, vector<TriangleBSP>& facesA, vector<TriangleBSP>& facesB);

	int Size();
	pair<int,int> GetRealSize();
};