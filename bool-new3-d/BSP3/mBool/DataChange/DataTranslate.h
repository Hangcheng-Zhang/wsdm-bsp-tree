#pragma once

class DataTranslate {

public:
	static Triangle_Set TranslateMesh(Mesh& p_Mesh);
	static void TranslateMeshWithlable(Mesh& p_Mesh, unordered_map<FaceId, bool> p_classify, pair<Triangle_Set, Triangle_Set>& r_labeledTriangles);


	static Edge_Set TranslateEdge(Mesh& p_Mesh);

	static Vertex_Set TranslateVertex(Mesh& p_Mesh);

	static pair<Triangle_Set,Edge_Set> TranslatePartitionPlane(bspFaceRenderingInfo& partiPlanes, Mesh& r_PartionPlaneMesh);

	static Triangle_Set TranslateTriangleBSPs(BSPTreeNode* leafNode);
	static pair<Triangle_Set, Triangle_Set> TranslateTriangleBSPsFromTwoMeshs(BSPTreeNode* leafNode);

	static pair<Mesh, Mesh> TranslateNode2SeperatedTriangleMesh(BSPTreeNode* leafNode);

	static Mesh TranslateFaceID2Mesh(Mesh& m_OriginMesh, unordered_set<FaceId>& facesId, bool needReverseNormal = false);
};