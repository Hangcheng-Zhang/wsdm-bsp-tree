#pragma once


namespace TestFunction {

	void TestBool();

	void TestWrite();

	void TestReadObjWriteStl();

	void TestOpenMesh();

	void TestMath();

	void TestBspTree();

	void SingleTestClassicModelsForBSP();
	void SingleTestClassicModelsForBSP(int ID);

	void SingleTestThingi10K();

	void BatchTestClassicModels();

	void BatchTestThingi10K();


	pair<double, double> TestBspTimeConsum(Mesh& meshA, Mesh& meshB, Tolerance& toler, int maxDepth, int leafShreshold, string p_BaseAddress, bool outPutResultMesh, bool outPutITRA, BSPConstructType bspType = BSPConstructType::AABB_MIDDLE_SPLIT, int Id = -1);

	pair<double, double> TestBspNodeDistributionInfo(Mesh& meshA, Mesh& meshB, Tolerance& toler, int maxDepth, int leafShreshold, string p_BaseAddress, bool outPutResultMesh, BSPConstructType bspType = BSPConstructType::AABB_MIDDLE_SPLIT);
}