#include "mPch.h"
#include "BSPTreeHeader.h"

BSPTreeNode::BSPTreeNode(vector<TriangleBSP> p_FacesA, vector<TriangleBSP> p_FacesB, Tolerance& p_Tolerance, int p_Depth, int p_MaxDepth, int p_LeafShresholdNum, bool p_Isroot /*= false*/, BSPConstructType p_Type /*= BSPConstructType::AABB*/) :
	m_FacesA(p_FacesA), m_FacesB(p_FacesB), m_Type(p_Type), m_Tolerance(p_Tolerance), m_Depth(p_Depth), m_IsRoot(p_Isroot), m_MaxDepth(p_MaxDepth), m_LeafShresholdNum(p_LeafShresholdNum)
{
	//达到截止深度
	if (p_Depth == p_MaxDepth) {
		m_IsLeaf = true;
		m_PositiveSon = nullptr;
		m_NegativeSon = nullptr;
		return;
	}

	//达到叶节点面片数目条件，生成叶节点
	if (p_FacesA.size() + p_FacesB.size() <= p_LeafShresholdNum) {
		m_IsLeaf = true;
		m_PositiveSon = nullptr;
		m_NegativeSon = nullptr;

		return;
	}

	//某一类面片数目为0，生成叶节点
	if ((0 == p_FacesA.size()) || (0 == p_FacesB.size())) {

		m_IsLeaf = true;
		m_PositiveSon = nullptr;
		m_NegativeSon = nullptr;

		return;
	}


	switch (m_Type) {
	case BSPConstructType::AABB_MIDDLE_SPLIT:
		ConstructAABBSons();
		break;

	case BSPConstructType::SDM:
		ConstructSDMSons();
		break;

	case BSPConstructType::SAH:
		ConstructSAHSons();
		break;

	case BSPConstructType::ObbMiddel:
		ConstructObbMiddleSons();
		break;

	case BSPConstructType::Gravity_SPLIT:
		ConstructGravity_SPLITSons();
		break;

	case BSPConstructType::SDM_OBB:
		ConstructGravity_SPLITSons();
		break;
	default:
		assert(false);
	}

}

//BSPTreeNode::BSPTreeNode(vector<TriangleBSP> p_FacesA, vector<TriangleBSP> p_FacesB, Plane3d p_PartitionPlane, Tolerance& p_Tolerance, bool p_Isroot /*= false*/, BSPConstructType p_Type /*= BSPConstructType::AABB*/):
//	m_FacesA(p_FacesA), m_FacesB(p_FacesB), m_Type(p_Type), m_PartitionPlane(p_PartitionPlane), m_Tolerance(p_Tolerance), m_IsRoot(p_Isroot)
//{
//
//
//}

void BSPTreeNode::GetLeafNodes(vector<BSPTreeNode*>& r_AllLeafNodes)
{
	if (m_IsLeaf) {
		r_AllLeafNodes.push_back(this);
	}
	else {

		assert(nullptr != m_PositiveSon);
		assert(nullptr != m_NegativeSon);

		m_PositiveSon->GetLeafNodes(r_AllLeafNodes);
		m_NegativeSon->GetLeafNodes(r_AllLeafNodes);
	}

}

int BSPTreeNode::Size()
{
	return m_FacesA.size() + m_FacesB.size();
}

bool BSPTreeNode::IsLeaf()
{
	return m_IsLeaf;
}

void BSPTreeNode::GetIntersectPair(IntersectTriangleCheckList& p_IntersectCheckList)
{
	for (auto& fa : m_FacesA) {
		for (auto& fb : m_FacesB) {

			p_IntersectCheckList.push_back(make_pair(fa.m_HoldFace.idx(), fb.m_HoldFace.idx()));

		}
	}

}

bool BSPTreeNode::IsValidNode()
{
	if ((0 == m_FacesA.size()) || (0 == m_FacesB.size()))
		return true;

	return false;
}

void BSPTreeNode::GetValidFaces(unordered_set<FaceId>& facesA, unordered_set<FaceId>& facesB)
{
	if ((0 == m_FacesA.size()) || (0 == m_FacesB.size()))
		return;

	for (auto& fa : m_FacesA) {
		facesA.insert(fa.m_HoldFace.idx());
	}

	for (auto& fb : m_FacesB) {
		facesB.insert(fb.m_HoldFace.idx());
	}

}

void BSPTreeNode::GetPartitionPlane(bspFaceRenderingInfo& r_Pplanes, planeRecords p_CurrentRecords)
{
	if (m_IsLeaf) return;


	if (m_IsRoot) {
		assert(m_Depth == 0);

		p_CurrentRecords.clear();
	}


	r_Pplanes.push_back(make_pair(m_PartitionPlane, p_CurrentRecords));



	p_CurrentRecords.push_back(make_pair(m_PartitionPlane, true));
	m_PositiveSon->GetPartitionPlane(r_Pplanes, p_CurrentRecords);

	p_CurrentRecords.pop_back();
	p_CurrentRecords.push_back(make_pair(m_PartitionPlane, false));
	m_NegativeSon->GetPartitionPlane(r_Pplanes, p_CurrentRecords);


}

void BSPTreeNode::ConstructAABBSons()
{
	vector<TriangleBSP> allFaces(m_FacesA);
	allFaces.insert(allFaces.end(), m_FacesB.begin(), m_FacesB.end());

	AABB aabb(allFaces);
	Plane3d partitionPlane = aabb.GetLongestAxisEqualParitionPlane();

	vector<TriangleBSP> rPositive, rNegative;
	GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


	vector<TriangleBSP> facesAPositive, facesANegative;
	vector<TriangleBSP> facesBPositive, facesBNegative;

	SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
	SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

	//建立positive子树
	m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum);

	m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum);

	m_IsLeaf = false;

	m_PartitionPlane = partitionPlane;
}

void BSPTreeNode::ConstructSDMSons()
{
	auto fFirst = m_FacesA;
	auto fSecond = m_FacesB;

	//ofstream fileA("PointInfoA.txt");
	//for (auto& fa : fFirst) {
	//	Triangle3d tri = fa.m_Triangle;

	//	for (int i = 0;i < 3;i++) {
	//		Point3d point = tri.VertexAt(i);

	//		fileA << "(" << point[0] << "," << point[1] << "," << point[2] << ")" << endl;
	//	}

	//}

	//ofstream fileB("PointInfoB.txt");
	//for (auto& fb : fSecond) {
	//	Triangle3d tri = fb.m_Triangle;

	//	for (int i = 0;i < 3;i++) {
	//		Point3d point = tri.VertexAt(i);

	//		fileB << "(" << point[0] << "," << point[1] << "," << point[2] << ")" << endl;
	//	}

	//}

	double sumXA = 0;
	double sumYA = 0;
	double sumZA = 0;

	double sumXB = 0;
	double sumYB = 0;
	double sumZB = 0;


	double standVarXA = 0;
	double standVarYA = 0;
	double standVarZA = 0;

	for (auto& fa : fFirst) {
		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			sumXA += point[0];
			sumYA += point[1];
			sumZA += point[2];
		}
	}

	double xBarA = sumXA / (fFirst.size() * 3);
	double yBarA = sumYA / (fFirst.size() * 3);
	double zBarA = sumZA / (fFirst.size() * 3);

	for (auto& fa : fFirst) {
		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			standVarXA += (point[0] - xBarA) * (point[0] - xBarA);
			standVarYA += (point[1] - yBarA) * (point[1] - yBarA);
			standVarZA += (point[2] - zBarA) * (point[2] - zBarA);
		}

	}

	standVarXA = sqrt(standVarXA / (fFirst.size() * 3 - 1));
	standVarYA = sqrt(standVarYA / (fFirst.size() * 3 - 1));
	standVarZA = sqrt(standVarZA / (fFirst.size() * 3 - 1));

	double scaleA = sqrt(standVarXA * standVarYA * standVarZA);



	//meshB
	double standVarXB = 0;
	double standVarYB = 0;
	double standVarZB = 0;

	for (auto& fb : fSecond) {
		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			sumXB += point[0];
			sumYB += point[1];
			sumZB += point[2];
		}
	}

	double xBarB = sumXB / (fSecond.size() * 3);
	double yBarB = sumYB / (fSecond.size() * 3);
	double zBarB = sumZB / (fSecond.size() * 3);

	for (auto& fb : fSecond) {

		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			standVarXB += (point[0] - xBarB) * (point[0] - xBarB);
			standVarYB += (point[1] - yBarB) * (point[1] - yBarB);
			standVarZB += (point[2] - zBarB) * (point[2] - zBarB);
		}


	}

	standVarXB = sqrt(standVarXB / (fSecond.size() * 3 - 1));
	standVarYB = sqrt(standVarYB / (fSecond.size() * 3 - 1));
	standVarZB = sqrt(standVarZB / (fSecond.size() * 3 - 1));
	double scaleB = sqrt(standVarXB * standVarYB * standVarZB);


	//weight
	double disWeightA = scaleB / (scaleA + scaleB);
	double disWeightB = scaleA / (scaleA + scaleB);

	double numWeightA = (double)fSecond.size() / (double)(fFirst.size() + fSecond.size());
	double numWeightB = (double)fFirst.size() / (double)(fFirst.size() + fSecond.size());

	ofstream outFile("..//TestModel//thingi10k//res2.txt",ios::app);
	if (m_Depth == 0) {
		if (scaleA / scaleB <= 1) {
			outFile << "distribution compare " << scaleB / scaleA << endl;
		}
		else {
			outFile << "distribution compare " << scaleA / scaleB << endl;
		}
		
		if (fFirst.size() >= fSecond.size()) {
			outFile << "size compare " << (double)fFirst.size() / (double)fSecond.size() << endl;
		}
		else {
			outFile << "size compare " << (double)fSecond.size() / (double)fFirst.size() << endl;
		}

	}
	outFile.close();

	//
	double lambdaS = 1, lambdaD = 0.5;

	double weightA = pow(disWeightA, lambdaD) * pow(numWeightA, lambdaS);
	double weightB = pow(disWeightB, lambdaD) * pow(numWeightB, lambdaS);

	double partialA = -1, partialB = -1, partialC = -1;
	if (disWeightA > disWeightB)
		partialA = disWeightA / disWeightB;
	else
		partialA = disWeightB / disWeightA;

	if (numWeightA > numWeightB)
		partialB = numWeightA / numWeightB;
	else
		partialB = numWeightB / numWeightA;

	if (weightA > weightB)
		partialC = weightA / weightB;
	else
		partialC = weightB / weightA;

	//double weightA = disWeightA * numWeightA;
	//double weightB = disWeightB * numWeightB;




	//double weightA = numWeightA;
	//double weightB = numWeightB;
	//double weightA = disWeightA;
	//double weightB = disWeightB;

	double centerX = sumXA * weightA + sumXB * weightB;
	double centerY = sumYA * weightA + sumYB * weightB;
	double centerZ = sumZA * weightA + sumZB * weightB;


	centerX /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);
	centerY /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);
	centerZ /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);



	Eigen::Matrix3d paraMatrix;
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			paraMatrix(i, j) = 0;
		}
	}


	//cout << "==========" << endl;
	//cout << paraMatrix<< endl;

	for (auto& fa : fFirst) {

		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			paraMatrix(0, 0) += weightA * pow(point[0] - centerX, 2);
			paraMatrix(1, 1) += weightA * pow(point[1] - centerY, 2);
			paraMatrix(2, 2) += weightA * pow(point[2] - centerZ, 2);

			paraMatrix(0, 1) += weightA * (point[0] - centerX) * (point[1] - centerY);
			paraMatrix(1, 0) += weightA * (point[0] - centerX) * (point[1] - centerY);

			paraMatrix(0, 2) += weightA * (point[0] - centerX) * (point[2] - centerZ);
			paraMatrix(2, 0) += weightA * (point[0] - centerX) * (point[2] - centerZ);

			paraMatrix(1, 2) += weightA * (point[1] - centerY) * (point[2] - centerZ);
			paraMatrix(2, 1) += weightA * (point[1] - centerY) * (point[2] - centerZ);
		}


	}

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;

	for (auto& fb : fSecond) {

		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			paraMatrix(0, 0) += weightB * pow(point[0] - centerX, 2);
			paraMatrix(1, 1) += weightB * pow(point[1] - centerY, 2);
			paraMatrix(2, 2) += weightB * pow(point[2] - centerZ, 2);

			paraMatrix(0, 1) += weightB * (point[0] - centerX) * (point[1] - centerY);
			paraMatrix(1, 0) += weightB * (point[0] - centerX) * (point[1] - centerY);

			paraMatrix(0, 2) += weightB * (point[0] - centerX) * (point[2] - centerZ);
			paraMatrix(2, 0) += weightB * (point[0] - centerX) * (point[2] - centerZ);

			paraMatrix(1, 2) += weightB * (point[1] - centerY) * (point[2] - centerZ);
			paraMatrix(2, 1) += weightB * (point[1] - centerY) * (point[2] - centerZ);
		}


	}

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;


	Eigen::EigenSolver<Eigen::MatrixXd> eSolver(paraMatrix);
	Eigen::MatrixXcd cEvecs = eSolver.eigenvectors();
	Eigen::MatrixXcd cEvals = eSolver.eigenvalues();

	Eigen::MatrixXd rEvecs = cEvecs.real();
	Eigen::MatrixXd rEvals = cEvals.real();

	Eigen::MatrixXf::Index maxEigvalusPos;
	rEvals.rowwise().sum().maxCoeff(&maxEigvalusPos);

	int maxEValus = int(maxEigvalusPos);

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;

	Vector3d planeNormal(rEvecs(0, maxEValus), rEvecs(1, maxEValus), rEvecs(2, maxEValus));


	//cout << "normal: " << rEvecs(0, maxEValus) << ", " << rEvecs(1, maxEValus) << ", " << rEvecs(2, maxEValus) << endl;
	//cout << "center: " << centerX << ", " << centerY << ", " << centerZ << endl;

	Plane3d partitionPlane(Point3d(centerX, centerY, centerZ), planeNormal);


	vector<TriangleBSP> rPositive, rNegative;
	GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


	vector<TriangleBSP> facesAPositive, facesANegative;
	vector<TriangleBSP> facesBPositive, facesBNegative;

	SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
	SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

	//建立positive子树
	m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SDM);

	m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SDM);

	m_IsLeaf = false;

	m_PartitionPlane = partitionPlane;

}

void BSPTreeNode::ConstructSAHSons()
{
	const size_t smallNodeThreshold = 256;


	if (m_FacesA.size() + m_FacesB.size() > smallNodeThreshold) {

		//large node middle split

		vector<TriangleBSP> allFaces(m_FacesA);
		allFaces.insert(allFaces.end(), m_FacesB.begin(), m_FacesB.end());

		AABB aabb(allFaces);
		Plane3d partitionPlane = aabb.GetLongestAxisEqualParitionPlane();

		vector<TriangleBSP> rPositive, rNegative;
		GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


		vector<TriangleBSP> facesAPositive, facesANegative;
		vector<TriangleBSP> facesBPositive, facesBNegative;

		SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
		SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

		//建立positive子树
		m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SAH);

		m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SAH);

		m_IsLeaf = false;

		m_PartitionPlane = partitionPlane;

	}
	else {

		
		//small node SAH
		//not split triangles
		vector<TriangleBSP> allFaces(m_FacesA);
		allFaces.insert(allFaces.end(), m_FacesB.begin(), m_FacesB.end());

		AABB aabb(allFaces);

		double currentSAHCost = double(m_Depth) + m_FacesA.size() + m_FacesB.size();

		//Plane3d partitionPlane = aabb.GetLongestAxisEqualParitionPlane();

		double minSAHCost = currentSAHCost;


		bool hasBetterSAH = false;
		Plane3d sahPlane;
		vector<TriangleBSP> rSahPositive, rSahNegative;


		for (auto& triCandidate : allFaces) {

			vector<TriangleBSP> rPositive, rNegative;

			bool sd = triCandidate.m_Triangle.isValid();
			double l1 = (triCandidate.m_Triangle.m_Vertices[0] - triCandidate.m_Triangle.m_Vertices[1]).length();
			double l2 = (triCandidate.m_Triangle.m_Vertices[1] - triCandidate.m_Triangle.m_Vertices[2]).length();
			double l3 = (triCandidate.m_Triangle.m_Vertices[2] - triCandidate.m_Triangle.m_Vertices[0]).length();

			Plane3d sahCandidatePlane = triCandidate.m_Triangle;

			GetNegativePositiveLabelWithOutSplitTriangle(rPositive, rNegative, sahCandidatePlane);

			if (rPositive.empty() || rNegative.empty()) continue;

			AABB positiveAABB(rPositive), negativeAABB(rNegative);


			double sAHCost = double(m_Depth) + rPositive.size() * positiveAABB.GetSurfaceArea() / aabb.GetSurfaceArea()
				+ rNegative.size() * negativeAABB.GetSurfaceArea() / aabb.GetSurfaceArea();

			if (sAHCost < minSAHCost) {
				minSAHCost = sAHCost;
				sahPlane = sahCandidatePlane;
				rSahPositive = rPositive;
				rSahNegative = rNegative;
				hasBetterSAH = true;
			}


		}

		if (hasBetterSAH) {

			vector<TriangleBSP> facesAPositive, facesANegative;
			vector<TriangleBSP> facesBPositive, facesBNegative;

			SeparationTrianglesFromDifferentMesh(rSahPositive, facesAPositive, facesBPositive);
			SeparationTrianglesFromDifferentMesh(rSahNegative, facesANegative, facesBNegative);

			//建立positive子树
			m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SAH);

			m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SAH);

			m_IsLeaf = false;

			m_PartitionPlane = sahPlane;

		}
		else {
			//all sah cost increase
			//not construct depth+1 son node

			m_IsLeaf = true;
			m_PositiveSon = nullptr;
			m_NegativeSon = nullptr;

			return;
		}

	}



}

void BSPTreeNode::ConstructObbMiddleSons()
{
	auto fFirst = m_FacesA;
	auto fSecond = m_FacesB;

	double sumXA = 0;
	double sumYA = 0;
	double sumZA = 0;

	for (auto& fa : fFirst) {
		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			sumXA += point[0];
			sumYA += point[1];
			sumZA += point[2];
		}
	}

	//meshB
	double sumXB = 0;
	double sumYB = 0;
	double sumZB = 0;

	for (auto& fb : fSecond) {
		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			sumXB += point[0];
			sumYB += point[1];
			sumZB += point[2];
		}
	}

	double xBarB = sumXB / (fSecond.size() * 3);
	double yBarB = sumYB / (fSecond.size() * 3);
	double zBarB = sumZB / (fSecond.size() * 3);

	double xBar = (sumXA + sumXB) / (fFirst.size() * 3 + fSecond.size() * 3);
	double yBar = (sumYA + sumYB) / (fFirst.size() * 3 + fSecond.size() * 3);
	double zBar = (sumZA + sumZB) / (fFirst.size() * 3 + fSecond.size() * 3);

	//初始化协方差矩阵
	Eigen::Matrix3d paraMatrix;
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			paraMatrix(i, j) = 0;
		}
	}

	//计算协方差矩阵中来自MeshA的分量
	for (auto& fa : fFirst) {

		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			paraMatrix(0, 0) += pow(point[0] - xBar, 2);
			paraMatrix(1, 1) += pow(point[1] - yBar, 2);
			paraMatrix(2, 2) += pow(point[2] - zBar, 2);

			paraMatrix(0, 1) += (point[0] - xBar) * (point[1] - yBar);
			paraMatrix(1, 0) += (point[0] - xBar) * (point[1] - yBar);

			paraMatrix(0, 2) += (point[0] - xBar) * (point[2] - zBar);
			paraMatrix(2, 0) += (point[0] - xBar) * (point[2] - zBar);

			paraMatrix(1, 2) += (point[1] - yBar) * (point[2] - zBar);
			paraMatrix(2, 1) += (point[1] - yBar) * (point[2] - zBar);
		}


	}

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;

	for (auto& fb : fSecond) {

		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			paraMatrix(0, 0) += pow(point[0] - xBar, 2);
			paraMatrix(1, 1) += pow(point[1] - yBar, 2);
			paraMatrix(2, 2) += pow(point[2] - zBar, 2);

			paraMatrix(0, 1) += (point[0] - xBar) * (point[1] - yBar);
			paraMatrix(1, 0) += (point[0] - xBar) * (point[1] - yBar);

			paraMatrix(0, 2) += (point[0] - xBar) * (point[2] - zBar);
			paraMatrix(2, 0) += (point[0] - xBar) * (point[2] - zBar);

			paraMatrix(1, 2) += (point[1] - yBar) * (point[2] - zBar);
			paraMatrix(2, 1) += (point[1] - yBar) * (point[2] - zBar);
		}


	}

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;


	Eigen::EigenSolver<Eigen::MatrixXd> eSolver(paraMatrix);
	Eigen::MatrixXcd cEvecs = eSolver.eigenvectors();
	Eigen::MatrixXcd cEvals = eSolver.eigenvalues();

	Eigen::MatrixXd rEvecs = cEvecs.real();
	Eigen::MatrixXd rEvals = cEvals.real();

	Eigen::MatrixXf::Index maxEigvalusPos;
	rEvals.rowwise().sum().maxCoeff(&maxEigvalusPos);

	int maxEValus = int(maxEigvalusPos);

	//cout << "==========" << endl;
	//cout << paraMatrix << endl;

	


	//cout << "normal: " << rEvecs(0, maxEValus) << ", " << rEvecs(1, maxEValus) << ", " << rEvecs(2, maxEValus) << endl;
	//cout << "center: " << centerX << ", " << centerY << ", " << centerZ << endl;

	Vector3d eigenVec1(rEvecs(0, 0), rEvecs(1, 0), rEvecs(2, 0));
	Vector3d eigenVec2(rEvecs(0, 1), rEvecs(1, 1), rEvecs(2, 1));
	Vector3d eigenVec3(rEvecs(0, 2), rEvecs(1, 2), rEvecs(2, 2));

	//计算投影在三个主方向上的中点
	

	Point3d initialPoint = fFirst.front().m_Triangle.VertexAt(0);

	double minProjectX = eigenVec1.dot(initialPoint);
	double maxProjectX = eigenVec1.dot(initialPoint);

	double minProjectY = eigenVec2.dot(initialPoint);
	double maxProjectY = eigenVec2.dot(initialPoint);

	double minProjectZ = eigenVec3.dot(initialPoint);
	double maxProjectZ = eigenVec3.dot(initialPoint);


	for (auto& fa : fFirst) {

		Triangle3d tri = fa.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			minProjectX = min(eigenVec1.dot(point), minProjectX);
			maxProjectX = max(eigenVec1.dot(point), maxProjectX);

			minProjectY = min(eigenVec2.dot(point), minProjectY);
			maxProjectY = max(eigenVec2.dot(point), maxProjectY);

			minProjectZ = min(eigenVec3.dot(point), minProjectZ);
			maxProjectZ = max(eigenVec3.dot(point), maxProjectZ);

		}
	}

	for (auto& fb : fSecond) {

		Triangle3d tri = fb.m_Triangle;

		for (int i = 0;i < 3;i++) {
			Point3d point = tri.VertexAt(i);

			minProjectX = min(eigenVec1.dot(point), minProjectX);
			maxProjectX = max(eigenVec1.dot(point), maxProjectX);

			minProjectY = min(eigenVec2.dot(point), minProjectY);
			maxProjectY = max(eigenVec2.dot(point), maxProjectY);

			minProjectZ = min(eigenVec3.dot(point), minProjectZ);
			maxProjectZ = max(eigenVec3.dot(point), maxProjectZ);
		}
	}

	Point3d obbCenter( (minProjectX + maxProjectX)/2 * eigenVec1 + (minProjectY + maxProjectY) / 2 * eigenVec2 + (minProjectZ + maxProjectZ) / 2 * eigenVec3);

	Vector3d planeNormal;
	//Vector3d planeNormal(rEvecs(0, maxEValus), rEvecs(1, maxEValus), rEvecs(2, maxEValus));

	if ((maxProjectX - minProjectX) > (maxProjectY - minProjectY)) {
		if ((maxProjectX - minProjectX) > (maxProjectZ - minProjectZ)) {
			//max x
			planeNormal = Vector3d(rEvecs(0, 0), rEvecs(1, 0), rEvecs(2, 0));
		}
		else {
			//max z
			planeNormal = Vector3d(rEvecs(0, 2), rEvecs(1, 2), rEvecs(2, 2));
		}
	
	}
	else {

		if ((maxProjectY - minProjectY) > (maxProjectZ - minProjectZ)) {
			//max y
			planeNormal = Vector3d(rEvecs(0, 1), rEvecs(1, 1), rEvecs(2, 1));
		}
		else {
			//max z
			planeNormal = Vector3d(rEvecs(0, 0), rEvecs(1, 0), rEvecs(2, 0));
		}
	
	
	}

	Plane3d partitionPlane(obbCenter, planeNormal);


	vector<TriangleBSP> rPositive, rNegative;
	GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


	vector<TriangleBSP> facesAPositive, facesANegative;
	vector<TriangleBSP> facesBPositive, facesBNegative;

	SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
	SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

	//if (double(rPositive.size() + rNegative.size()) >= 10 * double(m_FacesA.size() + m_FacesB.size())) {
	//	m_PositiveSon = nullptr;
	//	m_NegativeSon = nullptr;
	//	m_IsLeaf = true;

	//	return;
	//}

	//建立positive子树
	m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::ObbMiddel);

	m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::ObbMiddel);

	m_IsLeaf = false;

	m_PartitionPlane = partitionPlane;
}

void BSPTreeNode::ConstructGravity_SPLITSons()
{
	vector<TriangleBSP> allFaces(m_FacesA);
	allFaces.insert(allFaces.end(), m_FacesB.begin(), m_FacesB.end());

	AABB aabb(allFaces);
	Vector3d directionLongest = aabb.GetLongestAxis();

	double areaSum = 0;
	Point3d gravityPoint(0, 0, 0);

	for (auto& face : allFaces) {
		Vector3d v1 = face.m_Triangle.VertexAt(1) - face.m_Triangle.VertexAt(0);
		Vector3d v2 = face.m_Triangle.VertexAt(2) - face.m_Triangle.VertexAt(0);

		double triArea = v1.cross(v2).length() / 2;
		areaSum += triArea;

		gravityPoint += (face.m_Triangle.VertexAt(0) + face.m_Triangle.VertexAt(1) + face.m_Triangle.VertexAt(2)) * triArea;


	}

	gravityPoint /= 3;
	gravityPoint /= areaSum;

	Plane3d partitionPlane(gravityPoint, directionLongest);

	vector<TriangleBSP> rPositive, rNegative;
	GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


	vector<TriangleBSP> facesAPositive, facesANegative;
	vector<TriangleBSP> facesBPositive, facesBNegative;

	SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
	SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

	//建立positive子树
	m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum);

	m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum);

	m_IsLeaf = false;

	m_PartitionPlane = partitionPlane;
}

void BSPTreeNode::ConstructSDM_OBBSons()
{
	if (m_Depth <= 8) {
		auto fFirst = m_FacesA;
		auto fSecond = m_FacesB;

		//ofstream fileA("PointInfoA.txt");
		//for (auto& fa : fFirst) {
		//	Triangle3d tri = fa.m_Triangle;

		//	for (int i = 0;i < 3;i++) {
		//		Point3d point = tri.VertexAt(i);

		//		fileA << "(" << point[0] << "," << point[1] << "," << point[2] << ")" << endl;
		//	}

		//}

		//ofstream fileB("PointInfoB.txt");
		//for (auto& fb : fSecond) {
		//	Triangle3d tri = fb.m_Triangle;

		//	for (int i = 0;i < 3;i++) {
		//		Point3d point = tri.VertexAt(i);

		//		fileB << "(" << point[0] << "," << point[1] << "," << point[2] << ")" << endl;
		//	}

		//}

		double sumXA = 0;
		double sumYA = 0;
		double sumZA = 0;

		double sumXB = 0;
		double sumYB = 0;
		double sumZB = 0;


		double standVarXA = 0;
		double standVarYA = 0;
		double standVarZA = 0;

		for (auto& fa : fFirst) {
			Triangle3d tri = fa.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				sumXA += point[0];
				sumYA += point[1];
				sumZA += point[2];
			}
		}

		double xBarA = sumXA / (fFirst.size() * 3);
		double yBarA = sumYA / (fFirst.size() * 3);
		double zBarA = sumZA / (fFirst.size() * 3);

		for (auto& fa : fFirst) {
			Triangle3d tri = fa.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				standVarXA += (point[0] - xBarA) * (point[0] - xBarA);
				standVarYA += (point[1] - yBarA) * (point[1] - yBarA);
				standVarZA += (point[2] - zBarA) * (point[2] - zBarA);
			}

		}

		standVarXA = sqrt(standVarXA / (fFirst.size() * 3 - 1));
		standVarYA = sqrt(standVarYA / (fFirst.size() * 3 - 1));
		standVarZA = sqrt(standVarZA / (fFirst.size() * 3 - 1));

		double scaleA = sqrt(standVarXA * standVarYA * standVarZA);



		//meshB
		double standVarXB = 0;
		double standVarYB = 0;
		double standVarZB = 0;

		for (auto& fb : fSecond) {
			Triangle3d tri = fb.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				sumXB += point[0];
				sumYB += point[1];
				sumZB += point[2];
			}
		}

		double xBarB = sumXB / (fSecond.size() * 3);
		double yBarB = sumYB / (fSecond.size() * 3);
		double zBarB = sumZB / (fSecond.size() * 3);

		for (auto& fb : fSecond) {

			Triangle3d tri = fb.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				standVarXB += (point[0] - xBarB) * (point[0] - xBarB);
				standVarYB += (point[1] - yBarB) * (point[1] - yBarB);
				standVarZB += (point[2] - zBarB) * (point[2] - zBarB);
			}


		}

		standVarXB = sqrt(standVarXB / (fSecond.size() * 3 - 1));
		standVarYB = sqrt(standVarYB / (fSecond.size() * 3 - 1));
		standVarZB = sqrt(standVarZB / (fSecond.size() * 3 - 1));
		double scaleB = sqrt(standVarXB * standVarYB * standVarZB);


		//weight
		double disWeightA = scaleB / (scaleA + scaleB);
		double disWeightB = scaleA / (scaleA + scaleB);

		double numWeightA = (double)fSecond.size() / (double)(fFirst.size() + fSecond.size());
		double numWeightB = (double)fFirst.size() / (double)(fFirst.size() + fSecond.size());



		//
		double lambdaS = 1, lambdaD = 0.5;

		double weightA = pow(disWeightA, lambdaD) * pow(numWeightA, lambdaS);
		double weightB = pow(disWeightB, lambdaD) * pow(numWeightB, lambdaS);

		double partialA = -1, partialB = -1, partialC = -1;
		if (disWeightA > disWeightB)
			partialA = disWeightA / disWeightB;
		else
			partialA = disWeightB / disWeightA;

		if (numWeightA > numWeightB)
			partialB = numWeightA / numWeightB;
		else
			partialB = numWeightB / numWeightA;

		if (weightA > weightB)
			partialC = weightA / weightB;
		else
			partialC = weightB / weightA;

		//double weightA = disWeightA * numWeightA;
		//double weightB = disWeightB * numWeightB;




		//double weightA = numWeightA;
		//double weightB = numWeightB;
		//double weightA = disWeightA;
		//double weightB = disWeightB;

		double centerX = sumXA * weightA + sumXB * weightB;
		double centerY = sumYA * weightA + sumYB * weightB;
		double centerZ = sumZA * weightA + sumZB * weightB;


		centerX /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);
		centerY /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);
		centerZ /= (weightA * fFirst.size() * 3 + weightB * fSecond.size() * 3);



		Eigen::Matrix3d paraMatrix;
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				paraMatrix(i, j) = 0;
			}
		}


		//cout << "==========" << endl;
		//cout << paraMatrix<< endl;

		for (auto& fa : fFirst) {

			Triangle3d tri = fa.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				paraMatrix(0, 0) += weightA * pow(point[0] - centerX, 2);
				paraMatrix(1, 1) += weightA * pow(point[1] - centerY, 2);
				paraMatrix(2, 2) += weightA * pow(point[2] - centerZ, 2);

				paraMatrix(0, 1) += weightA * (point[0] - centerX) * (point[1] - centerY);
				paraMatrix(1, 0) += weightA * (point[0] - centerX) * (point[1] - centerY);

				paraMatrix(0, 2) += weightA * (point[0] - centerX) * (point[2] - centerZ);
				paraMatrix(2, 0) += weightA * (point[0] - centerX) * (point[2] - centerZ);

				paraMatrix(1, 2) += weightA * (point[1] - centerY) * (point[2] - centerZ);
				paraMatrix(2, 1) += weightA * (point[1] - centerY) * (point[2] - centerZ);
			}


		}

		//cout << "==========" << endl;
		//cout << paraMatrix << endl;

		for (auto& fb : fSecond) {

			Triangle3d tri = fb.m_Triangle;

			for (int i = 0;i < 3;i++) {
				Point3d point = tri.VertexAt(i);

				paraMatrix(0, 0) += weightB * pow(point[0] - centerX, 2);
				paraMatrix(1, 1) += weightB * pow(point[1] - centerY, 2);
				paraMatrix(2, 2) += weightB * pow(point[2] - centerZ, 2);

				paraMatrix(0, 1) += weightB * (point[0] - centerX) * (point[1] - centerY);
				paraMatrix(1, 0) += weightB * (point[0] - centerX) * (point[1] - centerY);

				paraMatrix(0, 2) += weightB * (point[0] - centerX) * (point[2] - centerZ);
				paraMatrix(2, 0) += weightB * (point[0] - centerX) * (point[2] - centerZ);

				paraMatrix(1, 2) += weightB * (point[1] - centerY) * (point[2] - centerZ);
				paraMatrix(2, 1) += weightB * (point[1] - centerY) * (point[2] - centerZ);
			}


		}

		//cout << "==========" << endl;
		//cout << paraMatrix << endl;


		Eigen::EigenSolver<Eigen::MatrixXd> eSolver(paraMatrix);
		Eigen::MatrixXcd cEvecs = eSolver.eigenvectors();
		Eigen::MatrixXcd cEvals = eSolver.eigenvalues();

		Eigen::MatrixXd rEvecs = cEvecs.real();
		Eigen::MatrixXd rEvals = cEvals.real();

		Eigen::MatrixXf::Index maxEigvalusPos;
		rEvals.rowwise().sum().maxCoeff(&maxEigvalusPos);

		int maxEValus = int(maxEigvalusPos);

		//cout << "==========" << endl;
		//cout << paraMatrix << endl;

		Vector3d planeNormal(rEvecs(0, maxEValus), rEvecs(1, maxEValus), rEvecs(2, maxEValus));


		//cout << "normal: " << rEvecs(0, maxEValus) << ", " << rEvecs(1, maxEValus) << ", " << rEvecs(2, maxEValus) << endl;
		//cout << "center: " << centerX << ", " << centerY << ", " << centerZ << endl;

		Plane3d partitionPlane(Point3d(centerX, centerY, centerZ), planeNormal);


		vector<TriangleBSP> rPositive, rNegative;
		GetNegativePositiveLabelWithSplitTriangle(rPositive, rNegative, partitionPlane);


		vector<TriangleBSP> facesAPositive, facesANegative;
		vector<TriangleBSP> facesBPositive, facesBNegative;

		SeparationTrianglesFromDifferentMesh(rPositive, facesAPositive, facesBPositive);
		SeparationTrianglesFromDifferentMesh(rNegative, facesANegative, facesBNegative);

		//建立positive子树
		m_PositiveSon = new BSPTreeNode(facesAPositive, facesBPositive, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SDM_OBB);

		m_NegativeSon = new BSPTreeNode(facesANegative, facesBNegative, m_Tolerance, m_Depth + 1, m_MaxDepth, m_LeafShresholdNum, false, BSPConstructType::SDM_OBB);

		m_IsLeaf = false;

		m_PartitionPlane = partitionPlane;
	}
	else {
	
	
	}

}

void BSPTreeNode::GetNegativePositiveLabelWithSplitTriangle(vector<TriangleBSP>& negativeFaces, vector<TriangleBSP>& positiveFaces, Plane3d partitionPlane)
{
	negativeFaces.clear();
	positiveFaces.clear();

	for (auto& face : m_FacesA) {
		TrianglePlaneIsIntersectInfo info;

		Triangle3d tri = face.m_Triangle;

		bool canJudgeFaceA = true;
		if (IsZero(tri.GetNormal().cross(partitionPlane.GetNormal()))) {
			if (IsPositive((tri.GetCenter() - partitionPlane.GetOrigin()).dot(partitionPlane.GetNormal()))) {
				info.isIntersect = false;
				info.direction = true;
			}
			else {
				info.isIntersect = false;
				info.direction = false;
			}

		}
		else {
			canJudgeFaceA = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);
		}


		if (!canJudgeFaceA) {
			//canJudgeFaceA = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);
			canJudgeFaceA = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);
			positiveFaces.push_back(face);
			negativeFaces.push_back(face);
		}
		else {
			if (info.isIntersect) {
				for (auto posiTri : info.newPositiveTriangles) {

					positiveFaces.push_back({ posiTri, face.m_HoldFace,face.m_holdMesh });
				}

				for (auto negaTri : info.newNegativeTriangles) {

					negativeFaces.push_back({ negaTri, face.m_HoldFace,face.m_holdMesh });
				}

				//negativeFaces.push_back(face);
				//positiveFaces.push_back(face);
			}
			else {
				if (false == info.direction) {
					negativeFaces.push_back(face);
				}
				else {
					positiveFaces.push_back(face);
				}
			}
		
		}


	}

	//AABB aabb1(positiveFaces);
	//AABB aabb2(negativeFaces);

	for (auto& face : m_FacesB) {

		TrianglePlaneIsIntersectInfo info;

		Triangle3d tri = face.m_Triangle;

		bool canJudgeFaceB = true;
		canJudgeFaceB = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);

		if (!canJudgeFaceB) {
			positiveFaces.push_back(face);
			negativeFaces.push_back(face);
		}	
		else {
			if (info.isIntersect) {
				for (auto posiTri : info.newPositiveTriangles) {

					positiveFaces.push_back({ posiTri, face.m_HoldFace,face.m_holdMesh });
				}

				for (auto negaTri : info.newNegativeTriangles) {

					negativeFaces.push_back({ negaTri, face.m_HoldFace,face.m_holdMesh });
				}

			}
			else
			{

				if (false == info.direction) {
					negativeFaces.push_back(face);
				}
				else {
					positiveFaces.push_back(face);
				}
			}
		}


	}
}

void BSPTreeNode::GetNegativePositiveLabelWithOutSplitTriangle(vector<TriangleBSP>& negativeFaces, vector<TriangleBSP>& positiveFaces, Plane3d partitionPlane)
{
	negativeFaces.clear();
	positiveFaces.clear();

	for (auto& face : m_FacesA) {
		TrianglePlaneIsIntersectInfo info;

		Triangle3d tri = face.m_Triangle;

		bool canJudgeFaceA = true;
		if (IsZero(tri.GetNormal().cross(partitionPlane.GetNormal()))) {
			if (IsPositive((tri.GetCenter() - partitionPlane.GetOrigin()).dot(partitionPlane.GetNormal()))) {
				info.isIntersect = false;
				info.direction = true;
			}
			else {
				info.isIntersect = false;
				info.direction = false;
			}

		}
		else {
			canJudgeFaceA = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);
		}


		if (!canJudgeFaceA) {
			//canJudgeFaceA = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);

			positiveFaces.push_back(face);
			negativeFaces.push_back(face);
		}
		else {
			if (info.isIntersect) {

				positiveFaces.push_back(face);
				negativeFaces.push_back(face);

				//negativeFaces.push_back(face);
				//positiveFaces.push_back(face);
			}
			else {
				if (false == info.direction) {
					negativeFaces.push_back(face);
				}
				else {
					positiveFaces.push_back(face);
				}
			}

		}


	}

	//AABB aabb1(positiveFaces);
	//AABB aabb2(negativeFaces);

	for (auto& face : m_FacesB) {

		TrianglePlaneIsIntersectInfo info;

		Triangle3d tri = face.m_Triangle;

		bool canJudgeFaceB = true;
		canJudgeFaceB = TrianglePlaneIsIntersect(tri, partitionPlane, info, m_Tolerance);

		if (!canJudgeFaceB) {

			positiveFaces.push_back(face);
			negativeFaces.push_back(face);
		}
		else {
			if (info.isIntersect) {

				positiveFaces.push_back(face);
				negativeFaces.push_back(face);

			}
			else
			{

				if (false == info.direction) {
					negativeFaces.push_back(face);
				}
				else {
					positiveFaces.push_back(face);
				}
			}
		}


	}
}


void BSPTreeNode::SeparationTrianglesFromDifferentMesh(vector<TriangleBSP> facesOrigin, vector<TriangleBSP>& facesA, vector<TriangleBSP>& facesB)
{
	assert(!m_FacesA.empty());
	assert(!m_FacesB.empty());

	Mesh* meshA = m_FacesA.front().m_holdMesh;
	Mesh* meshB = m_FacesB.front().m_holdMesh;

	for (auto& face : facesOrigin) {

		if (face.m_holdMesh == meshA) {
			facesA.push_back(face);
		}
		else if (face.m_holdMesh == meshB) {
			facesB.push_back(face);
		}
		else {
			assert(false);
		}
	}


}

pair<int, int> BSPTreeNode::GetRealSize()
{
	unordered_set<FaceId> fa, fb;

	for (auto& f : m_FacesA) {
		fa.insert(f.m_HoldFace.idx());
	}

	for (auto& f : m_FacesB) {
		fb.insert(f.m_HoldFace.idx());
	}

	return pair<int, int>(fa.size(),fb.size());
}
