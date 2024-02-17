#include "mPch.h"
#include "TestBool.h"

namespace TestFunction {

	void TestBool()
	{
		Mesh meshA, meshB;

		//GenerateMesh::CubicMesh(1, Mesh::Point(0, 0, 0), meshA);
		//GenerateMesh::CubicMesh(1, Mesh::Point(0.5, 0.5, 0.5), meshB);

		//if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bunny1.obj")){

		//	std::cerr << "Cannot read mesh" << std::endl;

		//	return;

		//}
		//if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bunny2.obj")) {

		//	std::cerr << "Cannot read mesh" << std::endl;

		//	return;

		//}

		//if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bool operate//1//cow.obj")) {

		//	std::cerr << "Cannot read mesh" << std::endl;

		//	return;

		//}
		//if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bool operate//1//cowTrans.obj")) {

		//	std::cerr << "Cannot read mesh" << std::endl;

		//	return;

		//}

		string m1 = "..//TestModel//bool operate//8//bunnyA.obj";
		string m2 = "..//TestModel//bool operate//8//bunnyB.obj";

		if (!OpenMesh::IO::read_mesh(meshB, m1)) {

			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}
		if (!OpenMesh::IO::read_mesh(meshA, m2)) {

			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}


		//ofstream of1("..//TestModel//bool operate//meanTime.txt");

		//for (int i = 0; i < 4;i++) {
		//	BSPConstructType bt;
		//	string s;
		//	switch (i)
		//	{
		//	case 0:
		//		bt = BSPConstructType::SDM;
		//		s = "SDM";
		//		break;
		//	case 1:
		//		bt = BSPConstructType::SAH;
		//		s = "SAH";
		//		break;
		//	case 2:
		//		bt = BSPConstructType::ObbMiddel;
		//		s = "OBB";
		//		break;
		//	case 3:
		//		bt = BSPConstructType::Gravity_SPLIT;
		//		s = "Gravity";
		//		break;
		//	default:
		//		break;
		//	}



		//	vector<double> ts;
		//	for (int i = 0; i < 5;i++) {


		//		if (!OpenMesh::IO::read_mesh(meshA, m1)) {

		//			std::cerr << "Cannot read mesh" << std::endl;

		//			return;

		//		}
		//		if (!OpenMesh::IO::read_mesh(meshB, m2)) {

		//			std::cerr << "Cannot read mesh" << std::endl;

		//			return;

		//		}


		//		auto t0 = std::chrono::steady_clock::now();

		//		//进行布尔运算
		//		Mesh result;
		//		MeshBoolOperate mBo(meshA, meshB, bt);

		//		mBo.Run(result, MeshBoolOperateType::Difference);
		//		//Mesh r_MeshA, r_MeshB;
		//		//mBo.RunTest(result, MeshBoolOperateType::Difference, r_MeshA, r_MeshB);

		//		auto t1 = std::chrono::steady_clock::now();
		//		//cout << "boolean time：" << std::chrono::duration<double, std::milli>(t1 - t0).count() << "毫秒" << endl;
		//		ts.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
		//		
		//	}

		//	double sum = 0;
		//	for (auto t : ts) {
		//		sum += t;
		//	}

		//	of1 << s << "mean time: " << sum / ts.size() << endl;


		//}

		//of1.close();




		//ofstream of("..//TestModel//bool operate//time.txt");

		//for (int i = 0; i < 4;i++) {
		//	BSPConstructType bt;
		//	string s;

		//	switch (i)
		//	{
		//	case 0:
		//		bt = BSPConstructType::SDM;
		//		s = "SDM";
		//		break;
		//	case 1:
		//		bt = BSPConstructType::SAH;
		//		s = "SAH";
		//		break;
		//	case 2:
		//		bt = BSPConstructType::ObbMiddel;
		//		s = "ObbMiddel";
		//		break;
		//	case 3:
		//		bt = BSPConstructType::Gravity_SPLIT;
		//		s = "Gravity_SPLIT";
		//		break;
		//	default:
		//		break;
		//	}
		//

		//	if (!OpenMesh::IO::read_mesh(meshA, m1)) {

		//		std::cerr << "Cannot read mesh" << std::endl;

		//		return;

		//	}
		//	if (!OpenMesh::IO::read_mesh(meshB, m2)) {

		//		std::cerr << "Cannot read mesh" << std::endl;

		//		return;

		//	}


		//	auto t0 = std::chrono::steady_clock::now();

		//	//进行布尔运算
		//	Mesh result;
		//	MeshBoolOperate mBo(meshA, meshB, bt);

		//	mBo.Run(result, MeshBoolOperateType::Difference);
		//	//Mesh r_MeshA, r_MeshB;
		//	//mBo.RunTest(result, MeshBoolOperateType::Difference, r_MeshA, r_MeshB);

		//	auto t1 = std::chrono::steady_clock::now();
		//	of << "boolean time：" << s <<" == " << std::chrono::duration<double, std::milli>(t1 - t0).count() << "毫秒" << endl;

		//}

		//of.close();

		auto t0 = std::chrono::steady_clock::now();

		//进行布尔运算
		Mesh result;
		MeshBoolOperate mBo(meshA, meshB, BSPConstructType::SDM);

		//mBo.Run(result, MeshBoolOperateType::Difference);
		Mesh r_MeshA, r_MeshB;
		mBo.RunTest(result, MeshBoolOperateType::Difference, r_MeshA, r_MeshB);

		auto t1 = std::chrono::steady_clock::now();
		cout << "boolean time：" << std::chrono::duration<double, std::milli>(t1 - t0).count() << "毫秒" << endl;




		if (!OpenMesh::IO::write_mesh(result, "..//TestModel//bool operate//Test.obj")) {

			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

			return;
		}

		if (!OpenMesh::IO::write_mesh(r_MeshA, "..//TestModel//bool operate//rMeshATest.obj")) {

			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

			return;

		}

		if (!OpenMesh::IO::write_mesh(r_MeshB, "..//TestModel//bool operate//rMeshBTest.obj")) {

			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

			return;

		}


		//结果显示

		bool showIntersect = false;
		bool showClassify = false;
		bool showResult = true;
		bool showA = false;
		bool showB = false;

		RenderWnd wnd;
		RenderScence scence;

		if (showIntersect) {
			auto& paths = mBo.m_IntersectResult.at(&meshB).GetIntersectPaths();
			auto& vers = mBo.m_IntersectResult.at(&meshB).GetIntersectVertices();
 
			//显示交点
			Vertex_Set intersectVertex;

			for (MeshIntersectVertex* p_ver : vers)
			{
				intersectVertex.push_back(*p_ver->m_Point3d);
			}
			{
				if (intersectVertex.size() > 0)
				{
					auto* ra = new VerticesRenderer<Vertex_Set>(intersectVertex);
					ra->SetVertexSize(15.0);
					ra->SetVertexColor(glm::vec3(1.0f, 0.0f, 0.0f));
					scence.AddUnit(ra);
				}
			}



			//{
			//	Vertex_Set intersectVertex;

			//	intersectVertex.push_back(Point3d(-0.1, 0.5, 0.2));

			//	{
			//		if (intersectVertex.size() > 0)
			//		{
			//			auto* ra = new VerticesRenderer<Vertex_Set>(intersectVertex);
			//			ra->SetVertexSize(15.0);
			//			ra->SetVertexColor(glm::vec3(0.0f, 0.0f, 1.0f));
			//			scence.AddUnit(ra);
			//		}
			//	}
			//
			//
			//}

			//显示相交路径
			for (MeshIntersectPath* p_Path : paths)
			{
				Edge_Set intersectPaths;

				size_t pathSize = p_Path->Size();

				for (size_t j = 0; j < pathSize; j++)
				{
					if (j > 0) {
						intersectPaths.push_back(*(*p_Path)[j - 1].m_Point3d);
						intersectPaths.push_back(*(*p_Path)[j].m_Point3d);
					}

				}

				{
					if (intersectPaths.size() > 0)
					{
						auto* ra = new EdgesRenderer<Edge_Set>(intersectPaths);
						ra->SetEdgeWidth(8.0);
						ra->SetEdgeColor(glm::vec3(1.0f, 0.0f, 0.0f));
						scence.AddUnit(ra);
					}
				}

			}
		
		
		}



		//显示分类结果
		if (showClassify) {
			Mesh& meshShow = meshA;

			ClassifyRecord& labelResult = mBo.m_TriangleRecord;
			unordered_map<FaceId, bool>& labelRecord = labelResult[&meshShow];

			pair<Triangle_Set, Triangle_Set> labeledTris;
			DataTranslate::TranslateMeshWithlable(meshShow, labelRecord, labeledTris);

			result.garbage_collection();
			result.request_face_normals();

			Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshShow);
			Edge_Set thisEdge = DataTranslate::TranslateEdge(meshShow);
			Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshShow);


			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(labeledTris.first);
				ra->SetTriangleColor(glm::vec3(0.0f, 1.0f, 1.0f));
				scence.AddUnit(ra);
			}

			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(labeledTris.second);
				ra->SetTriangleColor(glm::vec3(0.0f, 0.0f, 1.0f));
				scence.AddUnit(ra);
			}
		
		}



		//显示result
		if (showResult)
		{
			result.garbage_collection();

			Triangle_Set thisMesh = DataTranslate::TranslateMesh(result);
			Edge_Set thisEdge = DataTranslate::TranslateEdge(result);
			Vertex_Set thisVertex = DataTranslate::TranslateVertex(result);


			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
				ra->SetTriangleColor(glm::vec3(0.0f, 1.0f, 1.0f));
				scence.AddUnit(ra);
			}
 
 
			//{
			//	auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
			//	ra->SetVertexSize(15.0);
			//	ra->SetVertexColor(glm::vec3(0.0f, 1.0f, 0.0f));
			//	scence.AddUnit(ra);
			//}

			{
				auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

				//ra->SetVertexSize(15.0);
				ra->SetEdgeColor(glm::vec3(1.0f, 1.0f, 0.0f));
				scence.AddUnit(ra);
			}
		}


		//显示meshA
		if (showA)
		{
			Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshA);
			Edge_Set thisEdge = DataTranslate::TranslateEdge(meshA);
			Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshA);

			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);	
				ra->SetTriangleColor(glm::vec3(1.0f, 0.0f, 0.0f));
				scence.AddUnit(ra);
			}

			//{
			//	auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
			//	ra->SetVertexSize(15.0);
			//	ra->SetVertexColor(glm::vec3(0.0f, 1.0f, 1.0f));
			//	scence.AddUnit(ra);
			//}

			//{
			//	auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

			//	//ra->SetVertexSize(15.0);
			//	ra->SetEdgeColor(glm::vec3(0.0f, 1.0f, 0.0f));
			//	scence.AddUnit(ra);
			//}
		}

		//显示meshB
		if (showB)
		{
			Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshB);
			Edge_Set thisEdge = DataTranslate::TranslateEdge(meshB);
			Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshB);

			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
				ra->SetTriangleColor(glm::vec3(0.0f, 1.0f, 1.0f));
				scence.AddUnit(ra);
			}

			//{
			//	auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
			//	ra->SetVertexSize(15.0);
			//	ra->SetVertexColor(glm::vec3(1.0f, 1.0f, 0.0f));
			//	scence.AddUnit(ra);
			//}

			//{
			//	auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);
			//	ra->SetEdgeColor(glm::vec3(0.0f, 0.0f, 0.0f));
			//	scence.AddUnit(ra);
			//}
		}


		wnd.SetScence(&scence);
		wnd.Run();

		//try

		//{

		//	if (!OpenMesh::IO::write_mesh(result, "..//TestModel//out.off"))

		//	{

		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

		//		return;

		//	}

		//}

		//catch (std::exception& x)

		//{

		//	std::cerr << x.what() << std::endl;

		//	return;

		//}
	}

	void TestWrite()
	{
		Mesh meshA, meshB;

		//if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bunny.obj")) {

		//	std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

		//	return;

		//}

		//MeshNormalize(meshA,0.2);

		GenerateMesh::CubicMesh(0.9, Mesh::Point(0, 0, 0), meshB);
		MeshTransform(meshB, Vector3d(5, 0, 0));
		MeshTransform(meshB, Vector3d(-0.001, 0.1, -0.1));
		 
		//if (!OpenMesh::IO::write_mesh(meshA, "..//TestModel//bunnyA.obj")) {

		//	std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

		//	return;

		//}

		if (!OpenMesh::IO::write_mesh(meshB, "..//TestModel//cubicB.obj")) {

			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

			return;

		}

		//MeshTransform(meshA, Vector3d(0.5, 0, 0));
		//if (!OpenMesh::IO::write_mesh(meshA, "..//TestModel//bunny2.obj")) {

		//	std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

		//	return;
		//}
	}

	void TestReadObjWriteStl()
	{
		Mesh meshA,meshB;

		if (!OpenMesh::IO::read_mesh(meshA, "C:\\Users\\NUOSEN\\Desktop\\work\\BSP\\bool-new3-d\\TestModel\\bunny\\16\\bunnyA.STL")) {

			std::cerr << "Cannot write mesh to file 'bunny.obj'" << std::endl;

			return;

		}
		if (!OpenMesh::IO::read_mesh(meshB, "C:\\Users\\NUOSEN\\Desktop\\work\\BSP\\bool-new3-d\\TestModel\\bunny\\16\\bunnyB.STL")) {

			std::cerr << "Cannot write mesh to file 'bunny.obj'" << std::endl;

			return;

		}



		if (!OpenMesh::IO::write_mesh(meshA, "C:\\Users\\NUOSEN\\Desktop\\work\\BSP\\bool-new3-d\\TestModel\\bunny\\16\\bunnyA.obj")) {

			std::cerr << "Cannot write mesh to file 'bunny.stl'" << std::endl;

			return;

		}
		if (!OpenMesh::IO::write_mesh(meshB, "C:\\Users\\NUOSEN\\Desktop\\work\\BSP\\bool-new3-d\\TestModel\\bunny\\16\\bunnyB.obj")) {

			std::cerr << "Cannot write mesh to file 'bunny.stl'" << std::endl;

			return;

		}
	}
	
	void TestOpenMesh()
	{
		Mesh meshA, meshB;

		if (!OpenMesh::IO::read_mesh(meshA, "C:\\Users\\NUOSEN\\Desktop\\BSP\\bool-new3-d\\TestModel\\bunny\\6\\subtree1.obj"))
		{

			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

			return;

		}

		// write mesh to output.obj


		RenderWnd wnd;
		RenderScence scence;

		//显示meshA
		if (1)
		{
			Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshA);
			Edge_Set thisEdge = DataTranslate::TranslateEdge(meshA);
			Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshA);

			{
				auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
				ra->SetTriangleColor(glm::vec3(1.0f, 0.0f, 0.0f));
				scence.AddUnit(ra);
			}

			{
				auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
				ra->SetVertexSize(15.0);
				ra->SetVertexColor(glm::vec3(0.0f, 1.0f, 1.0f));
				scence.AddUnit(ra);
			}

			{
				auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

				//ra->SetVertexSize(15.0);
				ra->SetEdgeColor(glm::vec3(0.0f, 1.0f, 0.0f));
				scence.AddUnit(ra);
			}
		}

		////显示meshB
		//if (1)
		//{
		//	Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshB);
		//	Edge_Set thisEdge = DataTranslate::TranslateEdge(meshB);
		//	Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshB);

		//	{
		//		auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
		//		ra->SetTriangleColor(glm::vec3(1.0f, 0.0f, 0.0f));
		//		scence.AddUnit(ra);
		//	}

		//	{
		//		auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
		//		ra->SetVertexSize(15.0);
		//		ra->SetVertexColor(glm::vec3(0.0f, 1.0f, 1.0f));
		//		scence.AddUnit(ra);
		//	}

		//	{
		//		auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

		//		//ra->SetVertexSize(15.0);
		//		ra->SetEdgeColor(glm::vec3(0.0f, 1.0f, 0.0f));
		//		scence.AddUnit(ra);
		//	}
		//}


		wnd.SetScence(&scence);
		wnd.Run();


		//Mesh tMesh;

		//Mesh::VertexHandle vhandle[8];

		//vhandle[0] = tMesh.add_vertex(Mesh::Point(0, 0, 2));

		//vhandle[1] = tMesh.add_vertex(Mesh::Point(2, 0, 2));

		//vhandle[2] = tMesh.add_vertex(Mesh::Point(2, 2, 2));

		//vhandle[3] = tMesh.add_vertex(Mesh::Point(0, 2, 2));

		//vhandle[4] = tMesh.add_vertex(Mesh::Point(0, 0, 0));

		//vhandle[5] = tMesh.add_vertex(Mesh::Point(2, 0, 0));

		//vhandle[6] = tMesh.add_vertex(Mesh::Point(2, 2, 0));

		//vhandle[7] = tMesh.add_vertex(Mesh::Point(0, 2, 0));

		//// generate (quadrilateral) faces

		//std::vector<Mesh::VertexHandle>  face_vhandles;

		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[0]);

		//face_vhandles.push_back(vhandle[1]);

		//face_vhandles.push_back(vhandle[2]);

		//face_vhandles.push_back(vhandle[3]);

		//tMesh.add_face(face_vhandles);



		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[7]);

		//face_vhandles.push_back(vhandle[6]);

		//face_vhandles.push_back(vhandle[5]);

		//face_vhandles.push_back(vhandle[4]);

		//tMesh.add_face(face_vhandles);



		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[1]);

		//face_vhandles.push_back(vhandle[0]);

		//face_vhandles.push_back(vhandle[4]);

		//face_vhandles.push_back(vhandle[5]);

		//tMesh.add_face(face_vhandles);



		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[2]);

		//face_vhandles.push_back(vhandle[1]);

		//face_vhandles.push_back(vhandle[5]);

		//face_vhandles.push_back(vhandle[6]);

		//tMesh.add_face(face_vhandles);



		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[3]);

		//face_vhandles.push_back(vhandle[2]);

		//face_vhandles.push_back(vhandle[6]);

		//face_vhandles.push_back(vhandle[7]);

		//tMesh.add_face(face_vhandles);


		//face_vhandles.clear();

		//face_vhandles.push_back(vhandle[0]);

		//face_vhandles.push_back(vhandle[3]);

		//face_vhandles.push_back(vhandle[7]);

		//face_vhandles.push_back(vhandle[4]);

		//tMesh.add_face(face_vhandles);

		//tMesh.request_edge_status();
		//tMesh.request_vertex_status();
		//tMesh.request_face_status();



		/*
		//cout <<"==================" << endl;

		//{
		//	Mesh::FaceHandle fh = tMesh.face_handle(0);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	cout << "Face0 p1: " << p1 << endl;
		//	cout << "Face0 p2: " << p2 << endl;
		//	cout << "Face0 p3: " << p3 << endl;
		//}

		//{
		//	Mesh::FaceHandle fh = tMesh.face_handle(1);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	cout << "Face1 p1: " << p1 << endl;
		//	cout << "Face1 p2: " << p2 << endl;
		//	cout << "Face1 p3: " << p3 << endl;
		//}

		//{
		//	Mesh::FaceHandle fh = tMesh.face_handle(2);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	cout << "Face2 p1: " << p1 << endl;
		//	cout << "Face2 p2: " << p2 << endl;
		//	cout << "Face2 p3: " << p3 << endl;
		//}


		//{
		//	cout << "=====================" << endl;

		//	Mesh::FaceHandle fh = *tMesh.faces_begin();

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	Mesh::Point pC = (p1 + p2 + p3) / 3;

		//	Mesh::VertexHandle va = tMesh.add_vertex(pC);

		//	tMesh.delete_face(fh, false);


		//	face_vhandles.clear();
		//	face_vhandles.push_back(v0);
		//	face_vhandles.push_back(v1);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace1 = tMesh.add_face(face_vhandles);

		//	face_vhandles.clear();
		//	face_vhandles.push_back(v2);
		//	face_vhandles.push_back(v0);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace2 = tMesh.add_face(face_vhandles);

		//	face_vhandles.clear();
		//	face_vhandles.push_back(v1);
		//	face_vhandles.push_back(v2);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace3 = tMesh.add_face(face_vhandles);

		//	int j = 0;
		//	if (tMesh.is_valid_handle(newFace1)) j++;
		//	if (tMesh.is_valid_handle(newFace2)) j++;
		//	if (tMesh.is_valid_handle(newFace3)) j++;

		//	cout << newFace1.idx() << endl;
		//	cout << newFace2.idx() << endl;
		//	cout << newFace3.idx() << endl;

		//	cout << j << endl;

		//	{
		//		cout << "new face point 1" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace1);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}

		//	{
		//		cout << "new face point 2" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace2);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}

		//	{
		//		cout << "new face point 3" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace3);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}
		//}


		//{
		//	cout << "=====================" << endl;

		//	Mesh::FaceHandle fh = tMesh.face_handle(2);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	Mesh::Point pC = (p1 + p2 + p3) / 3;

		//	Mesh::VertexHandle va = tMesh.add_vertex(pC);

		//	tMesh.delete_face(fh, false);


		//	face_vhandles.clear();
		//	face_vhandles.push_back(v0);
		//	face_vhandles.push_back(v1);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace1 = tMesh.add_face(face_vhandles);

		//	face_vhandles.clear();
		//	face_vhandles.push_back(v2);
		//	face_vhandles.push_back(v0);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace2 = tMesh.add_face(face_vhandles);

		//	face_vhandles.clear();
		//	face_vhandles.push_back(v1);
		//	face_vhandles.push_back(v2);
		//	face_vhandles.push_back(va);
		//	Mesh::FaceHandle newFace3 = tMesh.add_face(face_vhandles);

		//	int j = 0;
		//	if (tMesh.is_valid_handle(newFace1)) j++;
		//	if (tMesh.is_valid_handle(newFace2)) j++;
		//	if (tMesh.is_valid_handle(newFace3)) j++;

		//	cout << newFace1.idx() << endl;
		//	cout << newFace2.idx() << endl;
		//	cout << newFace3.idx() << endl;

		//	cout << j << endl;

		//	{
		//		cout << "new face point 1" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace1);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}

		//	{
		//		cout << "new face point 2" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace2);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}

		//	{
		//		cout << "new face point 3" << endl;
		//		Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(newFace3);

		//		Mesh::VertexHandle v0 = *(fv++);
		//		Mesh::VertexHandle v1 = *(fv++);
		//		Mesh::VertexHandle v2 = *(fv++);

		//		Mesh::Point p1 = tMesh.point(v0);
		//		Mesh::Point p2 = tMesh.point(v1);
		//		Mesh::Point p3 = tMesh.point(v2);

		//		cout << "Face0 p1: " << p1 << endl;
		//		cout << "Face0 p2: " << p2 << endl;
		//		cout << "Face0 p3: " << p3 << endl;
		//	}
		//}


		//cout << "===============" << endl;

		//for (int i = 0;i < tMesh.n_faces();i++) {
		//
		//	Mesh::FaceHandle fh = tMesh.face_handle(i);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	cout << "Face" << i << " p1: " << p1 << endl;
		//	cout << "Face" << i << " p2: " << p2 << endl;
		//	cout << "Face" << i << " p3: " << p3 << endl;

		//	cout  << endl;
		//}

		//tMesh.garbage_collection();

		//cout << "****************" << endl;
		//for (int i = 0;i < tMesh.n_faces();i++) {

		//	Mesh::FaceHandle fh = tMesh.face_handle(i);

		//	Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//	Mesh::VertexHandle v0 = *(fv++);
		//	Mesh::VertexHandle v1 = *(fv++);
		//	Mesh::VertexHandle v2 = *(fv++);

		//	Mesh::Point p1 = tMesh.point(v0);
		//	Mesh::Point p2 = tMesh.point(v1);
		//	Mesh::Point p3 = tMesh.point(v2);

		//	cout << "Face" << i << " p1: " << p1 << endl;
		//	cout << "Face" << i << " p2: " << p2 << endl;
		//	cout << "Face" << i << " p3: " << p3 << endl;

		//	cout << endl;
		//}

		//cout << "===============" << endl;
		*/

		//Mesh::FaceHandle fh = *tMesh.faces_begin();

		//Mesh::FaceVertexCCWIter fv = tMesh.fv_ccwbegin(fh);

		//Mesh::VertexHandle v0 = *(fv++);
		//Mesh::VertexHandle v1 = *(fv++);

		//Mesh::Point p1 = tMesh.point(v0);
		//Mesh::Point p2 = tMesh.point(v1);

		//Mesh::Point pC = (p1 + p2) / 2;

		//for (Mesh::FaceIter fi = tMesh.faces_begin(); fi != tMesh.faces_end(); ++fi) {
		//	cout << (*fi).idx() << endl;
		//}

		//Mesh::EdgeHandle eh = tMesh.edge_handle(tMesh.find_halfedge(v0, v1));
		//Mesh::VertexHandle va = tMesh.add_vertex(pC);
		// 
		//tMesh.split_edge(eh, va);

		//for (Mesh::FaceIter fi = tMesh.faces_begin(); fi != tMesh.faces_end(); ++fi) {
		//	cout << (*fi).idx() << endl;
		//}

		//Mesh::EdgeHandle eh = *mesh.edges_begin();
		//Mesh::HalfedgeHandle he1 = mesh.halfedge_handle(eh, 0);
		//Mesh::HalfedgeHandle he2 = mesh.halfedge_handle(eh, 1);

		//cout << mesh.point(mesh.from_vertex_handle(he1)) << endl;
		//cout << mesh.point(mesh.to_vertex_handle(he1)) << endl;

		//cout << mesh.point(mesh.from_vertex_handle(he2)) << endl;
		//cout << mesh.point(mesh.to_vertex_handle(he2)) << endl;

		//mesh.flip(eh);

		//cout << mesh.point(mesh.from_vertex_handle(he1)) << endl;
		//cout << mesh.point(mesh.to_vertex_handle(he1)) << endl;

		//cout << mesh.point(mesh.from_vertex_handle(he2)) << endl;
		//cout << mesh.point(mesh.to_vertex_handle(he2)) << endl;




		//// write mesh to output.obj

		//try

		//{

		//	if (!OpenMesh::IO::write_mesh(tMesh, "output.off"))

		//	{

		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;

		//		return;

		//	}

		//}

		//catch (std::exception& x)

		//{

		//	std::cerr << x.what() << std::endl;

		//	return;

		//}



		return;
	}

	void TestMath()
	{
		Vector3d v(0.7, 0.7, 0.7);
		
		auto v2 = v.normalized();
		cout << v.length() << endl;
		cout << v2.length() << endl;

		v.normalize();
		cout << v.length() << endl;

		
	}

	void TestBspTree()
	{
		//BatchTestClassicModels();
		//BatchTestThingi10K();

		SingleTestClassicModelsForBSP();
		//SingleTestClassicModelsForBSP(23);
		//SingleTestClassicModelsForBSP(25);
		//SingleTestClassicModelsForBSP(28);



		//SingleTestThingi10K();

		//

		//{
		//	auto* ra = new TrianglesRenderer<Triangle_Set>(planeTris);
		//	ra->SetTriangleColor(glm::vec3(0.0f, 1.0f, 0.0f));
		//	scence.AddUnit(ra);
		//}

		//{
		//	auto* ra = new EdgesRenderer<Edge_Set>(planeEdges);

		//	ra->SetEdgeWidth(10);
		//	ra->SetEdgeColor(glm::vec3(1.0f, 1.0f, 0.0f));
		//	scence.AddUnit(ra);
		//}

		//if (!OpenMesh::IO::write_mesh(planeMesh, outputAddress)) {
		//	std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//	return;
		//}








		//if (false) {
		//	//显示meshA
		//	if (true)
		//	{

		//		Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshA);
		//		Edge_Set thisEdge = DataTranslate::TranslateEdge(meshA);
		//		Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshA);

		//		{
		//			auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
		//			ra->SetTriangleColor(glm::vec3(0.0f, 0.0f, 1.0f));
		//			scence.AddUnit(ra);
		//		}


		//		{
		//			auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
		//			ra->SetVertexSize(15.0);
		//			ra->SetVertexColor(glm::vec3(1.0f, 0.0f, 0.0f));
		//			scence.AddUnit(ra);
		//		}

		//		{
		//			auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

		//			//ra->SetVertexSize(15.0);
		//			ra->SetEdgeColor(glm::vec3(1.0f, 1.0f, 0.0f));
		//			scence.AddUnit(ra);
		//		}
		//	}

		//	//显示meshB
		//	if (true)
		//	{

		//		Triangle_Set thisMesh = DataTranslate::TranslateMesh(meshB);
		//		Edge_Set thisEdge = DataTranslate::TranslateEdge(meshB);
		//		Vertex_Set thisVertex = DataTranslate::TranslateVertex(meshB);

		//		{
		//			auto* ra = new TrianglesRenderer<Triangle_Set>(thisMesh);
		//			ra->SetTriangleColor(glm::vec3(0.0f, 0.0f, 1.0f));
		//			scence.AddUnit(ra);
		//		}


		//		{
		//			auto* ra = new VerticesRenderer<Vertex_Set>(thisVertex);
		//			ra->SetVertexSize(15.0);
		//			ra->SetVertexColor(glm::vec3(0.0f, 1.0f, 0.0f));
		//			scence.AddUnit(ra);
		//		}

		//		{
		//			auto* ra = new EdgesRenderer<Edge_Set>(thisEdge);

		//			//ra->SetVertexSize(15.0);
		//			ra->SetEdgeColor(glm::vec3(1.0f, 1.0f, 0.0f));
		//			scence.AddUnit(ra);
		//		}
		//	}
		//}




	}

	void SingleTestClassicModelsForBSP()
	{

		Mesh meshA, meshB;

		string testResult_Num = "28";
		string fileType = "obj";
		string outputPartitionPlaneAddress1 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane1.obj";
		string outputPartitionPlaneAddress2 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane2.obj";
		string outputPartitionPlaneAddress3 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane3.obj";
		string outputPartitionPlaneAddress4 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane4.obj";


		string outputPartitionModelAddress = "..//TestModel//bunny//" + testResult_Num + "//TestTreeResult//";


		if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bunny//" + testResult_Num + "//bunnyA." + fileType)) {
		//if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bool operate//" + testResult_Num + "//bunnyA." + fileType)) {
			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}
		if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bunny//" + testResult_Num + "//bunnyB." + fileType)) {
		//if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bool operate//" + testResult_Num + "//bunnyB." + fileType)) {
			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}



		Tolerance toler;

		{
			meshA.request_vertex_status();
			meshA.request_edge_status();
			meshA.request_face_status();

			meshB.request_vertex_status();
			meshB.request_edge_status();
			meshB.request_face_status();

			meshA.request_face_normals();
			meshB.request_face_normals();

			meshA.update_face_normals();
			meshB.update_face_normals();
		}


		PreProcess pr({ &meshA,&meshB }, toler);
		pr.run();



		//for (int i = 8; i < 9;i++) {

		//	bool needOutputMesh = false;
		//	if ((8 == i) || (10 == i) || (15 == i)) needOutputMesh = true;

		//	bool needITRA = false;
		//	if ((8 == i) || (10 == i) || (15 == i)) needITRA = true;

		//	double SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::SDM);
		//	double SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::SAH);
		//	double ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::ObbMiddel);
		//	double Gravity = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::Gravity_SPLIT);

		//	//cout << "========" << endl;
		//	//cout << "SDMUTFR:     " << SDMUTFR << endl;
		//	//cout << "SAHUTFR:     " << SAHUTFR << endl;
		//	//cout << "ObbMiddelUTFR:  " << ObbMiddel << endl;
		//	//cout << "GravityUTFR: " << Gravity << endl;
		//	//cout << "========" << endl;		
		//}
		
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::SDM);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::ObbMiddel);

		//for (int i = 3;i >= 1;i--) {
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::SDM);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::ObbMiddel);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::SAH);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::Gravity_SPLIT);
		//}

		//TestBspTimeConsum(meshA, meshB, toler, 1, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		//TestBspTimeConsum(meshA, meshB, toler, 2, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);
		
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::SDM);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::ObbMiddel);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::Gravity_SPLIT);


		auto SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		auto SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
		auto ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
		auto Gravity = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);

		//for (int i = 4; i < 16;i++) {
		//	if (5 == i) continue;
		//	if (7 == i) continue;

		//	bool neadwrite = false;
		//	if (10 == i) {
		//		neadwrite = true;
		//	}

		//	auto SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, neadwrite, true, BSPConstructType::SDM);
		//	auto SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, neadwrite, true, BSPConstructType::SAH);
		//	auto ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, neadwrite, true, BSPConstructType::ObbMiddel);
		//	auto Gravity = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, neadwrite, true, BSPConstructType::Gravity_SPLIT);
		//
		//}



		//cout << "========" << endl;
		//cout << "SDMUTFR:     " << SDMUTFR.first << endl;
		//cout << "SAHUTFR:     " << SAHUTFR.first << endl;
		//cout << "ObbMiddelUTFR:  " << ObbMiddel.first << endl;
		//cout << "GravityUTFR: " << Gravity.first << endl;
		//cout << "SDMITRA:     " << SDMUTFR.second << endl;
		//cout << "SAHITRA:     " << SAHUTFR.second << endl;
		//cout << "ObbMiddelITRA:  " << ObbMiddel.second << endl;
		//cout << "GravityITRA: " << Gravity.second << endl;
		//cout << "========" << endl;	

	
		//////显示分割平面	 
		//{
		//	BSPTree* faceTree = new BSPTree(meshA, meshB, toler, 2, 25, BSPConstructType::SDM);
		//	bspFaceRenderingInfo partiPlanes;
		//	faceTree->GetPartitionPlane(partiPlanes);

		//	Mesh planeMesh;
		//	Triangle_Set planeTris;
		//	Edge_Set planeEdges;
		//	tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

		//	if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress1)) {
		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//		return;
		//	}
		//}


		//{
		//	BSPTree* faceTree1 = new BSPTree(meshA, meshB, toler, 1, 25, BSPConstructType::ObbMiddel);
		//	//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, BSPConstructType::Gravity_SPLIT);

		//	bspFaceRenderingInfo partiPlanes;
		//	faceTree1->GetPartitionPlane(partiPlanes);

		//	Mesh planeMesh;
		//	Triangle_Set planeTris;
		//	Edge_Set planeEdges;
		//	tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

		//	if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress1)) {
		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//		return;
		//	}


			////显示子树划分出的三角网格
			//vector<BSPTreeNode*> allLeafNodes;
			//faceTree1->GetAllLeafNode(allLeafNodes);

			//cout << meshA.n_faces() << endl;
			//cout << meshB.n_faces() << endl;

			//int ii = 0;
			//for (BSPTreeNode* noded : allLeafNodes) {
			//	ii++;
			//	Mesh tempMeshA, tempMeshB;
			//	tie(tempMeshA, tempMeshB) = DataTranslate::TranslateNode2SeperatedTriangleMesh(noded);

			//	cout << tempMeshA.n_faces() << endl;
			//	cout << tempMeshB.n_faces() << endl;

			//	if (!OpenMesh::IO::write_mesh(tempMeshA, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) + "subMeshA.obj")) {
			//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			//		return;
			//	}
			//	if (!OpenMesh::IO::write_mesh(tempMeshB, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) + "subMeshB.obj")) {
			//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			//		return;
			//	}

			//}


		//	delete faceTree1;
		//}

	/*	{
			BSPTree* faceTree2 = new BSPTree(meshA, meshB, toler, 1, 5, BSPConstructType::SAH);

			bspFaceRenderingInfo partiPlanes;
			faceTree2->GetPartitionPlane(partiPlanes);

			Mesh planeMesh;
			Triangle_Set planeTris;
			Edge_Set planeEdges;
			tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

			if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress2)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return;
			}
			delete faceTree2;
		}
		{
			BSPTree* faceTree3 = new BSPTree(meshA, meshB, toler, 1, 5, BSPConstructType::ObbMiddel);

			bspFaceRenderingInfo partiPlanes;
			faceTree3->GetPartitionPlane(partiPlanes);

			Mesh planeMesh;
			Triangle_Set planeTris;
			Edge_Set planeEdges;
			tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

			if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress3)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return;
			}
			delete faceTree3;
		}
		{
			BSPTree* faceTree4 = new BSPTree(meshA, meshB, toler, 1, 5, BSPConstructType::Gravity_SPLIT);


			bspFaceRenderingInfo partiPlanes;
			faceTree4->GetPartitionPlane(partiPlanes);

			Mesh planeMesh;
			Triangle_Set planeTris;
			Edge_Set planeEdges;
			tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

			if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress4)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return;
			}

			delete faceTree4;
		} */


		////显示分割后网格
		//{
		//	//显示子树划分出的三角网格
		//	vector<BSPTreeNode*> allLeafNodes;
		//	faceTree1->GetAllLeafNode(allLeafNodes);
		//	
		//	int ii = 0;
		//	for (BSPTreeNode* noded : allLeafNodes) {
		//		ii++;
		//		Mesh tempMeshA, tempMeshB;
		//		tie(tempMeshA,tempMeshB) = DataTranslate::TranslateNode2SeperatedTriangleMesh(noded);

		//		if (!OpenMesh::IO::write_mesh(tempMeshA, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) +"subMeshA.obj")) {
		//			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//			return;
		//		}
		//		if (!OpenMesh::IO::write_mesh(tempMeshB, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) + "subMeshB.obj")) {
		//			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//			return;
		//		}

		//	}

		//}



	}

	void SingleTestClassicModelsForBSP(int ID)
	{

		Mesh meshA, meshB;

		string testResult_Num = to_string(ID);
		string fileType = "obj";
		string outputPartitionPlaneAddress1 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane1.obj";
		string outputPartitionPlaneAddress2 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane2.obj";
		string outputPartitionPlaneAddress3 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane3.obj";
		string outputPartitionPlaneAddress4 = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane4.obj";


		string outputPartitionModelAddress = "..//TestModel//bunny//" + testResult_Num + "//TestTreeResult//";


		if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bunny//" + testResult_Num + "//bunnyA." + fileType)) {
			//if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bool operate//" + testResult_Num + "//bunnyA." + fileType)) {
			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}
		if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bunny//" + testResult_Num + "//bunnyB." + fileType)) {
			//if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bool operate//" + testResult_Num + "//bunnyB." + fileType)) {
			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}



		Tolerance toler;

		{
			meshA.request_vertex_status();
			meshA.request_edge_status();
			meshA.request_face_status();

			meshB.request_vertex_status();
			meshB.request_edge_status();
			meshB.request_face_status();

			meshA.request_face_normals();
			meshB.request_face_normals();

			meshA.update_face_normals();
			meshB.update_face_normals();
		}


		PreProcess pr({ &meshA,&meshB }, toler);
		pr.run();



		ofstream of("..//TestModel//bunny//" + to_string(ID) + "//time.txt",ios::app);


		for (int i = 4; i < 16;i++) {

			of << "========" << endl;
			of << "SDM: SAH: ObbMiddel: Gravity: " << endl;

			if (i == 5) continue;
			if (i == 7) continue;
			if (i == 9) continue;
			if (i == 11) continue;
			if (i == 13) continue;

			bool needOutputMesh = false;

			//if ((8 == i) || (10 == i) || (15 == i)) needOutputMesh = true;

			bool needITRA = true;

			//if ((8 == i) || (10 == i) || (15 == i)) needITRA = true;

			auto SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::SDM, ID);
			auto SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::SAH, ID);
			auto ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::ObbMiddel, ID);
			auto Gravity = TestBspTimeConsum(meshA, meshB, toler, i, 25, outputPartitionModelAddress, needOutputMesh, needITRA, BSPConstructType::Gravity_SPLIT, ID);


			of << "SDM: SAH: ObbMiddel: Gravity: " << endl;
			of << SDMUTFR.first << "  " << SDMUTFR.second << endl;
			of << SAHUTFR.first << "  " << SAHUTFR.second << endl; 
			of << ObbMiddel.first << "  " << ObbMiddel.second << endl;
			of << Gravity.first << "  " << Gravity.second << endl;
			of << "========" << endl;
		}


		of.close();

		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::SDM);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::ObbMiddel);

		//for (int i = 3;i >= 1;i--) {
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::SDM);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::ObbMiddel);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::SAH);
			//TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, true, true, BSPConstructType::Gravity_SPLIT);
		//}

		//TestBspTimeConsum(meshA, meshB, toler, 1, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		//TestBspTimeConsum(meshA, meshB, toler, 2, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);

		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
		//TestBspTimeConsum(meshA, meshB, toler, 3, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);

		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::SDM);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::ObbMiddel);
		//TestBspNodeDistributionInfo(meshA, meshB, toler, 15, 5, outputPartitionModelAddress, true, BSPConstructType::Gravity_SPLIT);


		//auto SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		//auto SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
		//auto ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
		//auto Gravity = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);


	}

	void SingleTestThingi10K()
	{

		int num = 20;
		Mesh meshA, meshB;

		string baseAddress = "..//TestModel//thingi10k//" + to_string(num);

		vector<string>  files;
		getFileNames(baseAddress, files);

		if (!OpenMesh::IO::read_mesh(meshA, files[0])) {

			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}
		if (!OpenMesh::IO::read_mesh(meshB, files[1])) {

			std::cerr << "Cannot read mesh" << std::endl;

			return;

		}




		Tolerance toler;

		{
			meshA.request_vertex_status();
			meshA.request_edge_status();
			meshA.request_face_status();

			meshB.request_vertex_status();
			meshB.request_edge_status();
			meshB.request_face_status();

			meshA.request_face_normals();
			meshB.request_face_normals();

			meshA.update_face_normals();
			meshB.update_face_normals();
		}


		PreProcess pr({ &meshA,&meshB }, toler);
		pr.run();

		string outputPartitionModelAddress = "";

		auto SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);
		auto SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
		auto ObbMiddel = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
		auto Gravity = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);

		cout << "========" << endl;
		cout << "SDMUTFR:     " << SDMUTFR.first << endl;
		cout << "SAHUTFR:     " << SAHUTFR.first << endl;
		cout << "ObbMiddelUTFR:  " << ObbMiddel.first << endl;
		cout << "GravityUTFR: " << Gravity.first << endl;
		cout << "========" << endl;

		//BSPTree* faceTree1 = new BSPTree(meshA, meshB, toler, 1, 5, BSPConstructType::SDM);
		//BSPTree* faceTree2 = new BSPTree(meshA, meshB, toler, 1, 5, BSPConstructType::ObbMiddel);

		////显示分割平面	 
		//{
		//	bspFaceRenderingInfo partiPlanes;
		//	faceTree1->GetPartitionPlane(partiPlanes);

		//	Mesh planeMesh;
		//	Triangle_Set planeTris;
		//	Edge_Set planeEdges;
		//	tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

		//	if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress1)) {
		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//		return;
		//	}
		//}

		////显示分割后网格
		//{
		//	//显示子树划分出的三角网格
		//	vector<BSPTreeNode*> allLeafNodes;
		//	faceTree1->GetAllLeafNode(allLeafNodes);
		//	
		//	int ii = 0;
		//	for (BSPTreeNode* noded : allLeafNodes) {
		//		ii++;
		//		Mesh tempMeshA, tempMeshB;
		//		tie(tempMeshA,tempMeshB) = DataTranslate::TranslateNode2SeperatedTriangleMesh(noded);

		//		if (!OpenMesh::IO::write_mesh(tempMeshA, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) +"subMeshA.obj")) {
		//			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//			return;
		//		}
		//		if (!OpenMesh::IO::write_mesh(tempMeshB, "..//TestModel//bunny//" + testResult_Num + "//" + to_string(ii) + "subMeshB.obj")) {
		//			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//			return;
		//		}

		//	}

		//}


		//{
		//	bspFaceRenderingInfo partiPlanes;
		//	faceTree2->GetPartitionPlane(partiPlanes);

		//	Mesh planeMesh;
		//	Triangle_Set planeTris;
		//	Edge_Set planeEdges;
		//	tie(planeTris, planeEdges) = DataTranslate::TranslatePartitionPlane(partiPlanes, planeMesh);

		//	if (!OpenMesh::IO::write_mesh(planeMesh, outputPartitionPlaneAddress2)) {
		//		std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		//		return;
		//	}
		//}

	}

	void BatchTestClassicModels()
	{
		Mesh meshA, meshB;

		//ofstream outFile("..//TestModel//bunny//res.txt");

		bool needWriteResult = true;

		for (int i = 7;i < 14;i++) {

			string testResult_Num = to_string(i);
			string outputPartitionPlaneAddress = "..//TestModel//bunny//" + testResult_Num + "//PartitionPlane.obj";
			string outputPartitionModelAddress = "..//TestModel//bunny//" + testResult_Num + "//TestTreeResult//";

			if (!OpenMesh::IO::read_mesh(meshA, "..//TestModel//bunny//" + to_string(i) + "//bunnyA.obj")) {

				std::cerr << "Cannot read mesh" << std::endl;

				return;
				
			}
			if (!OpenMesh::IO::read_mesh(meshB, "..//TestModel//bunny//" + to_string(i) + "//bunnyB.obj")) {

				std::cerr << "Cannot read mesh" << std::endl;

				return;

			}

			Tolerance toler;

			{
				meshA.request_vertex_status();
				meshA.request_edge_status();
				meshA.request_face_status();

				meshB.request_vertex_status();
				meshB.request_edge_status();
				meshB.request_face_status();

				meshA.request_face_normals();
				meshB.request_face_normals();

				meshA.update_face_normals();
				meshB.update_face_normals();
			}


			PreProcess pr({ &meshA,&meshB }, toler);
			pr.run();

			cout << "classic:========= " << i << endl;
			auto SDMITRA = TestBspNodeDistributionInfo(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, BSPConstructType::SDM);

			//double SDMUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 5, outputPartitionModelAddress, needWriteResult, true, BSPConstructType::SDM);
			//double SAHUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 5, outputPartitionModelAddress, needWriteResult, true, BSPConstructType::SAH);
			//double ObbMiddelUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 5, outputPartitionModelAddress, needWriteResult, true, BSPConstructType::ObbMiddel);
			//double GravityUTFR = TestBspTimeConsum(meshA, meshB, toler, 10, 5, outputPartitionModelAddress, needWriteResult, true, BSPConstructType::Gravity_SPLIT);

			//outFile << "========" << endl;
			//outFile << "model num:   " << i << endl;
			//outFile << "SDMUTFR:     " << SDMUTFR << endl;
			//outFile << "SAHUTFR:     " << SAHUTFR << endl;
			//outFile << "ObbMiddelUTFR:  " << ObbMiddelUTFR << endl;
			//outFile << "GravityUTFR: " << GravityUTFR << endl;
			//outFile << "========" << endl;
		}


		//outFile.close();
	}

	void BatchTestThingi10K()
	{

		string outputPartitionPlaneAddress = "..//TestModel//thingi10k//PartitionPlane.obj";
		string outputPartitionModelAddress = "..//TestModel//thingi10k//TestTreeResult//";

		//ofstream outFile("..//TestModel//thingi10k//res.txt");
		
		for (int i = 21;i < 51;i++) {

			//if (i == 17)
			//	continue;

			cout << "thingi10k:============== " << i << endl;
			Mesh meshA, meshB;

			string baseAddress = "..//TestModel//thingi10k//" + to_string(i);

			vector<string>  files;
			getFileNames(baseAddress, files);

			if (!OpenMesh::IO::read_mesh(meshA, files[0])) {

				std::cerr << "Cannot read mesh" << std::endl;

				return;

			}
			if (!OpenMesh::IO::read_mesh(meshB, files[1])) {

				std::cerr << "Cannot read mesh" << std::endl;

				return;

			}

			Tolerance toler;

			{
				meshA.request_vertex_status();
				meshA.request_edge_status();
				meshA.request_face_status();

				meshB.request_vertex_status();
				meshB.request_edge_status();
				meshB.request_face_status();

				meshA.request_face_normals();
				meshB.request_face_normals();

				meshA.update_face_normals();
				meshB.update_face_normals();
			}


			PreProcess pr({ &meshA,&meshB }, toler);
			pr.run();

			//auto SDMITRA = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SDM);


			auto SDMITRA = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, false, BSPConstructType::SDM);
			//auto SAHITRA = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::SAH);
			//auto ObbMiddelITRA = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::ObbMiddel);
			//auto GravityITRA = TestBspTimeConsum(meshA, meshB, toler, 10, 25, outputPartitionModelAddress, false, true, BSPConstructType::Gravity_SPLIT);

			//outFile << "========" << endl;
			//outFile << "model num:   " << i << endl;

			//outFile << "SDMNTRR:     " << SDMITRA.first << endl;
			//outFile << "SAHNTRR:     " << SAHITRA.first << endl;
			//outFile << "ObbMiddelNTRR:  " << ObbMiddelITRA.first << endl;
			//outFile << "GravityNTRR: " << GravityITRA.first << endl;

			//outFile << "SDMITRA:     " << SDMITRA.second << endl;
			//outFile << "SAHITRA:     " << SAHITRA.second << endl;
			//outFile << "ObbMiddelITRA:  " << ObbMiddelITRA.second << endl;
			//outFile << "GravityITRA: " << GravityITRA.second << endl;

			//outFile << "========" << endl;
		}

		//outFile.close();


	}

	pair<double, double> TestBspTimeConsum(Mesh& meshA, Mesh& meshB, Tolerance& toler, int maxDepth, int leafShreshold, string p_BaseAddress, bool outPutResultMesh, bool outPutITRA, BSPConstructType bspType /*= BSPConstructType::AABB*/, int Id)
	{
		double r_UTFR = -1;
		double r_ITRA = -1;

		string baseAddress = p_BaseAddress + to_string(maxDepth) + "//";
		cout << "maxDepth " << maxDepth << endl;
		cout << "leafShreshold " << leafShreshold << endl;

		switch (bspType)
		{
		case BSPConstructType::AABB_MIDDLE_SPLIT:
			cout << "BSPConstructType " << "AABB"  << endl;
			break;
		case BSPConstructType::SDM:
			cout << "BSPConstructType " << "SDM" << endl;
			break;
		case BSPConstructType::SAH:
			cout << "BSPConstructType " << "SAH" << endl;
			break;
		case BSPConstructType::ObbMiddel:
			cout << "BSPConstructType " << "ObbMiddel" << endl;
			break;
		case BSPConstructType::Gravity_SPLIT:
			cout << "BSPConstructType " << "Gravity_SPLIT" << endl;
			break;
		default:
			break;
		}


		auto tBuildTreeBegin = std::chrono::steady_clock::now();

		BSPTree* faceTree = new BSPTree(meshA, meshB, toler, maxDepth, leafShreshold, bspType);

		auto tBuildTreeEnd = std::chrono::steady_clock::now();
		cout << "建树耗时： " << std::chrono::duration<double, std::milli>(tBuildTreeEnd - tBuildTreeBegin).count() << "毫秒" << endl;


		//显示子树划分出的三角网格
		vector<BSPTreeNode*> allLeafNodes;
		faceTree->GetAllLeafNode(allLeafNodes);



		Mesh resultValidMeshA, resultValidMeshB;

		Mesh realIntersectFacesMeshA, realIntersectFacesMeshB;
		Mesh fakeIntersectFacesMeshA, fakeIntersectFacesMeshB;
		Mesh resultInValidMeshA, resultInValidMeshB;
		{
			IntersectTriangleCheckList intersectCheckList;

			// lambda for user-defined hash function
			auto hash = [](const pair<FaceId, FaceId>& c) {
				return hash_val(c.first, c.second);
			};

			// lambda for user-defined equality criterion
			auto eq = [](const pair<FaceId, FaceId>& c1, const pair<FaceId, FaceId>& c2) {
				return c1 == c2;
			};

			// create unordered set with user-defined behavior
			std::unordered_set<pair<FaceId, FaceId>, decltype(hash), decltype(eq)> nowProcessFacePairsSet(intersectCheckList.size(), hash, eq);

			//std::unordered_set<pair<EdgeId, FaceId>, decltype(hash), decltype(eq)> nowProcessEdgeAFaceBPairSet(intersectCheckList.size() * 2, hash, eq);

			//std::unordered_set<pair<EdgeId, FaceId>, decltype(hash), decltype(eq)> nowProcessEdgeBFaceAPairSet(intersectCheckList.size() * 2, hash, eq);


			auto t1 = std::chrono::steady_clock::now();

			//for (int i = 0;(i < 5) && (nowNodes != m_AllLeafNodes.end()); i++, nowNodes++) {

			//cout << "nodes count: " << allLeafNodes.size() << endl;

			for (auto nowNodes = allLeafNodes.begin(); nowNodes != allLeafNodes.end(); nowNodes++) {
				(*nowNodes)->GetIntersectPair(intersectCheckList);
			}

			//auto t2 = std::chrono::steady_clock::now();
			//std::cout << "获取初步数据时间：" << std::chrono::duration<double, std::milli>(t2 - t1).count() << "毫秒" << endl;

			//先对采集面片对去重
			for (auto& checkPair : intersectCheckList) {
				nowProcessFacePairsSet.insert(checkPair);
			}

			auto t3 = std::chrono::steady_clock::now();
			std::cout << "查询时间：" << std::chrono::duration<double, std::milli>(t3 - t1).count() << "毫秒" << endl;


			//test/////////////////////////////////////////
			if (Id >= 0) {
				ofstream of("..//TestModel//bunny//" + to_string(Id) + "//time.txt", ios::app);

				of << "建树时间， 查询时间" << endl;
				of << std::chrono::duration<double, std::milli>(tBuildTreeEnd - tBuildTreeBegin).count() << "   ";
				of << std::chrono::duration<double, std::milli>(t3 - t1).count() << endl;
				of.close();		
			}

			//test/////////////////////////////////////////


			int validCount = 0;
			unordered_set<FaceId> realIntersectFacesA, fakeIntersectFacesA;
			unordered_set<FaceId> realIntersectFacesB, fakeIntersectFacesB;

			if (outPutITRA) {
				for (auto& verifyPair : nowProcessFacePairsSet) {

					Mesh::FaceVertexCCWIter fvit = meshA.fv_ccwbegin(meshA.face_handle(verifyPair.first));
					Mesh::VertexHandle triV0 = *fvit;
					Mesh::VertexHandle triV1 = *(++fvit);
					Mesh::VertexHandle triV2 = *(++fvit);

					Triangle3d triMathA(array<Point3d, 3>{
						Point3d(meshA.point(triV0)),
							Point3d(meshA.point(triV1)),
							Point3d(meshA.point(triV2))
					});

					Mesh::FaceVertexCCWIter fvitB = meshB.fv_ccwbegin(meshB.face_handle(verifyPair.second));
					Mesh::VertexHandle triV0B = *fvitB;
					Mesh::VertexHandle triV1B = *(++fvitB);
					Mesh::VertexHandle triV2B = *(++fvitB);

					Triangle3d triMathB(array<Point3d, 3>{
						Point3d(meshB.point(triV0B)),
							Point3d(meshB.point(triV1B)),
							Point3d(meshB.point(triV2B))
					});


					if (TriangleTriangleIsIntersect(triMathA, triMathB)) {
						validCount++;
						realIntersectFacesA.insert(verifyPair.first);
						realIntersectFacesB.insert(verifyPair.second);
					}
					else {
						fakeIntersectFacesA.insert(verifyPair.first);
						fakeIntersectFacesB.insert(verifyPair.second);

					}
				}

			}

			for (auto& fa : realIntersectFacesA) {
				fakeIntersectFacesA.erase(fa);
			}

			for (auto& fa : realIntersectFacesB) {
				fakeIntersectFacesB.erase(fa);
			}

			unordered_set<FaceId> valid_FacesA, valid_FacesB;
			unordered_set<FaceId> inValid_FacesA, inValid_FacesB;

			for (auto nowNodes = allLeafNodes.begin(); nowNodes != allLeafNodes.end(); nowNodes++) {
				(*nowNodes)->GetValidFaces(valid_FacesA, valid_FacesB);
			}

			r_UTFR = (1 - double(valid_FacesA.size() + valid_FacesB.size()) / double(meshA.n_faces() + meshB.n_faces())) * 100;
			cout << "UTFR :" << r_UTFR << endl;
			
			cout << "new UTFR1: " << (1 - double(valid_FacesA.size()) / double(meshA.n_faces())) * (1 - double(valid_FacesB.size()) / double(meshB.n_faces())) << endl;
			cout << "new UTFR2: " << 1 - (double(valid_FacesA.size()) * double(valid_FacesB.size())) / ( double(meshA.n_faces())  * double( meshB.n_faces() ) ) << endl;
			cout << "new UTFR3: " << (double(valid_FacesA.size()) * double(valid_FacesB.size())) / (double(meshA.n_faces()) * double(meshB.n_faces())) << endl;
			cout << "real intersect: " << double(nowProcessFacePairsSet.size()) / double(meshA.n_faces() * meshB.n_faces()) << endl;


			//cout << "ITRA" << double(validCount)/double(nowProcessFacePairsSet.size()) << endl;
			//r_ITRA = double(validCount) / double(nowProcessFacePairsSet.size());
			r_ITRA = double(realIntersectFacesA.size() + realIntersectFacesB.size()) / double(valid_FacesA.size() + valid_FacesB.size())*100.0;
			cout << "ITRA :" << r_ITRA  << endl;

			//cout << "ITRA2:" << double(validCount) / double(nowProcessFacePairsSet.size()) << endl;
			//cout << "fakeITRA3:" << double(fakeIntersectFacesA.size() + fakeIntersectFacesB.size()) / double(valid_FacesA.size() + valid_FacesB.size()) << endl;
			//cout << "fakefaces:" << fakeIntersectFacesA.size() + fakeIntersectFacesB.size() << endl;
			//cout << "report faces:" << valid_FacesA.size() + valid_FacesB.size() << endl;
			//cout << "report real faces:" << realIntersectFacesA.size() + realIntersectFacesB.size() << endl;
			//cout << "efficiency :" << double(nowProcessFacePairsSet.size())/ double(meshA.n_faces() * meshB.n_faces()) << endl;


			if (outPutResultMesh) {

				realIntersectFacesMeshA = DataTranslate::TranslateFaceID2Mesh(meshA, realIntersectFacesA);
				realIntersectFacesMeshB = DataTranslate::TranslateFaceID2Mesh(meshB, realIntersectFacesB);

				fakeIntersectFacesMeshA = DataTranslate::TranslateFaceID2Mesh(meshA, fakeIntersectFacesA);
				fakeIntersectFacesMeshB = DataTranslate::TranslateFaceID2Mesh(meshB, fakeIntersectFacesB);


				for (Mesh::FaceIter fI = meshA.faces_begin(); fI != meshA.faces_end(); ++fI) {
					FaceId faceId = (*fI).idx();

					if (valid_FacesA.find(faceId) == valid_FacesA.end())
						inValid_FacesA.insert(faceId);
				}

				for (Mesh::FaceIter fI = meshB.faces_begin(); fI != meshB.faces_end(); ++fI) {
					FaceId faceId = (*fI).idx();

					if (valid_FacesB.find(faceId) == valid_FacesB.end())
						inValid_FacesB.insert(faceId);
				}
			
				resultInValidMeshA = DataTranslate::TranslateFaceID2Mesh(meshA, inValid_FacesA);
				resultInValidMeshB = DataTranslate::TranslateFaceID2Mesh(meshB, inValid_FacesB);
			
			}
		}

		if (outPutResultMesh) {
			string endAddressA1;
			string endAddressB1;

			string endAddressA2;
			string endAddressB2;

			string endAddressA3;
			string endAddressB3;
			switch (bspType)
			{
			case BSPConstructType::AABB_MIDDLE_SPLIT:
				endAddressA1 = "AABBResultMeshA1Red.obj";
				endAddressB1 = "AABBResultMeshB1Red.obj";

				endAddressA2 = "AABBResultMeshA2Pink.obj";
				endAddressB2 = "AABBResultMeshB2Pink.obj";

				endAddressA3 = "AABBResultMeshA3Blue.obj";
				endAddressB3 = "AABBResultMeshB3Yellow.obj";

				break;
			case BSPConstructType::SDM:
				endAddressA1 = "SDMResultMeshA1Red.obj";
				endAddressB1 = "SDMResultMeshB1Red.obj";

				endAddressA2 = "SDMResultMeshA2Pink.obj";
				endAddressB2 = "SDMResultMeshB2Pink.obj";

				endAddressA3 = "SDMResultMeshA3Blue.obj";
				endAddressB3 = "SDMResultMeshB3Yellow.obj";
				break;
			case BSPConstructType::SAH:
				endAddressA1 = "SAHResultMeshA1Red.obj";
				endAddressB1 = "SAHResultMeshB1Red.obj";

				endAddressA2 = "SAHResultMeshA2Pink.obj";
				endAddressB2 = "SAHResultMeshB2Pink.obj";

				endAddressA3 = "SAHResultMeshA3Blue.obj";
				endAddressB3 = "SAHResultMeshB3Yellow.obj";
				break;
			case BSPConstructType::ObbMiddel:
				endAddressA1 = "ObbMiddelResultMeshA1Red.obj";
				endAddressB1 = "ObbMiddelResultMeshB1Red.obj";

				endAddressA2 = "ObbMiddelResultMeshA2Pink.obj";
				endAddressB2 = "ObbMiddelResultMeshB2Pink.obj";

				endAddressA3 = "ObbMiddelResultMeshA3Blue.obj";
				endAddressB3 = "ObbMiddelResultMeshB3Yellow.obj";
				break;
			case BSPConstructType::Gravity_SPLIT:
				endAddressA1 = "Gravity_SPLITResultMeshA1Red.obj";
				endAddressB1 = "Gravity_SPLITResultMeshB1Red.obj";

				endAddressA2 = "Gravity_SPLITResultMeshA2Pink.obj";
				endAddressB2 = "Gravity_SPLITResultMeshB2Pink.obj";

				endAddressA3 = "Gravity_SPLITResultMeshA3Blue.obj";
				endAddressB3 = "Gravity_SPLITResultMeshB3Yellow.obj";
				break;
			default:
				break;
			}

			if (!OpenMesh::IO::write_mesh(realIntersectFacesMeshA, baseAddress + endAddressA1)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}

			if (!OpenMesh::IO::write_mesh(realIntersectFacesMeshB, baseAddress + endAddressB1)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}

			if (!OpenMesh::IO::write_mesh(fakeIntersectFacesMeshA, baseAddress + endAddressA2)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}

			if (!OpenMesh::IO::write_mesh(fakeIntersectFacesMeshB, baseAddress + endAddressB2)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}

			if (!OpenMesh::IO::write_mesh(resultInValidMeshA, baseAddress + endAddressA3)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}

			if (!OpenMesh::IO::write_mesh(resultInValidMeshB, baseAddress + endAddressB3)) {
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
				return make_pair(r_UTFR, r_ITRA);
			}
		}


		delete faceTree;

		return make_pair(r_UTFR, r_ITRA);
	}

	pair<double,double> TestBspNodeDistributionInfo(Mesh& meshA, Mesh& meshB, Tolerance& toler, int maxDepth, int leafShreshold, string p_BaseAddress, bool outPutResultMesh, BSPConstructType bspType)
	{
		double r_UTFR = -1;
		double r_ITRA = -1;

		string baseAddress = p_BaseAddress + to_string(maxDepth) + "//";
		cout << "maxDepth " << maxDepth << endl;
		cout << "leafShreshold " << leafShreshold << endl;

		switch (bspType)
		{
		case BSPConstructType::AABB_MIDDLE_SPLIT:
			cout << "BSPConstructType " << "AABB" << endl;
			break;
		case BSPConstructType::SDM:
			cout << "BSPConstructType " << "SDM" << endl;
			break;
		case BSPConstructType::SAH:
			cout << "BSPConstructType " << "SAH" << endl;
			break;
		case BSPConstructType::ObbMiddel:
			cout << "BSPConstructType " << "ObbMiddel" << endl;
			break;
		case BSPConstructType::Gravity_SPLIT:
			cout << "BSPConstructType " << "Gravity_SPLIT" << endl;
			break;
		case BSPConstructType::SDM_OBB:
			cout << "BSPConstructType " << "SDM_OBB" << endl;
			break;
		default:
			break;
		}


		auto tBuildTreeBegin = std::chrono::steady_clock::now();

		BSPTree* faceTree = new BSPTree(meshA, meshB, toler, maxDepth, leafShreshold, bspType);

		auto tBuildTreeEnd = std::chrono::steady_clock::now();
		cout << "建树耗时： " << std::chrono::duration<double, std::milli>(tBuildTreeEnd - tBuildTreeBegin).count() << "毫秒" << endl;


		//显示子树划分出的三角网格
		vector<BSPTreeNode*> allLeafNodes;
		faceTree->GetAllLeafNode(allLeafNodes);



		Mesh resultValidMeshA, resultValidMeshB;

		Mesh realIntersectFacesMeshA, realIntersectFacesMeshB;
		Mesh fakeIntersectFacesMeshA, fakeIntersectFacesMeshB;
		Mesh resultInValidMeshA, resultInValidMeshB;
		{
			IntersectTriangleCheckList intersectCheckList;

			// lambda for user-defined hash function
			auto hash = [](const pair<FaceId, FaceId>& c) {
				return hash_val(c.first, c.second);
			};

			// lambda for user-defined equality criterion
			auto eq = [](const pair<FaceId, FaceId>& c1, const pair<FaceId, FaceId>& c2) {
				return c1 == c2;
			};

			// create unordered set with user-defined behavior
			std::unordered_set<pair<FaceId, FaceId>, decltype(hash), decltype(eq)> nowProcessFacePairsSet(intersectCheckList.size(), hash, eq);

			auto t1 = std::chrono::steady_clock::now();


			unordered_map<int, int> depthCount;
			for (int i = 0; i <= maxDepth; i++) {
				depthCount[i] = 0;
			}

			int validNodeCount = 0;

			//for (auto nowNodes = allLeafNodes.begin(); nowNodes != allLeafNodes.end(); nowNodes++) {

			//	int& j = depthCount[(*nowNodes)->m_Depth];
			//	j++;

			//	cout << "Depth: " << (*nowNodes)->m_Depth << "nodes faces: " << (*nowNodes)->m_FacesA.size() + (*nowNodes)->m_FacesB.size() << "::: faces A size: " << (*nowNodes)->m_FacesA.size() << "::: faces B size: " << (*nowNodes)->m_FacesB.size() << endl;
			//	if ((*nowNodes)->IsValidNode()) validNodeCount++;

			//	if (((*nowNodes)->m_Depth >= 15) && ((*nowNodes)->Size() >= 2000)) {
			//		cout << "Real SizeA: " << (*nowNodes)->GetRealSize().first << "Real SizeB: " << (*nowNodes)->GetRealSize().second << endl;
			//	}

			//	(*nowNodes)->GetIntersectPair(intersectCheckList);
			//}

			//cout << "节点数为：" << allLeafNodes.size() <<" 有效节点数为：" << validNodeCount << " 占比："<<double(validNodeCount)/double(allLeafNodes.size())*100 <<"%" << endl;
			//for (int i = 0; i <= maxDepth; i++) {
			//	cout << "深度为" << i << "的节点占比：" << double(depthCount[i]) / double(allLeafNodes.size()) * 100 << "%" << endl;
			//}



			//auto t2 = std::chrono::steady_clock::now();
			//std::cout << "获取初步数据时间：" << std::chrono::duration<double, std::milli>(t2 - t1).count() << "毫秒" << endl;

			//先对采集面片对去重
			//for (auto& checkPair : intersectCheckList) {
			//	nowProcessFacePairsSet.insert(checkPair);
			//}

			auto t3 = std::chrono::steady_clock::now();
			std::cout << "查询时间：" << std::chrono::duration<double, std::milli>(t3 - t1).count() << "毫秒" << endl;


			int validCount = 0;
			unordered_set<FaceId> realIntersectFacesA, fakeIntersectFacesA;
			unordered_set<FaceId> realIntersectFacesB, fakeIntersectFacesB;

			for (auto& verifyPair : nowProcessFacePairsSet) {

				Mesh::FaceVertexCCWIter fvit = meshA.fv_ccwbegin(meshA.face_handle(verifyPair.first));
				Mesh::VertexHandle triV0 = *fvit;
				Mesh::VertexHandle triV1 = *(++fvit);
				Mesh::VertexHandle triV2 = *(++fvit);

				Triangle3d triMathA(array<Point3d, 3>{
					Point3d(meshA.point(triV0)),
						Point3d(meshA.point(triV1)),
						Point3d(meshA.point(triV2))
				});

				Mesh::FaceVertexCCWIter fvitB = meshB.fv_ccwbegin(meshB.face_handle(verifyPair.second));
				Mesh::VertexHandle triV0B = *fvitB;
				Mesh::VertexHandle triV1B = *(++fvitB);
				Mesh::VertexHandle triV2B = *(++fvitB);

				Triangle3d triMathB(array<Point3d, 3>{
					Point3d(meshB.point(triV0B)),
						Point3d(meshB.point(triV1B)),
						Point3d(meshB.point(triV2B))
				});


				if (TriangleTriangleIsIntersect(triMathA, triMathB)) {
					validCount++;
					realIntersectFacesA.insert(verifyPair.first);
					realIntersectFacesB.insert(verifyPair.second);
				}
				else {
					fakeIntersectFacesA.insert(verifyPair.first);
					fakeIntersectFacesB.insert(verifyPair.second);

				}
			}

			for (auto& fa : realIntersectFacesA) {
				fakeIntersectFacesA.erase(fa);
			}

			for (auto& fa : realIntersectFacesB) {
				fakeIntersectFacesB.erase(fa);
			}

			unordered_set<FaceId> valid_FacesA, valid_FacesB;
			unordered_set<FaceId> inValid_FacesA, inValid_FacesB;

			for (auto nowNodes = allLeafNodes.begin(); nowNodes != allLeafNodes.end(); nowNodes++) {
				(*nowNodes)->GetValidFaces(valid_FacesA, valid_FacesB);
			}

			cout << "UTFR :" << 1 - double(valid_FacesA.size() + valid_FacesB.size()) / double(meshA.n_faces() + meshB.n_faces()) << endl;
			r_UTFR = 1 - double(valid_FacesA.size() + valid_FacesB.size()) / double(meshA.n_faces() + meshB.n_faces());
			//cout << "ITRA" << double(validCount)/double(nowProcessFacePairsSet.size()) << endl;
			//r_ITRA = double(validCount) / double(nowProcessFacePairsSet.size());
			//r_ITRA = double(realIntersectFacesA.size() + realIntersectFacesB.size()) / double(valid_FacesA.size() + valid_FacesB.size());
			//auto s1 = realIntersectFacesA.size() + realIntersectFacesB.size();
			//auto s2 = valid_FacesA.size() + valid_FacesB.size();
			//cout << "ITRA :" << r_ITRA << endl;
			//cout << "1/ITRA :" << 1.0/r_ITRA << endl;

			//cout << "ITRA2:" << double(validCount) / double(nowProcessFacePairsSet.size()) << endl;
			//cout << "fakeITRA3:" << double(fakeIntersectFacesA.size() + fakeIntersectFacesB.size()) / double(valid_FacesA.size() + valid_FacesB.size()) << endl;
			//cout << "fakefaces:" << fakeIntersectFacesA.size() + fakeIntersectFacesB.size() << endl;
			//cout << "report faces:" << valid_FacesA.size() + valid_FacesB.size() << endl;
			//cout << "report real faces:" << realIntersectFacesA.size() + realIntersectFacesB.size() << endl;

			//auto t3 = std::chrono::steady_clock::now();
			//std::cout << "查询时间：" << std::chrono::duration<double, std::milli>(t3 - t1).count() << "毫秒" << endl;
		}

		delete faceTree;

		return make_pair(r_UTFR, r_ITRA);
	}



}

