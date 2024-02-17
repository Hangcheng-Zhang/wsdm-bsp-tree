#pragma once

#define PI 3.14159265358979

	class Tolerance
	{
	public:
		double m_Linear;
		double m_Squared;

	public:
		Tolerance();
		Tolerance(double p_L, double p_S);
		~Tolerance();

	public:
		double Linear() const { return m_Linear; }
		double Squared() const { return m_Squared; }

		static const Tolerance& GetTolerance() { 
			static Tolerance GLOBAL_TOLERANCE; 
			return GLOBAL_TOLERANCE;
		}

		static const Tolerance& GetToleranceHighResolution() {
			static Tolerance GLOBAL_TOLERANCE_HighResolution(1e-9,1e-9);
			return GLOBAL_TOLERANCE_HighResolution;
		}
	};

	//extern const Tolerance& GetTolerance()
	//{
	//	static Tolerance GLOBAL_TOLERANCE;
	//	return GLOBAL_TOLERANCE;
	//}


	array<double,3> GetColorRGB(int num);

	class Parameter2d
	{

	public://��Ա����
		//����������ȡֵû��Ҫ�����������¹������ã��ֱ��Ӧ���������е� 0-1 λԪ�أ�ʹ������׶�
		double u;
		double v;

	public: //��������������
		Parameter2d(); //�޲ι��캯�������𹹽�һ����Ԫ��Ϊ 0.0 �Ķ�Ԫ����
		explicit Parameter2d(double p_DefaultValue); //�Դ���Ĳ������ ��Ԫ�����Ĺ��캯��
		Parameter2d(double p_u, double p_v); //�Դ���� 2 ��ֵ���� ��Ԫ�����Ĺ��캯��

		Parameter2d(const Parameter2d& p_Parameter); //�������캯��

		~Parameter2d();

	public:
		Parameter2d& operator=(const Parameter2d& p_Parameter);

	public:
		//static Parameter2d NEGATIVE_INFINITY;
		//static Parameter2d POSITIVE_INFINITY;
	};

	class Line3d
	{
	protected:
		Vector3d m_Direction;  //!< ֱ�߷���
		Point3d m_Origin;  //!< ֱ��ԭ��

		enum class ConstructType
		{
			BY_ORIGIN_DIRECTION,
			BY_POINTA_POINTB
		};
	public:
		Line3d(const Point3d& p_Origin, const Vector3d& p_Direction, ConstructType p_CT);  //!< ����ά�����ά��������Ϊ�����Ĺ��캯��
		//Line3d(const Point3d& p_PointA, const Point3d& p_PointB);  //!< ��������ά��Ϊ�����Ĺ��캯��
		~Line3d();  //!< ��������

	public:

		const Point3d& GetOrigin() const { return m_Origin; }
		Point3d& GetOrigin() { return m_Origin; }

		Vector3d& GetDirection() { return m_Direction; }
		const Vector3d& GetDirection() const { return m_Direction; }

	public:
		virtual Point3d CalculateValue(double p_Parameter);  //!< �Դ�������߲��������������ά������
		virtual  const Point3d CalculateValue(double p_Parameter) const;  //!< �Դ�������߲��������������ά������

	};

	class LineSegment3d : public Line3d
	{
	protected:
		Point3d m_Begin;
		Point3d m_End;

	public:
		LineSegment3d(Point3d p_Begin, Point3d p_End);
		~LineSegment3d();

	public:
		Point3d& GetBegin() { return m_Begin; }
		const Point3d& GetBegin() const { return m_Begin; }

		Point3d& GetEnd() { return m_End; }
		const Point3d& GetEnd() const { return m_End; }

		double GetLength() const;

		virtual Point3d CalculateValue(double p_Parameter) override;  //!< �Դ�������߲��������������ά������
		virtual  const Point3d CalculateValue(double p_Parameter) const override;  //!< �Դ�������߲��������������ά������
	};




	class Plane3d {

	public:
		Point3d m_Origin;
		Vector3d m_Normal;

		bool m_hasCoordinate;

		Vector3d m_DirectionU;
		Vector3d m_DirectionV;
	public:
		Plane3d();
		Plane3d(Point3d p_Origin, Vector3d p_Normal);
		Plane3d(Point3d p_Origin, Vector3d p_Normal, Vector3d p_DirectionU, Vector3d p_DirectionV);
		~Plane3d();

		void MakePlane(Point3d p_Origin, Vector3d p_Normal, Vector3d p_DirectionU, Vector3d p_DirectionV);

		Point3d GetOrigin() { return m_Origin; }
		const Point3d GetOrigin() const { return m_Origin; }

		Vector3d GetNormal() { return m_Normal; }
		const Vector3d GetNormal() const { return m_Normal; }

		bool hasCoordinate() { return m_hasCoordinate; }

		const Vector3d GetDirectionU() const;
		Vector3d GetDirectionU();

		const Vector3d GetDirectionV() const;
		Vector3d GetDirectionV();
	};


	class Triangle3d
	{
	public:
		vector<Point3d> m_Vertices;
		Point3d m_Center;
		Vector3d m_Normal;
	public:
		Triangle3d(array<Point3d, 3> p_Points);
		~Triangle3d() {}


		Point3d GetCenter() { return m_Center; }
		const Point3d GetCenter() const { return m_Center; }

		Vector3d GetNormal() { return m_Normal; }
		const Vector3d GetNormal() const { return m_Normal; }

		const Point3d VertexAt(int i)  const;
		Point3d VertexAt(int i);
		
		LineSegment3d GetLinesegment(int i);
		const LineSegment3d GetLinesegment(int i) const;

		bool isValid();
	public:
		operator Plane3d() const;
	};

	Triangle3d TranslateToTriangle(MeshTriangleHandle p_Tri);


	struct TriangleBSP {
		Triangle3d m_Triangle;
		Mesh::FaceHandle m_HoldFace;
		Mesh* m_holdMesh;
	};

	class AABB {

	public:
		AABB(Point3d p1, Point3d p2);
		//AABB(vector<MeshTriangleHandle> p_Triangles);
		AABB(vector<TriangleBSP> p_Triangles);

	public:
		double xMin;
		double xMax;

		double yMin;
		double yMax;

		double zMin;
		double zMax;

	public:
		Plane3d GetLongestAxisEqualParitionPlane();
		Vector3d GetLongestAxis();

		double GetSurfaceArea();
	};


