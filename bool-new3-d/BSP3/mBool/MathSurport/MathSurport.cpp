#include "mPch.h"

#include "MathSurport.h"


	Tolerance::Tolerance() :
		m_Linear(1e-6),
		m_Squared(1e-6)
	{

	}

	Tolerance::Tolerance(double p_L, double p_S):
		m_Linear(p_L),
		m_Squared(p_S)
	{
	}

	Tolerance::~Tolerance()
	{

	}





	Parameter2d::Parameter2d() :
		u(0), v(0)
	{

	}

	Parameter2d::Parameter2d(double p_DefaultValue) :
		u(p_DefaultValue), v(p_DefaultValue)
	{

	}

	Parameter2d::Parameter2d(double p_u, double p_v) :
		u(p_u), v(p_v)
	{

	}

	Parameter2d::Parameter2d(const Parameter2d& p_Parameter) :
		u(p_Parameter.u), v(p_Parameter.v)
	{

	}

	Parameter2d::~Parameter2d()
	{

	}

	Parameter2d& Parameter2d::operator=(const Parameter2d& p_Parameter)
	{
		u = p_Parameter.u;
		v = p_Parameter.v;

		return *this;
	}


	Triangle3d::Triangle3d(array<Point3d, 3> p_Points)
	{
		m_Vertices.assign(p_Points.begin(), p_Points.end());

		auto& v = m_Vertices;
		auto& v0 = v[0];
		auto& v1 = v[1];
		auto& v2 = v[2];

		m_Center = (v0 + v1 + v2) / 3;

		Vector3d direc1 = (v0 - m_Center).normalize();
		Vector3d direc2 = (v1 - m_Center).normalize();

		m_Normal = direc1.cross(direc2);

		if (fabs(m_Normal.length() - 1) < Tolerance::GetTolerance().Linear()) {
			m_Normal.normalize();
		}

		m_Normal.normalize();

		double l1 = (m_Vertices[0] - m_Vertices[1]).length();
		double l2 = (m_Vertices[2] - m_Vertices[1]).length();
		double l3 = (m_Vertices[2] - m_Vertices[0]).length();
		assert(m_Normal.length() > Tolerance::GetTolerance().Linear());
		
	}

	const Point3d Triangle3d::VertexAt(int i) const
	{
		assert((i >= 0) && (i <= 2));

		return m_Vertices[i];
	}

	Point3d Triangle3d::VertexAt(int i)
	{
		assert((i >= 0) && (i <= 2));


		return m_Vertices[i];
	}

	LineSegment3d Triangle3d::GetLinesegment(int i)
	{
		assert(i >= 0);
		assert(i <= 2);

		if(0 == i) return LineSegment3d(m_Vertices[0], m_Vertices[1]);

		if (1 == i) return LineSegment3d(m_Vertices[1], m_Vertices[2]);

		if (2 == i) return LineSegment3d(m_Vertices[2], m_Vertices[0]);

		assert(false);

		return LineSegment3d(m_Vertices[0], m_Vertices[1]);
	}

	const LineSegment3d Triangle3d::GetLinesegment(int i) const
	{	
		assert(i >= 0);
		assert(i <= 2);

		if (0 == i) return LineSegment3d(m_Vertices[0], m_Vertices[1]);

		if (1 == i) return LineSegment3d(m_Vertices[1], m_Vertices[2]);

		if (2 == i) return LineSegment3d(m_Vertices[2], m_Vertices[0]);

		assert(false);

		return LineSegment3d(m_Vertices[0], m_Vertices[1]);
	}

	bool Triangle3d::isValid()
	{
		return !(IsZero(m_Vertices[0] - m_Vertices[1]) || IsZero(m_Vertices[1] - m_Vertices[2]) || IsZero(m_Vertices[2] - m_Vertices[0]));
	}

	Triangle3d::operator Plane3d() const
	{
		auto& v0 = m_Vertices[0];
		const Vector3d& v1 = m_Normal;

		Vector3d Ux(v0 - m_Center);
		Ux.normalize();

		Vector3d Uy = Ux.cross(v1);
		Uy.normalize();


		return Plane3d(m_Center, m_Normal, Ux, Uy);

	}

	Plane3d::Plane3d()
	{
	}

	Plane3d::Plane3d(Point3d p_Origin, Vector3d p_Normal):
		m_Origin(p_Origin), m_Normal(p_Normal), m_hasCoordinate(false)
	{
		assert(m_Normal.length() > Tolerance::GetTolerance().Linear());
		m_Normal.normalize();
	}

	Plane3d::Plane3d(Point3d p_Origin, Vector3d p_Normal, Vector3d p_DirectionU, Vector3d p_DirectionV) :
		m_Origin(p_Origin), m_Normal(p_Normal), m_DirectionU(p_DirectionU), m_DirectionV(p_DirectionV), m_hasCoordinate(true)
	{
		assert(m_Normal.length() > Tolerance::GetTolerance().Linear());
		m_Normal.normalize();
	}

	Plane3d::~Plane3d()
	{

	}

	void Plane3d::MakePlane(Point3d p_Origin, Vector3d p_Normal, Vector3d p_DirectionU, Vector3d p_DirectionV)
	{
		m_Origin = p_Origin;

		m_Normal = p_Normal;
		assert(m_Normal.length() > Tolerance::GetTolerance().Linear());
		m_Normal.normalize();

		m_DirectionU = p_DirectionU;
		m_DirectionV = p_DirectionV;
		m_hasCoordinate = true;
	}

	const Vector3d Plane3d::GetDirectionU() const
	{
		if (!m_hasCoordinate) assert(false);
		return m_DirectionU;
	}

	Vector3d Plane3d::GetDirectionU()
	{
		if (!m_hasCoordinate) assert(false);
		return m_DirectionU;
	}

	const Vector3d Plane3d::GetDirectionV() const
	{
		if (!m_hasCoordinate) assert(false);
		return m_DirectionV;
	}

	Vector3d Plane3d::GetDirectionV()
	{
		if (!m_hasCoordinate) assert(false);
		return m_DirectionV;
	}

	Line3d::Line3d(const Point3d& p_Origin, const Vector3d& p_Direction, ConstructType p_CT)
	{
		switch (p_CT) {
		case ConstructType::BY_ORIGIN_DIRECTION:
		{
			m_Origin = p_Origin;

			assert(p_Direction.length() > Tolerance::GetTolerance().Linear());
			p_Direction.normalized();
			m_Direction = p_Direction;
			break;
		}
		case ConstructType::BY_POINTA_POINTB:
		{
			m_Origin = p_Origin;
			m_Direction = (p_Direction - p_Origin).normalize();
			break;
		}
		default:

			break;
		}

	}

	Line3d::~Line3d()
	{

	}

	Point3d Line3d::CalculateValue(double p_Parameter)
	{
		Point3d r_Point = m_Origin + m_Direction * p_Parameter;

		return r_Point;
	}

	const Point3d Line3d::CalculateValue(double p_Parameter) const
	{
		Point3d r_Point = m_Origin + m_Direction * p_Parameter;

		return r_Point;
	}

	LineSegment3d::LineSegment3d(Point3d p_Begin, Point3d p_End) :
		Line3d(p_Begin, p_End, Line3d::ConstructType::BY_POINTA_POINTB), m_Begin(p_Begin), m_End(p_End)
	{

	}

	LineSegment3d::~LineSegment3d()
	{

	}

	double LineSegment3d::GetLength() const
	{

		return (m_Begin - m_End).length();
	}

	Point3d LineSegment3d::CalculateValue(double p_Parameter)
	{
		return Line3d::CalculateValue(p_Parameter);
		//return m_Begin * p_Parameter + (1 - p_Parameter) * m_End;
	}

	const Point3d LineSegment3d::CalculateValue(double p_Parameter) const
	{

		return Line3d::CalculateValue(p_Parameter);

		//return m_Begin * p_Parameter + (1 - p_Parameter) * m_End;
	}

	array<double, 3> GetColorRGB(int num)
	{
		num = abs(num) % 1000;

		assert(num >= 0);
		assert(num <= 999);

		double r = 999 - num;
		double g = abs(num - 500);
		double b = num;

		r /= 999;
		g /= 500;
		b /= 999;

		array<double, 3> r_RGB{r,g,b};

		return r_RGB;
	}



	Triangle3d TranslateToTriangle(MeshTriangleHandle p_Triangle)
	{
		Mesh::FaceVertexCCWIter fvit = p_Triangle.second->fv_ccwbegin(p_Triangle.first);
		Mesh::VertexHandle triV0 = *fvit;
		Mesh::VertexHandle triV1 = *(++fvit);
		Mesh::VertexHandle triV2 = *(++fvit);

		assert(triV0.is_valid());
		assert(triV1.is_valid());
		assert(triV2.is_valid());

		Point3d p1(p_Triangle.second->point(triV0));
		Point3d p2(p_Triangle.second->point(triV1));
		Point3d p3(p_Triangle.second->point(triV1));

		return Triangle3d(array<Point3d, 3>{
				Point3d(p_Triangle.second->point(triV0)),
				Point3d(p_Triangle.second->point(triV1)),
				Point3d(p_Triangle.second->point(triV2))
		});
	}

	AABB::AABB(Point3d p1, Point3d p2)
	{
		//get xmin xmax
		if (p1[0] <= p2[0]) {	
			xMin = p1[0];
			xMax = p2[0];
		}
		else {
			xMin = p2[0];
			xMax = p1[0];		
		}

		//get ymin ymax
		if (p1[1] <= p2[1]) {
			yMin = p1[1];
			yMax = p2[1];
		}
		else {
			yMin = p2[1];
			yMax = p1[1];
		}

		//get zmin zmax
		if (p1[2] <= p2[2]) {
			zMin = p1[2];
			zMax = p2[2];
		}
		else {
			zMin = p2[2];
			zMax = p1[2];
		}
	}

	AABB::AABB(vector<TriangleBSP> p_Triangles)
	{
		assert(!p_Triangles.empty());

		//初始化aabb6个顶点值
		Triangle3d tInitial = p_Triangles.front().m_Triangle;

		Point3d pInitial = tInitial.VertexAt(0);

		double txMin = pInitial[0];
		double txMax = pInitial[0];

		double tyMin = pInitial[1];
		double tyMax = pInitial[1];

		double tzMin = pInitial[2];
		double tzMax = pInitial[2];

		for (auto& face : p_Triangles) {
			Triangle3d tri = face.m_Triangle;

			for (int i = 0; i < 3;i++) {
				Point3d v = tri.VertexAt(i);
				if (v[0] < txMin) txMin = v[0];
				if (v[0] > txMax) txMax = v[0];

				if (v[1] < tyMin) tyMin = v[1];
				if (v[1] > tyMax) tyMax = v[1];

				if (v[2] < tzMin) tzMin = v[2];
				if (v[2] > tzMax) tzMax = v[2];
			}

		}

		xMin = txMin;
		xMax = txMax;

		yMin = tyMin;
		yMax = tyMax;

		zMin = tzMin;
		zMax = tzMax;
	}



	Plane3d AABB::GetLongestAxisEqualParitionPlane()
	{
		Point3d center((xMax + xMin) / 2, (yMax + yMin) / 2, (zMax + zMin) / 2);

		double direcXLength = xMax - xMin;
		double direcYLength = yMax - yMin;
		double direcZLength = zMax - zMin;

		if ((direcXLength >= direcYLength) && (direcXLength >= direcZLength)) {
			//x方向最长

			return Plane3d(center, Vector3d(1, 0, 0));
		}
		else if ((direcYLength >= direcXLength) && (direcYLength >= direcZLength)) {
			//y方向最长

			return Plane3d(center, Vector3d(0, 1, 0));
		}
		else if ((direcZLength >= direcXLength) && (direcZLength >= direcYLength)) {
			//z方向最长

			return Plane3d(center, Vector3d(0, 0, 1));
		}


		assert(false);
		return Plane3d();
	}

	Vector3d AABB::GetLongestAxis()
	{

		double direcXLength = xMax - xMin;
		double direcYLength = yMax - yMin;
		double direcZLength = zMax - zMin;

		if ((direcXLength >= direcYLength) && (direcXLength >= direcZLength)) {
			//x方向最长

			return Vector3d(1, 0, 0);
		}
		else if ((direcYLength >= direcXLength) && (direcYLength >= direcZLength)) {
			//y方向最长

			return Vector3d(0, 1, 0);
		}
		else if ((direcZLength >= direcXLength) && (direcZLength >= direcYLength)) {
			//z方向最长

			return Vector3d(0, 0, 1);
		}


		assert(false);
		return Vector3d();

	}

	double AABB::GetSurfaceArea()
	{		
		return 2.0 * ((xMax - xMin) * (yMax-yMin) + (xMax - xMin) * (zMax - zMin) + (yMax - yMin) * (zMax - zMin));
	}
