#include "RayTracing.h"
#include <omp.h>
#include "../../Util/GEMS_Memory.h"

calculation::RayTracing::RayTracing(STLMirror* _mirror)
	:mirror(_mirror)
{
}

calculation::RayTracing::~RayTracing()
{
}

void calculation::RayTracing::setMirror(STLMirror * mirror)
{
	this->mirror = mirror;
}

void calculation::RayTracing::calcNormalOfLine_Mirror(const Vector3 & startPiont,
	const Vector3 & direction, Vector3 & normal, Vector3 & intersection,
	bool & isIntersect, double & t)
{
	//调用omp版本
	calcNormalOfLine_MirrorByPolyData_omp(startPiont, direction, normal,
		intersection, isIntersect, t);	
}


void calculation::RayTracing::calcReflect(const Vector3 & startPiont, const Vector3 & direction, 
	Vector3 & reflex, Vector3 & intersection, bool & isIntersect)
{
	calcReflectByPolyData(startPiont, direction, reflex, intersection, isIntersect);
}


void calculation::RayTracing::calcReflectByPolyDataBatch(const vector<vector<Vector3>>& startPiont,
	const vector<vector<Vector3>>& direction,
	vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
	vector<vector<bool>> &isIntersect)
{
	for (int i = 0; i < startPiont.size(); i++)
		for (int j = 0; j < startPiont[i].size(); j++)
		{
			bool isTmep = false;
			calcReflectByPolyData(startPiont[i][j], direction[i][j], reflex[i][j],
				intersection[i][j], isTmep);
			isIntersect[i][j] = isTmep;
		}
}

void calculation::RayTracing::calcReflectByPolyData(const Vector3 & startPiont, 
	const Vector3 & direction, Vector3 & reflex, 
	Vector3 & intersection, bool & isIntersect)
{
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();

	isIntersect = false;//没找就还没找到
	/*
	//Jin
	//threadNum
	int threadNum = 8;
	vector<int> Ns_omp(threadNum);
	vector<int> Nn_omp(threadNum);
	//omp Dividing
	for (int i = 0; i<threadNum; i++) {
		Nn_omp[i] = (EleNum) / threadNum;
		int rr = (EleNum) % threadNum;
		if (i<rr && rr != 0)  Nn_omp[i] += 1;
		if (i == 0) {
			Ns_omp[i] = 0;
		}
		else {
			Ns_omp[i] = Ns_omp[i - 1] + Nn_omp[i - 1];
		}
	}
	//omp version
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double t;
		int id = omp_get_thread_num();
		Vector3 startPiont_omp = startPiont;
		Vector3 direction_omp = direction;
		Vector3 intersection_omp;
		for (int i = Ns_omp[id]; i < Ns_omp[id] + Nn_omp[id]; i++)  //求与反射面的交点 OMP分段
		{
			if (isIntersect) continue;	//如果有找到了的就空转

			vtkIdList * p;
			p = polyData->GetCell(i)->GetPointIds();
			double * point;
			point = polyData->GetPoint(p->GetId(0));
			Vector3 NodesXYZ1(point[0], point[1], point[2]);
			point = polyData->GetPoint(p->GetId(1));
			Vector3 NodesXYZ2(point[0], point[1], point[2]);
			point = polyData->GetPoint(p->GetId(2));
			Vector3 NodesXYZ3(point[0], point[1], point[2]);

			if (this->isIntersect(startPiont_omp, direction_omp, NodesXYZ1,
				NodesXYZ2, NodesXYZ3, intersection_omp, t))
			{
				if (t >= 0)
				{
					Vector3 tempa = NodesXYZ1 - NodesXYZ2;
					Vector3 tempb = NodesXYZ1 - NodesXYZ3;
					Vector3 n_light = tempa.Cross(tempb);  //法向量

					isIntersect = true;	//Jin 找到了
					intersection = intersection_omp;
					reflex = reflectLight(direction_omp, n_light);
					//return;//不能return了
				}
			}
		}//for
	}//omp
	//isIntersect = false;
	*/
	
	double t;
	for (int i = 0; i < EleNum; i++)  //求与反射面的交点
	{
		vtkIdList * p;
		p = polyData->GetCell(i)->GetPointIds();
		double * point;
		point = polyData->GetPoint(p->GetId(0));
		Vector3 NodesXYZ1(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(1));
		Vector3 NodesXYZ2(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(2));
		Vector3 NodesXYZ3(point[0], point[1], point[2]);
		if (this->isIntersect(startPiont, direction, NodesXYZ1,
			NodesXYZ2, NodesXYZ3, intersection, t))
		{
			if (t >= 0)
			{
				Vector3 tempa = NodesXYZ1 - NodesXYZ2;
				Vector3 tempb = NodesXYZ1 - NodesXYZ3;
				Vector3 n_light = tempa.Cross(tempb);  //法向量

				isIntersect = true;
				reflex = reflectLight(direction, n_light);
				return;
			}
		}
	}
	isIntersect = false;
	
}
void calculation::RayTracing::SetupNodes(void)
{	//添加这个函数是为了写一个buffer,OMP中读PloyData老出错 _Jin
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();
	NodesXYZ1.resize(EleNum);
	NodesXYZ2.resize(EleNum);
	NodesXYZ3.resize(EleNum);
	for (int i = 0; i < EleNum; i++)  //求与反射面的交点
	{
		vtkIdList * p;
		p = polyData->GetCell(i)->GetPointIds();
		double * point;
		point = polyData->GetPoint(p->GetId(0));
		NodesXYZ1[i].set(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(1));
		NodesXYZ2[i].set(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(2));
		NodesXYZ3[i].set(point[0], point[1], point[2]);
	}
}


void calculation::RayTracing::calcNormalOfLine_MirrorByPolyData_omp(
	const Vector3 & startPiont,const Vector3 & direction, Vector3 & normal, 
	Vector3 & intersection, bool & isIntersect, double & t)
{
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();
	
	
	//Jin _omp
	isIntersect = false;//没找就还没找到	

	int threadNum = 8;	//默认8线程，依据机器还可以更多 超线程还是有些作用的
	vector<int> Ns_omp(threadNum);
	vector<int> Nn_omp(threadNum);
	//omp Dividing 分配网格点数
	for (int i = 0; i<threadNum; i++) {
		Nn_omp[i] = (EleNum) / threadNum;
		int rr = (EleNum) % threadNum;
		if (i<rr && rr != 0)  Nn_omp[i] += 1;
		if (i == 0) {
			Ns_omp[i] = 0;
		}
		else {
			Ns_omp[i] = Ns_omp[i - 1] + Nn_omp[i - 1];
		}
	}
	//omp version
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
	
		int id = omp_get_thread_num();
		//int id = omp_get_thread_num
		Vector3 startPiont_omp = startPiont;
		Vector3 direction_omp = direction;
		Vector3 intersection_omp = intersection;
		double t_omp = t;
		vtkIdList * p; 
		//double* point1=NULL;	point1 = Allocate_1D(point1, 3); Vector3 NodesXYZ1;
		//double* point2=NULL;	point2 = Allocate_1D(point2, 3); Vector3 NodesXYZ2;
		//double* point3=NULL;	point3 = Allocate_1D(point3, 3); Vector3 NodesXYZ3;
		for (int i = Ns_omp[id]; i < Ns_omp[id] + Nn_omp[id]; i++)  //求与反射面的交点 OMP分段
		{
			if (isIntersect) break;	//如果有其他线程找跳出当前循环

			//Omp里读polyData老出错，就写了个Buffer

			//p = polyData->GetCell(i)->GetPointIds();
			//point1 = polyData->GetPoint(p->GetId(0));
			//NodesXYZ1.set(point1[0], point1[1], point1[2]);
			//point2 = polyData->GetPoint(p->GetId(1));
			//NodesXYZ2.set(point2[0], point2[1], point2[2]);
			//point3 = polyData->GetPoint(p->GetId(2));
			//NodesXYZ3.set(point3[0], point3[1], point3[2]);

			//判断是否相交
			if (this->isIntersect(startPiont_omp, direction_omp, NodesXYZ1[i],
				NodesXYZ2[i], NodesXYZ3[i], intersection_omp, t_omp))
			{
				Vector3 tempa = NodesXYZ1[i] - NodesXYZ2[i];
				Vector3 tempb = NodesXYZ1[i] - NodesXYZ3[i];
				normal = tempa.Cross(tempb);  //法向量
				t = t_omp;					  //判断是镜前还是镜后
				intersection = intersection_omp;	//交点位置，用于后续计算传播距离

				isIntersect = true;	//Jin 找到了
				//return;//不能return了
				
				break;//结束整个for
			}
		}//for each OMP id

	}//omp
	//Jin//
	
}

void calculation::RayTracing::calcNormalOfLine_MirrorByPolyData(
	const Vector3 & startPiont, const Vector3 & direction, Vector3 & normal,
	Vector3 & intersection, bool & isIntersect, double & t)
{
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();

	 //原来单线程版本的
	 isIntersect = false;
	 for (int i = 0; i < EleNum; i++)  //求与反射面的交点
	 {
	 vtkIdList * p;
	 p = polyData->GetCell(i)->GetPointIds();
	 double * point;
	 point = polyData->GetPoint(p->GetId(0));
	 Vector3 NodesXYZ1(point[0], point[1], point[2]);
	 point = polyData->GetPoint(p->GetId(1));
	 Vector3 NodesXYZ2(point[0], point[1], point[2]);
	 point = polyData->GetPoint(p->GetId(2));
	 Vector3 NodesXYZ3(point[0], point[1], point[2]);
	 if (this->isIntersect(startPiont, direction, NodesXYZ1,
	 NodesXYZ2, NodesXYZ3, intersection, t))
	 {
	 Vector3 tempa = NodesXYZ1 - NodesXYZ2;
	 Vector3 tempb = NodesXYZ1 - NodesXYZ3;
	 normal = tempa.Cross(tempb);  //法向量

	 isIntersect = true;
	 return;
	 }
	 }
	 

}

bool calculation::RayTracing::isIntersect(const Vector3 & orig, const Vector3 & dir,
	const Vector3 & v0, const Vector3 & v1, const Vector3 & v2, 
	Vector3 & intersection, double & t)
{
	double u, v;
	// E1
	Vector3 E1 = v1 - v0;

	// E2
	Vector3 E2 = v2 - v0;

	// P
	Vector3 P = dir.Cross(E2);

	// determinant
	double det = E1.Dot(P);

	Vector3 T;
	T = orig - v0;

	// If determinant is near zero, ray lies in plane of triangle
	//if (det < 0.00000001 && det > -0.00000001)
	//	return false;

	// Calculate u and make sure u <= 1
	u = T.Dot(P);
	double fInvDet = 1.0f / det;
	u *= fInvDet;

	if (u < 0.0 || u > 1)
		return false;

	// Q
	Vector3 Q = T.Cross(E1);

	// Calculate v and make sure u + v <= 1
	v = dir.Dot(Q);
	v *= fInvDet;
	if (v < 0.0 || u + v > 1)
		return false;

	// Calculate t, scale parameters, ray intersects triangle
	t = E2.Dot(Q);
	t *= fInvDet;

	intersection = orig + dir * t;

	return true;
}


Vector3 calculation::RayTracing::reflectLight(const Vector3 & a, const Vector3 & n)
{
	//先单位化
	double absa = pow(a.Dot(a), 0.5);
	double absn = pow(n.Dot(n), 0.5);
	Vector3 tempa = a * (1 / absa);
	Vector3 tempn = n * (1 / absn);
	double I = 2 * tempn.Dot(tempa);
	if (I < 0)
		I = -I;
	else
		tempa = Vector3(0.0, 0.0, 0.0) - tempa;

	return tempn * I + tempa;
}

