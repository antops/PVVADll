#include "pvva.h"
#include "../../Util/Position3D.h"
#include "../../Util/Matrix4D.h"
#include "../../Util/Vector3D.h"
#include "../../Util/Rays.h"
//#include "../../CUDADLL/CUDARayTracing.h"
#include "../../FFTPropagatorDLL/FFTPropagator.h"
#include "../../FFTASRotationDLL/FFTASRotation.h"
#include <cmath>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkSTLWriter.h>
#include <omp.h>


using namespace calculation;


PVVA::PVVA()
{
	Initialization();
	CATRMode = false;
}

PVVA::~PVVA()
{
	Free_2D(Px);
	Free_2D(Py);
	Free_2D(Pz);

	Free_2D(Plane);
	Free_2D(n_Plane);
	Free_2D(R_Plane);
	
	Free_2D(Ex1);
	Free_2D(Ey1);
	Free_2D(Ez1);
	Free_2D(Hx1);
	Free_2D(Hy1);
	Free_2D(Hz1);

	Free_2D(Ex_R);
	Free_2D(Ey_R);
	Free_2D(Ez_R);
}

void PVVA::Initialization()
{
	N = 0; M = 0;
	f = 10.65e9;
	lamda = 0; 
	//k = 0;
	dx = 0;	dy = 0;	z0 = 0;
	unit_factor = 1;

	Ex_In = NULL;
	Ey_In = NULL;

	Px = NULL;
	Py = NULL;
	Pz = NULL;

	Plane = NULL;
	n_Plane = NULL;
	R_Plane = NULL;
	Rn_Plane = NULL;

	Ex1 = NULL;
	Ey1 = NULL;
	Ez1 = NULL;
	Hx1 = NULL;;
	Hy1 = NULL;
	Hz1 = NULL;

	Ex_R = NULL;
	Ey_R = NULL;
	Ez_R = NULL;

	returnFloat = NULL;
	user = NULL;

	threadNum = 4;
}

void PVVA::setUnit(double factor)
{
	unit_factor = factor;
}

double PVVA::IntersectPlane(const Vector3 & orig, const Vector3 & dir,
	const Vector3 & Plane_org, const Vector3 & Plane_n, Vector3 & intersection)
{
	double temp1 = dir.Dot(Plane_n);
	if(temp1 == 0)
		return 0.0;
	double temp2 = Plane_org.Dot(Plane_n) - orig.Dot(Plane_n);
	double t = temp2 / temp1;
	intersection = orig + dir * t;
	if (t > 10)
	{
		double l = 0;
	}
	return t;

}

double PVVA::IntersectPoint(const Vector3 &orig, const Vector3 &dir,
	const Vector3 &v0, const Vector3 &E1, const Vector3 &E2, Vector3 &intersection)
{
	// P
	Vector3 P = dir.Cross(E2);

	// determinant
	double det = E1.Dot(P);

	double u, v, t;
	Vector3 T;
	T = orig - v0;

	// If determinant is near zero, ray lies in plane of triangle
	if (det < 10e-10 && det > -10e-10)
		return 0;

	u = T.Dot(P);
	double fInvDet = 1.0f / det;
	u *= fInvDet;

	// Q
	Vector3 Q = T.Cross(E1);

	v = dir.Dot(Q);
	v *= fInvDet;

	// Calculate t, scale parameters, ray intersects triangle
	t = E2.Dot(Q);
	t *= fInvDet;

	intersection = orig + dir * t;
	return t;
}

void PVVA::Poynting()
{
	//计算全局坐标系中，镜前虚拟面场的射线描述
	double absx, absy, absz;
	//double sum = 0;

	// 源的坐标系转换到绝对坐标系（只求旋转）
	Vector3D RotateAsixSou(SourceGraphTrans.getRotate_x(), 
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrixSou = Matrix4D::getRotateMatrix(
		SourceGraphTrans.getRotate_theta(), RotateAsixSou);


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Px[i][j] = ConjugateMul(Ey1[i][j], Hz1[i][j]) - ConjugateMul(Ez1[i][j], Hy1[i][j]);
			Py[i][j] = ConjugateMul(Ez1[i][j], Hx1[i][j]) - ConjugateMul(Ex1[i][j], Hz1[i][j]);
			Pz[i][j] = ConjugateMul(Ex1[i][j], Hy1[i][j]) - ConjugateMul(Ey1[i][j], Hx1[i][j]);

			absx = Px[i][j].real();
			absy = Py[i][j].real();
			absz = Pz[i][j].real();

			n_Plane[i][j].set(absx, absy, absz);

			n_Plane[i][j] = rotatMatrixSou * n_Plane[i][j];
			n_Plane[i][j].Normalization();
		}
	}
}

void PVVA::getPlane(Vector3 **& _org, Vector3 **& _n) const
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			_org[i][j] = Plane[i][j];
			_n[i][j] = n_Plane[i][j];
		}
}

void PVVA::SetReturnFloat(void(*returnFloat)(float, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}

void PVVA::updateSource_n(const Vector3& new_n)
{
	//new_n.Normalization();
	if (new_n.x != 0 || new_n.y != 0 || new_n.z != 1)
	{
		Vector3 rotate_axis = Vector3(0, 0, 1).Cross(new_n); // 旋转轴
		double rotate_theta = acos(Vector3(0, 0, 1).Dot(new_n));
		rotate_theta = rotate_theta / Pi * 180;
		SourceGraphTrans.updateRotate(rotate_axis, rotate_theta);
	}
	else
	{
		SourceGraphTrans.updateRotate(Vector3(0, 0, 1), 0);
	}
}

void PVVA::Result()
//这里
{
	CalArbitraryPlane(SourceGraphTrans, TargetGraphTrans);
	
}

void PVVA::CalAmplitude()
{
	Vector3 tempA1, tempA2;
	Vector3 tempB1, tempB2;
	Vector3 tempC1, tempC2;
	Vector3 tempD1, tempD2;
	double tempcos1, tempcos2;
	double PreSquare, PostSquare;
	double tempratio;
	double LAE, TAE;
	n_Source.Normalization();
	//中心部分四个顶点
	for (int i = 1; i < N - 2; i++)
	{
		for (int j = 1; j < M - 2; j++)
		{
			// 反射前投影面积
			tempA1 = Plane[i - 1][j];
			tempB1 = Plane[i][j + 1];
			tempC1 = Plane[i + 1][j];
			tempD1 = Plane[i][j - 1];

			tempcos1 = n_Plane[i][j].Dot(n_Source);
		
			PreSquare = tempcos1 * dx * dy * 2;
			LAE = pow((Ex1[i][j].real() * Ex1[i][j].real() + Ex1[i][j].imag() * Ex1[i][j].imag()
				+ Ey1[i][j].real() * Ey1[i][j].real() + Ey1[i][j].imag() * Ey1[i][j].imag()
				+ Ez1[i][j].real() * Ez1[i][j].real() + Ez1[i][j].imag() * Ez1[i][j].imag()), 0.5);
			// 反射后投影面积
			tempA2 = R_Plane[i - 1][j];
			tempB2 = R_Plane[i][j + 1];
			tempC2 = R_Plane[i + 1][j];
			tempD2 = R_Plane[i][j - 1];
			tempcos2 = Rn_Plane[i][j].Dot(R_Source);

			TAE = pow((Ex_R[i][j].real() * Ex_R[i][j].real() + Ex_R[i][j].imag() * Ex_R[i][j].imag()
				+ Ey_R[i][j].real() * Ey_R[i][j].real() + Ey_R[i][j].imag() * Ey_R[i][j].imag()
				+ Ez_R[i][j].real() * Ez_R[i][j].real() + Ez_R[i][j].imag() * Ez_R[i][j].imag()), 0.5);

			//PostSquare = tempcos2 * 
			double SquareTest1 = CalSquare(tempA2, tempB2, tempC2, tempD2);
			double SquareTest2 = CalSquare(tempB2, tempC2, tempD2, tempA2);
			if (abs(SquareTest1 - SquareTest2) < 0.00000001)
				PostSquare = tempcos2 * SquareTest1;
			else
			{
				PostSquare = PreSquare;
			}
			if (PostSquare == 0)
				PostSquare = 1;
			if (TAE == 0)
				TAE = 1;
			tempratio = pow(fabs(PreSquare / PostSquare), 0.5) * LAE / TAE;

			Ex_R[i][j] = Ex_R[i][j] * tempratio;
			Ey_R[i][j] = Ey_R[i][j] * tempratio;
			Ez_R[i][j] = Ez_R[i][j] * tempratio;
		}
	}
}

void PVVA::InterVal()
{

	FILE* writefile;
	string fileAddress = "./PVVD/out0.dat";
	errno_t err;
	err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	double tx, ty, tz, rx, ry, rz, rth;
	int flag = 0;
	fwrite(&f, sizeof(double), 1, writefile);
	tx = SourceGraphTrans.getTrans_x(); ty = SourceGraphTrans.getTrans_y(); tz = SourceGraphTrans.getTrans_z();
	rx = SourceGraphTrans.getRotate_x(); ry = SourceGraphTrans.getRotate_y(); rz = SourceGraphTrans.getRotate_z();
	rth = SourceGraphTrans.getRotate_theta();
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex_R[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey_R[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);

	// 逆矩阵
	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(),
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D R_rotatMatrix = Matrix4D::getRotateMatrix(
		-SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D R_translateMatrix = Matrix4D::getTranslateMatrix(
		-SourceGraphTrans.getTrans_x(),
		-SourceGraphTrans.getTrans_y(), -SourceGraphTrans.getTrans_z());
	Matrix4D R_Matrix = R_rotatMatrix * R_translateMatrix;
		
	bool ** hit = NULL;//判断相没相交
	hit = Allocate_2D(hit, N, M);

	Vector3 tempRPlane;
	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			if (R_Plane[i][j].x == 100)
				hit[i][j] = false;
			else
				hit[i][j] = true;
			tempRPlane = R_Matrix * R_Plane[i][j];
			R_Plane[i][j] = tempRPlane;
		}
	

		for(int i = 0; i<N;i++)
			for (int j = 0; j < M; j++) {
				Ex_In[i][j] = complex<double>(0, 0);
				Ey_In[i][j] = complex<double>(0, 0);
				Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
				Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
			}
		int NumSlice = N / 100;
		if (NumSlice == 0) NumSlice = 1;
		vector<int> NStart;
		vector<int> NNum;
		NStart.resize(NumSlice);		NNum.resize(NumSlice);

		for (int ii = 0; ii < NumSlice; ii++) {
			int nn = N / NumSlice;
			int resn = N%NumSlice;
			if ((ii < resn) && (resn != 0)) NNum[ii] = nn + 1;
			else NNum[ii] = nn;
		}
		for (int ii = 0; ii < NumSlice; ii++) {
			if (ii == 0) NStart[ii] = 0;
			else	NStart[ii] = NStart[ii - 1] + NNum[ii - 1];
		}
		//omp acceleration Preparation
		int totalNum = NumSlice*NumSlice;
		vector<int> startS;	startS.resize(threadNum);
		vector<int> numS;	numS.resize(threadNum);
		vector<int> stopS;	stopS.resize(threadNum);
		for (int i = 0; i < threadNum; i++) {
			int num = totalNum / threadNum;
			int rr = (totalNum) % threadNum;
			if (i<rr && rr != 0) num += 1;
			numS[i] = num;
			//starti
			if (i == 0) startS[i] = 0;
			else startS[i] = startS[i - 1] + numS[i - 1];
			//stopi
			stopS[i] = startS[i] + numS[i];
		}

		omp_set_num_threads(threadNum);
		#pragma omp parallel
		{	
			int id = omp_get_thread_num();
			for (int ip = 0; ip < numS[id]; ip++) {
				int io = ip*threadNum + id;

				int is = io / NumSlice;
				int js = io % NumSlice;

				if (id == 0) { if (returnFloat) returnFloat(((ip-startS[id]+0.5)/numS[id] * 20 + 60), user); }
					vtkFloatArray *scalars = vtkFloatArray::New();
					vtkSmartPointer<vtkPoints> points =
						vtkSmartPointer<vtkPoints>::New();

					double MAXX = dx * (N - 1) / 2;
					double MAXY = dy * (M - 1) / 2;


					int NxStart = NStart[is];	if (is > 0) NxStart = NxStart - 4;
					int NxEnd = NStart[is] + NNum[is];	if (is < NumSlice - 1) NxEnd = NxEnd + 4;
					//int NxNum = NxEnd - NxStart;

					int NyStart = NStart[js];	if (js > 0) NyStart = NyStart - 4;
					int NyEnd = NStart[js] + NNum[js];	if (js < NumSlice - 1) NyEnd = NyEnd + 4;
					//int NyNum = NyEnd - NyStart;

					int num = 0;
					//绘制点集-反射光线位置 将反射光线落点转换到投影面坐标系上
					for (int i = NxStart; i < NxEnd; i++) {
						for (int j = NyStart; j < NyEnd; j++) {
							if (!hit[i][j])
								continue;
							if (R_Plane[i][j].x < -MAXX || R_Plane[i][j].x > MAXX ||
								R_Plane[i][j].y < -MAXY || R_Plane[i][j].y > MAXY)
								continue;
							points->InsertNextPoint(R_Plane[i][j].x, R_Plane[i][j].y, 0);
							scalars->InsertTuple1(num, i*M + j + 0.01);//这里巧用Scalar其实是为了标记点的index
							num++;
						}
					}

					vtkSmartPointer<vtkPolyData> polydata =
						vtkSmartPointer<vtkPolyData>::New();
					polydata->SetPoints(points);
					polydata->GetPointData()->SetScalars(scalars);

					double scalarsRange[2];
					scalars->GetRange(scalarsRange);
					scalars->Delete();

					vtkSmartPointer<vtkDelaunay2D> delaunay =
						vtkSmartPointer<vtkDelaunay2D>::New();
					delaunay->SetInputData(polydata);
					delaunay->Update();
					polydata = delaunay->GetOutput();


					int EleNum = polydata->GetNumberOfCells();

					if (CATRMode) {//紧缩场模式，//去掉长边cell
						double ds = (dx + dy) / 0.5;

						float Maxlength;
						vector<int> LongTriangleId;	LongTriangleId.resize(0);
						vector<Vector2> Trip; Trip.resize(3);
						vector<double> ll;	ll.resize(3);
						for (int i = 0; i < EleNum; i++) {
							vtkIdList * p;
							p = polydata->GetCell(i)->GetPointIds();
							double * point;
							point = polydata->GetPoint(p->GetId(0));	Trip[0] = Vector2(point[0], point[1]);
							point = polydata->GetPoint(p->GetId(1));	Trip[1] = Vector2(point[0], point[1]);
							point = polydata->GetPoint(p->GetId(2));	Trip[2] = Vector2(point[0], point[1]);
							ll[0] = Vector2(Trip[0].x - Trip[1].x, Trip[0].y - Trip[1].y).Length();
							ll[1] = Vector2(Trip[1].x - Trip[2].x, Trip[1].y - Trip[2].y).Length();
							ll[2] = Vector2(Trip[2].x - Trip[0].x, Trip[2].y - Trip[0].y).Length();
							Maxlength = 0;
							for (int k = 0; k < 3; k++) {
								if (ll[k] > Maxlength) Maxlength = ll[k];
							}
							if (Maxlength > 2 * ds) {
								LongTriangleId.insert(LongTriangleId.end(), i);
							}
						}
						int indexNum = LongTriangleId.size();
						for (int i = 0; i < indexNum; i++) {
							polydata->DeleteCell(LongTriangleId[indexNum - 1 - i]);
						}
						polydata->RemoveDeletedCells();
						polydata->Modified();
					}

					EleNum = polydata->GetNumberOfCells();

					for (int i = 0; i < EleNum; i++)  //求与反射面的交点
					{
						double maxX = -MAXX;
						double minX = MAXX;
						double maxY = -MAXY;
						double minY = MAXY;

						vtkIdList * p;
						p = polydata->GetCell(i)->GetPointIds();
						double * point;

						vector<Vector2> A(3);
						vector<int> dataIJ(3);
						Vector2 P;
						for (int k = 0; k < 3; k++)
						{
							dataIJ[k] = polydata->GetPointData()
								->GetScalars()->GetVariantValue(p->GetId(k)).ToInt();//取出Scalar
							point = polydata->GetPoint(p->GetId(k));//找到点的位置
							if (maxX < point[0])
								maxX = point[0];
							if (minX > point[0])
								minX = point[0];
							if (maxY < point[1])
								maxY = point[1];
							if (minY > point[1])
								minY = point[1];
							A[k].set(point[0], point[1]);
						}
						int iiMin = ceil(minX / dx) + (N - 1) / 2;
						int iiMax = floor(maxX / dx) + (N - 1) / 2;
						int jjMin = ceil(minY / dy) + (M - 1) / 2;
						int jjMax = floor(maxY / dy) + (M - 1) / 2;
						if (iiMin < 0)	iiMin = 0;
						if (jjMin < 0)	jjMin = 0;
						register complex<double> tempEx;
						register complex<double> tempEy;
						for (int ii = iiMin; ii <= iiMax&&ii < N; ii++)
							for (int jj = jjMin; jj <= jjMax&&jj < M; jj++)
							{
								P.set(Plane[ii][jj].x, Plane[ii][jj].y);
								if (InterVal_PointinTriangle(A[0], A[1], A[2], P))//如果均匀格点位置落在三角型里面，然后在三角型里面进行线性插值
								{
									double a = (P - A[0]).area();
									double b = (P - A[1]).area();
									double c = (P - A[2]).area();

									tempEx = (
										b * c * Ex_R[dataIJ[0] / M][dataIJ[0] % M] +
										a * c * Ex_R[dataIJ[1] / M][dataIJ[1] % M] +
										b * a * Ex_R[dataIJ[2] / M][dataIJ[2] % M]) /
										(b * c + a * c + b * a);

									tempEy = (
										b * c * Ey_R[dataIJ[0] / M][dataIJ[0] % M] +
										a * c * Ey_R[dataIJ[1] / M][dataIJ[1] % M] +
										b * a * Ey_R[dataIJ[2] / M][dataIJ[2] % M]) /
										(b * c + a * c + b * a);
									if (abs(tempEx) + abs(tempEy) > 1e-5) {
										Ex_In[ii][jj] = tempEx;
										Ey_In[ii][jj] = tempEy;
									}
									else {
										ii = ii;
									}
								}//if
							}//two for
					}//for i  EleNum

			}
		}//omp
	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(80, user);
	}
	
	Free_2D(hit);	hit = NULL;

	fileAddress = "./PVVD/out1.dat";
	err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	tx = SourceGraphTrans.getTrans_x(); ty = SourceGraphTrans.getTrans_y(); tz = SourceGraphTrans.getTrans_z();
	rx = SourceGraphTrans.getRotate_x(); ry = SourceGraphTrans.getRotate_y(); rz = SourceGraphTrans.getRotate_z();
	rth = SourceGraphTrans.getRotate_theta();
	fwrite(&f, sizeof(double), 1, writefile);
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex_In[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey_In[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);
	

	RotateAsix = Vector3D(SourceGraphTrans.getRotate_x(),SourceGraphTrans.getRotate_y(),SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(
		SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(
		SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());

	// 设置源各个点的位置
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Plane[i][j] = translateMatrix * rotatMatrix * Plane[i][j];
		}
	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(90, user);
	}
}

void PVVA::CalReflectExyz(const Vector3 & n, const complex<double>& Ex,
	const complex<double>& Ey, const complex<double>& Ez, complex<double>& Ex_out,
	complex<double>& Ey_out, complex<double>& Ez_out)
{
	complex<double> ne(n.x * Ex.real() + n.y * Ey.real() + n.z * Ez.real(),
		n.x * Ex.imag() + n.y * Ey.imag() + n.z * Ez.imag());
	complex<double> tempx = 2 * n.x * ne;
	complex<double> tempy = 2 * n.y * ne;
	complex<double> tempz = 2 * n.z * ne;

	Ex_out = tempx - Ex;
	Ey_out = tempy - Ey;
	Ez_out = tempz - Ez;
}

double PVVA::CalSquare(const Vector3 & A, const Vector3 & B, const Vector3 & C, const Vector3 & D) const
{
	double AB = (B - A).Length();
	double AC = (C - A).Length();
	double AD = (D - A).Length();
	double BC = (C - B).Length();
	double DC = (C - D).Length();
	double p1 = (AB + AC + BC) / 2;
	double p2 = (AD + AC + DC) / 2;

	double tempS1 = sqrt(p1*(p1 - AB)*(p1 - BC)*(p1 - AC));
	double tempS2 = sqrt(p2*(p2 - AC)*(p2 - AD)*(p2 - DC));

	return tempS1 + tempS2;
}

//计算三角型面积
double PVVA::CalSquare(const Vector3 & A, const Vector3 & B,
	const Vector3 & C) const
{
	double AB = (B - A).Length();
	double AC = (C - A).Length();
	double BC = (B - C).Length();
	double p1 = (AB + AC + BC) / 2;
	
	return sqrt(p1*(p1 - AB)*(p1 - BC)*(p1 - AC));
}

bool PVVA::get_line_intersection(
	const Vector2 & A, const Vector2 & B, 
	const Vector2 & C, const Vector2 & D, Vector2 & O)
{
	double s10_x = B.x - A.x;
	double s10_y = B.y - A.y;
	double s32_x = D.x - C.x;
	double s32_y = D.y - C.y;

	double denom = s10_x * s32_y - s32_x * s10_y;
	if (abs(denom) < 0.000001)//平行或共线
		return 0; // Collinear
	bool denomPositive = denom > 0.0;

	double s02_x = A.x - C.x;
	double s02_y = A.y - C.y;
	double s_numer = s10_x * s02_y - s10_y * s02_x;
	if ((s_numer < 0.0) == denomPositive)
		//参数是大于等于0且小于等于1的，分子分母必须同号且分子小于等于分母
		return false; // No collision

	double t_numer = s32_x * s02_y - s32_y * s02_x;
	if ((t_numer < 0.0) == denomPositive)
		return false; // No collision

	if (((s_numer > denom) == denomPositive) || ((t_numer > denom) == denomPositive))
		return false; // No collision
				  // Collision detected
	double t = t_numer / denom;
	
	O.set(A.x + (t * s10_x), A.y + (t * s10_y));

	return true;
}

void PVVA::AllocateMemory()
{
	//坡印廷矢量
	Px = Allocate_2D(Px, N, M);
	Py = Allocate_2D(Py, N, M);
	Pz = Allocate_2D(Pz, N, M);

	//源
	Ex_In = Allocate_2D(Ex_In, N, M);
	Ey_In = Allocate_2D(Ey_In, N, M);

	//传播后的电磁场
	Ex1 = Allocate_2D(Ex1, N, M);
	Ey1 = Allocate_2D(Ey1, N, M);
	Ez1 = Allocate_2D(Ez1, N, M);

	Hx1 = Allocate_2D(Hx1, N, M);
	Hy1 = Allocate_2D(Hy1, N, M);
	Hz1 = Allocate_2D(Hz1, N, M);

	//反射后的电场
	Ex_R = Allocate_2D(Ex_R, N, M);
	Ey_R = Allocate_2D(Ey_R, N, M);
	Ez_R = Allocate_2D(Ez_R, N, M);

	// 平面位置（反射先，反射后）和各点的法向量
	n_Plane = Allocate_2D(n_Plane, N, M);
	Plane = Allocate_2D(Plane, N, M);
	R_Plane = Allocate_2D(R_Plane, N, M);
	Rn_Plane = Allocate_2D(Rn_Plane, N, M);

}

void PVVA::setSource(const FieldBase* _in) {
	const vector<vector<complex<double>>>&Ex_temp = _in->Ex;
	const vector<vector<complex<double>>>&Ey_temp = _in->Ey;

	SourceGraphTrans = _in->graphTransField;
	SourceGraphTrans.normalization(unit_factor);
	dx = _in->ds_x*unit_factor;	dy = _in->ds_y*unit_factor;
	N = _in->N_width;			M = _in->M_depth;
	AllocateMemory();

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Ex_In[i][j] = Ex_temp[i][j];
			Ey_In[i][j] = Ey_temp[i][j];
		}

	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());


	// 设置源各个点的位置
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
			tempPoint = translateMatrix * rotatMatrix * tempPoint;
			Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
		}

	// 设置源的方向
	Position3D tempPoint(0, 0, 1);
	tempPoint = rotatMatrix * tempPoint;
	n_Source.set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
	n_Source.Normalization();

	// 平面的中点
	Position3D tempPoint1(0, 0, 0);
	tempPoint1 = translateMatrix * rotatMatrix * tempPoint1;
	Org_Source.set(tempPoint1.X(), tempPoint1.Y(), tempPoint1.Z());

	//	updateSource_n(Vector3(1, 1, 0));
}

void PVVA::setOut(const FieldBase* _out) {
	TargetGraphTrans = _out->graphTransField;
}

void PVVA::getSourcePoint(Vector3 & interPoint, Vector3 & n_Point) const
{
	interPoint = Org_Source;
	n_Point = n_Source;
}

void PVVA::setSourcePoint(const Vector3 & interPoint, const Vector3 & n_Point)
{
	Org_Source = interPoint;
	n_Source = n_Point;
}

GraphTrans PVVA::getSourceGraphTrans(const Vector3 & n_Point)
{
	Vector3 a = n_Point;
	GraphTrans res;
	a.Normalization();
	if (a.x != 0 || a.y != 0 || a.z != 1)
	{
		Vector3 rotate_axis = Vector3(0, 0, 1).Cross(a); // 旋转轴
		double rotate_theta = acos(Vector3(0, 0, 1).Dot(a));
		rotate_theta = rotate_theta / Pi * 180;
	
		res.updateRotate(rotate_axis, rotate_theta);
	}
	else
	{
		res.updateRotate(Vector3(0, 0, 1), 0);
	}
	return res;
}

void PVVA::setFre(double _fre)
{
	f = _fre;
	lamda = C_Speed / f;

}

void PVVA::setMirror(Mirror * mirror)
{
	this->mirror = mirror;
}

bool PVVA::CheckModle(Vector3 & interPoint, Vector3 & n_Point)
{
	return false;
}


void PVVA::CalZ0Theta() 
//这个函数计算入射面中心法向与镜面的交点，
//注意，这一版里，出射口面的法向由输入定义
{
	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(10, user);
	}

	Vector3 InterPoint, tempReflect;

	int i = 0;
	RayTracing rayTracing(mirror);
	bool tempIsIntersect = false;
	//寻找交点
	rayTracing.calcReflect(Org_Source, n_Source, tempReflect,
		InterPoint, tempIsIntersect);
	if (!tempIsIntersect)
		return; // 没有相交

	//计算平移的距离（此时传播方向没有改变）
	//这里是计算从口面到镜面的距离
	z0 = CalDistance(Org_Source, InterPoint);
	//把口面场从入射场计算到镜前场
	CalPlane(z0);

	FILE* writefile;
	string fileAddress = "./PVVD/outbref.dat";
	int err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	int flag = 0;
	double tx = SourceGraphTrans.getTrans_x(); double ty = SourceGraphTrans.getTrans_y(); double tz = SourceGraphTrans.getTrans_z();
	double rx = SourceGraphTrans.getRotate_x(); double ry = SourceGraphTrans.getRotate_y(); double rz = SourceGraphTrans.getRotate_z();
	double rth = SourceGraphTrans.getRotate_theta();
	fwrite(&f, sizeof(double), 1, writefile);
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex1[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey1[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);
	//
	Vector3 U,V,N;
	TargetGraphTrans.updateRotate(Vector3(TargetGraphTrans.getRotate_x(),
		                                  TargetGraphTrans.getRotate_y(),
		                                  TargetGraphTrans.getRotate_z()),
								  TargetGraphTrans.getRotate_theta());
	TargetGraphTrans.getUVN(U,V,N);
	R_Source = N;//镜反射后虚拟面的法向
	R_Source.Normalization();
	//在镜前场上计算坡印廷矢量
	Poynting();
}

void PVVA::CalPlane(double dis)
{

	FFTDI cal(f, 0, 0,dis, N, M);
	cal.SetthreadNum(threadNum);
	cal.Setds(dx, dy);
	cal.SetInput(Ex_In, Ey_In);
	
	cal.StartCal();
	cal.output(Ex1, Ey1, Ez1, Hx1, Hy1, Hz1); // 计算平面传播

	for (int i = 0; i < N; i++)     // 更新平面位置
		for (int j = 0; j < M; j++)
		{
			Plane[i][j] = Plane[i][j] + n_Source * dis;
		}

	Org_Source = Org_Source + n_Source * dis;
	SourceGraphTrans.updateTranslate(Org_Source);
}

void PVVA::CalArbitraryPlane(GraphTrans _GT0, GraphTrans _GTt, int mode)
{
	//mode == 0 先旋转再传播 适合镜前传播
	//mode == 1 先传播再旋转 适合镜后传播
	GraphTrans GT0 = _GT0;
	GraphTrans GTt = _GTt;

	//将入射场进行旋转至与目标1一致的坐标系
	Vector3 U0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(1.0, 0.0, 0.0);
	Vector3 V0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(0.0, 1.0, 0.0);
	Vector3 N0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(0.0, 0.0, 1.0);
	Vector3 Ut = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(1.0, 0.0, 0.0);
	Vector3 Vt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(0.0, 1.0, 0.0);
	Vector3 Nt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(0.0, 0.0, 1.0);

	if (mode == 0) {//先旋转再传播
		if (Ut.Dot(Vector3(U0)) > 0.999 && Vt.Dot(Vector3(V0)) > 0.999) {
			//定义完全相同不用转啦
		}
		else {
			FFTASRotation* RotField = new FFTASRotation;
			RotField->SetParas(this->f, this->N, this->M, this->dx,this->dy);
			RotField->SetRotationParas(GT0, GTt);
			RotField->SetInput(Ex_In, Ey_In);
			RotField->PerformRotate();
			RotField->output(Ex_In, Ey_In);
			delete RotField;
		}

		//将入射场传播至虚拟面1，注意，这里是在虚拟面1所在的坐标系下
		Vector3 dis; dis = Vector3(GTt.getTrans_x() - GT0.getTrans_x(), GTt.getTrans_y() - GT0.getTrans_y(), GTt.getTrans_z() - GT0.getTrans_z());
		//请出虚拟面的旋转矩阵把全局坐标系下的dis转到虚拟面1坐标系

		double Up, Vp, Np;
		Up = Ut.Dot(dis);	Vp = Vt.Dot(dis);	Np = Nt.Dot(dis);
		FFTDI* FDI = new FFTDI(this->f, Up, Vp, Np, N, M);
		FDI->SetInput(Ex_In, Ey_In);
		FDI->Setds(dx,dy);
		FDI->StartCal();
		FDI->output(Ex1, Ey1, Ez1, Hx1, Hy1, Hz1);
		delete FDI;
	}
	else if (mode == 1) {
		//将入射场传播至虚拟面1，注意，这里是在虚拟面1所在的坐标系下
		Vector3 dis; dis = Vector3(GTt.getTrans_x() - GT0.getTrans_x(), GTt.getTrans_y() - GT0.getTrans_y(), GTt.getTrans_z() - GT0.getTrans_z());
		//请出虚拟面的旋转矩阵把全局坐标系下的dis转到虚拟面1坐标系

		double Up, Vp, Np;
		Up = U0.Dot(dis);	Vp = V0.Dot(dis);	Np = N0.Dot(dis);
		FFTDI* FDI = new FFTDI(this->f, Up, Vp, Np, N, M);
		FDI->SetInput(Ex_In, Ey_In);
		FDI->Setds(dx,dy);
		FDI->StartCal();
		FDI->outsource(Ex_In, Ey_In);
		delete FDI;

		if (Ut.Dot(Vector3(U0)) > 0.999 && Vt.Dot(Vector3(V0)) > 0.999) {
			//定义完全相同不用转啦
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < M; j++) {
					Ex1[i][j] = Ex_In[i][j];
					Ey1[i][j] = Ey_In[i][j];
				}
			}
		}
		else {
			FFTASRotation* RotField = new FFTASRotation;
			RotField->SetParas(this->f, this->N, this->M, this->dx,this->dy);
			RotField->SetRotationParas(GT0, GTt);
			RotField->SetInput(Ex_In, Ey_In);
			RotField->PerformRotate();
			RotField->output(Ex1, Ey1, Ez1);
			delete RotField;
		}
	}
}

void PVVA::ReflectCUDA() {

	//加入功能，在这里产生出射光线
	RayVec rayvec;

	//cout << "Calculating Reflecting surface" << endl;
	// 遍历每一个点 N * M * EleNum 
	Vector3 InterPoint;			//交点
	Vector3 Reflight, n_light; // 反射光线 和 法向量
	Vector3 n_light_Plane;  // 相对于平面的法向量

	// 绝对坐标系转换到源的坐标系（只求旋转）
	Vector3D RotateAsixSou(SourceGraphTrans.getRotate_x(),
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrixSou = Matrix4D::getRotateMatrix(
		-SourceGraphTrans.getRotate_theta(), RotateAsixSou);

	Vector3 tempa, tempb;

	double dir_t;  // 坡印廷矢量方向的系数
	double plane_t;
	double d1, d2; //虚拟路径
	double tempphase; // 相位
	complex<double> tempejphase; // = exp(j*phase)
	RayTracing rayTracing(mirror);
	rayTracing.setthreadNum(threadNum);
	//rayTracing.setthreadNum(1);
	bool isInter = false;
	int num = N*M;
	vector<Vector3> startPiont(num);
	vector<Vector3> direction(num);
	vector<Vector3> nomal(num);
	vector<Vector3> intersection(num);
	vector<bool> isIntersect(num);
	vector<float> port(num);
	for (int i = 0; i < N; i++)//光线
	{
		for (int j = 0; j < M; j++)//光线
		{
			startPiont[i*N + j] = Plane[i][j];
			direction[i*N + j] = n_Plane[i][j];
		}
	}

	int returnNum = 0;
	int returnScale = N / 10;
	if (returnFloat) {		returnFloat(20, user);	}

	//这里把入射光线与镜面的反射计算完了
	rayTracing.calcReflectBatch(startPiont, direction, nomal, intersection, isIntersect, port);

	//确定需考虑的光线个数
	int RayNum = 0;
	for (int i = 0; i < N; i+=10) {
		for (int j = 0; j < N; j+=10) {
			if (isIntersect[i*N + j]) RayNum++;
		}
	}
	rayvec.Rays.resize(RayNum);
	//添加光线追踪起点
	int countRay;
	countRay = 0;
	for (int i = 0; i < N; i+=10) {
		for (int j = 0; j < N; j+=10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, startPiont[i*N + j]);
				countRay++;
			}
		}
	}
	//添加到反射镜交点
	countRay = 0;
	for (int i = 0; i < N; i+=10) {
		for (int j = 0; j < N; j+=10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, intersection[i*N + j]);
				countRay++;
			}
		}
	}

	//计算反射光线与出射虚拟面的相交，即反射光线落到出射虚拟面上
	for (int i = 0; i < N; i++)
	{
		if (returnFloat) // 如果没有注册则不回调
		{
			if (i >= returnNum * returnScale)
			{
				returnFloat(30 + 3 * returnNum, user);
				returnNum++;
			}
		}
		for (int j = 0; j < M; j++)
		{	
			int index = j + i*M;
			Ex_R[i][j] = 0;
			Ey_R[i][j] = 0;
			Ez_R[i][j] = 0;
			R_Plane[i][j] = 100;
			//这里开始，计算
			if (!isIntersect[index])	continue;//没相交就不用算了

			n_light = nomal[index];
			InterPoint = intersection[index];
			dir_t = port[index];
			n_light.Normalization();

			// 将反射面法向量转化到源坐标系上(先转换到绝对坐标系)
			n_light_Plane = rotatMatrixSou * n_light;
			n_light_Plane.Normalization();  // 单位化
			
			// 只做极化变换 这里计算场的反射
			CalReflectExyz(n_light_Plane, Ex1[i][j], Ey1[i][j],
				Ez1[i][j], Ex_R[i][j], Ey_R[i][j], Ez_R[i][j]);

			// 反射光线
			Reflight = RayTracing::reflectLight(n_Plane[i][j], n_light);
			Rn_Plane[i][j] = Reflight;

			if (dir_t > 0.0000000001)  // 平面在反射面的前面
			{
				//这里，InterPoint是入射射线与镜面的交点
				//Reflight是镜面反射光线的方向
				//Org_Source是虚拟面的中心，注意反射前反射后相同，这里用到的是反射后的
				//R_Source是镜反射后虚拟面的法向
				//R_Plane[i][j]是镜面反射光线与镜后虚拟面的交点
				plane_t = IntersectPlane(InterPoint, Reflight,
					Org_Source, R_Source, R_Plane[i][j]);
				d2 = CalDistance(InterPoint, R_Plane[i][j]);
				d1 = CalDistance(InterPoint, Plane[i][j]);

				if (plane_t > 0.0)  // 虚拟面2在反射面的前面
					tempphase = -(d1 + d2) / lamda * 2 * Pi;
				else
					tempphase = -(d1 - d2) / lamda * 2 * Pi;
				tempejphase = complex <double>(cos(tempphase), sin(tempphase));
				Ex_R[i][j] = ComMul(Ex_R[i][j], tempejphase);  // 只做相位变换
				Ey_R[i][j] = ComMul(Ey_R[i][j], tempejphase);
				Ez_R[i][j] = ComMul(Ez_R[i][j], tempejphase);
			}
			else if (dir_t < -0.0000000001)   // 平面在反射面的后面
			{
				plane_t = IntersectPlane(InterPoint, Reflight,
					Org_Source, R_Source, R_Plane[i][j]);
				d2 = CalDistance(InterPoint, R_Plane[i][j]);
				d1 = CalDistance(InterPoint, Plane[i][j]);

				if (plane_t < 0.0)  // 虚拟面2在反射面的后面
					tempphase = (d1 + d2) / lamda * 2 * Pi;
				else
					tempphase = (d1 - d2) / lamda * 2 * Pi;
				tempejphase = complex <double>(cos(tempphase), sin(tempphase));
				Ex_R[i][j] = ComMul(Ex_R[i][j], tempejphase); // 只做相位变换
				Ey_R[i][j] = ComMul(Ey_R[i][j], tempejphase);
				Ez_R[i][j] = ComMul(Ez_R[i][j], tempejphase);
			}
			else  // 平面与反射面相交
			{
				R_Plane[i][j] = InterPoint;
			}
		}
	} // endloop

	//添加到反射镜交点
	countRay = 0;
	for (int i = 0; i < N; i+=10) {
		for (int j = 0; j < N; j+=10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, R_Plane[i][j]);
				countRay++;
			}
		}
	}

	rayvec.writeRays2BiFile("./pvvarays.dat");

	CalAmplitude();  // 只做幅度变换

	Matrix4D rotatMatrixSou1 = Matrix4D::getRotateMatrix(
		SourceGraphTrans.getRotate_theta(), RotateAsixSou);

	//源的传播方向改变-其实镜面反射后虚拟面作为新的指向
	n_Source = R_Source;
	n_Source.Normalization();
	//
	SourceGraphTrans.updateRotate(Vector3(TargetGraphTrans.getRotate_x(),
										  TargetGraphTrans.getRotate_y(),
										  TargetGraphTrans.getRotate_z()),
								  TargetGraphTrans.getRotate_theta());
	//updateSource_n(n_Source);

	Vector3D RotateAsixSou2(SourceGraphTrans.getRotate_x(),
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrixSou2 = Matrix4D::getRotateMatrix(
		-SourceGraphTrans.getRotate_theta(), RotateAsixSou2);
	rotatMatrixSou1 = rotatMatrixSou2 * rotatMatrixSou1;

	//极化旋转 - 虚拟面上的场，将全局坐标的场极化变化到虚拟面2上的场极化，局部。
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Vector3 tempReal(Ex_R[i][j].real(), Ey_R[i][j].real(), Ez_R[i][j].real());
			Vector3 tempImag(Ex_R[i][j].imag(), Ey_R[i][j].imag(), Ez_R[i][j].imag());
			tempReal = rotatMatrixSou1 * tempReal;
			tempImag = rotatMatrixSou1 * tempImag;
			Ex_R[i][j] = complex<double>(tempReal.x, tempImag.x);
			Ey_R[i][j] = complex<double>(tempReal.y, tempImag.y);
			Ez_R[i][j] = complex<double>(tempReal.z, tempImag.z);
		}
	}

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(60, user);
	}

}

FieldBase PVVA::getFieldBaseout() {
	FieldBase result;
	result.ds_x = dx;	result.ds_y = dy;
	result.graphTransField = TargetGraphTrans;
	result.N_width = N;	result.M_depth = M;
	result.Ex.resize(N);	result.Ey.resize(N);	result.Ez.resize(N);
	result.Hx.resize(N);	result.Hy.resize(N);	result.Hz.resize(N);
	for (int i = 0; i < N; i++) {
		result.Ex[i].resize(M);	result.Ey[i].resize(M);	result.Ez[i].resize(M);
		result.Hx[i].resize(M);	result.Hy[i].resize(M);	result.Hz[i].resize(M);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			result.Ex[i][j] = Ex1[i][j];	result.Ey[i][j] = Ey1[i][j];	result.Ez[i][j] = Ez1[i][j];
			result.Hx[i][j] = Hx1[i][j];	result.Hy[i][j] = Hy1[i][j];	result.Hz[i][j] = Hz1[i][j];
		}
	}
	return result;
}

complex<double> PVVA::ConjugateMul(const complex<double> &A, const complex<double> &B) const
{
	return complex<double>(A.real() * B.real() + A.imag() * B.imag(),
		-A.real() * B.imag() + A.imag() * B.real());
}

complex<double> PVVA::ComMul(const complex<double>& A, const complex<double>& B) const
{
	return complex<double>(A.real() * B.real() - A.imag() * B.imag(),
		A.real() * B.imag() + A.imag() * B.real());
}

double PVVA::CalDistance(const Vector3 &a, const Vector3 &b) const
{
	return pow(pow((a.x - b.x), 2) + pow((a.y - b.y), 2)+ pow((a.z - b.z), 2), 0.5);
}

bool PVVA::InterVal_PointinTriangle(const Vector2 & A, 
	const Vector2 & B, const Vector2 & C, const Vector2 & P)
{
	Vector2 v0 = C - A;
	Vector2 v1 = B - A;
	Vector2 v2 = P - A;

	double dot00 = v0.Dot(v0);
	double dot01 = v0.Dot(v1);
	double dot02 = v0.Dot(v2);
	double dot11 = v1.Dot(v1);
	double dot12 = v1.Dot(v2);

	double inverDeno = 1.0 / (dot00 * dot11 - dot01 * dot01);

	double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < -0.0001 || u > 1.00001) // if u out of range, return directly
	{
		return false;
	}

	double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < -0.0001 || v > 1.00001) // if v out of range, return directly
	{
		return false;
	}

	return u + v <= 1.00001;
}



