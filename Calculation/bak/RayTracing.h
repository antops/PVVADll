/*
*	created by liyun 2017/11/27
*   function 计算光线追踪 输入镜子的指针 输出反射光线和交点以及是否相交
*   version 1.0
*/

#ifndef RAYTRACING_H  
#define RAYTRACING_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include "../VTK/STLMirror.h"
#include "../../Util/Vector3.h"
#include "../../Util/Matrix4D.h"
#include "../../Util/Vector3D.h"
#include <vector>

namespace calculation
{
	class RayTracing
	{
	public:
		RayTracing(STLMirror* _mirror);

		~RayTracing();

		void setMirror(STLMirror* mirror);

		// 计算直线与mirror相交，并输出交点与法线
		void calcNormalOfLine_Mirror(const Vector3& startPiont,
			const Vector3& direction,
			Vector3 &normal, Vector3 &intersection,
			bool &isIntersect, double &t);

		// 输出反射光线和交点以及是否相交
		void calcReflect(const Vector3& startPiont,
			const Vector3& direction,
			Vector3 &reflex, Vector3 &intersection,
			bool &isIntersect);


		// 根据入射光线和法线求反射光线
		static Vector3 reflectLight(const Vector3& a, const Vector3& n);
		void SetupNodes(void);

	private:

		// 根据模型剖分数据计算反射
		void calcReflectByPolyDataBatch(const vector<vector<Vector3>>& startPiont,
			const vector<vector<Vector3>>& direction,
			vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
			vector<vector<bool>> &isIntersect);

		void calcReflectByPolyData(const Vector3& startPiont,
			const Vector3& direction, Vector3 &reflex,
			Vector3 &intersection,
			bool &isIntersect);

		void calcNormalOfLine_MirrorByPolyData(const Vector3& startPiont,
			const Vector3& direction,
			Vector3 &normal, Vector3 &intersection,
			bool &isIntersect, double &t);

		void calcNormalOfLine_MirrorByPolyData_omp(const Vector3& startPiont,
			const Vector3& direction,
			Vector3 &normal, Vector3 &intersection,
			bool &isIntersect, double &t);


		// 三角形与直线相交判断
		bool isIntersect(const Vector3 &orig, const Vector3 &dir,
			const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
			Vector3 &intersection, double &t);


		STLMirror* mirror;
		vector<Vector3> NodesXYZ1;
		vector<Vector3> NodesXYZ2;
		vector<Vector3> NodesXYZ3;
	};

}

#endif // RAYTRACING_H

