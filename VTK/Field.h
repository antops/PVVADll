/*
*	created by liyun 2017/12/27
*   function 显示场 3D 2D 相位 幅度 dB
*   version 2.0 重构了保存数据的数据结构
*/
#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <fstream>
#include <complex>
#include <String>

#include "../../Util/GraphTrans.h"

#include "../../Util/Vector3.h"

class Field
{
public:
	Field();
	~Field();

	void init();

	void freeMemory();
	void allocateMemory();

	void setNM(int N, int M);
	void getNM(int& N, int &M) const;
	void setPlane(const GraphTrans & _graphTransPara, double _ds);

	void setField(complex<double> ** _Ex, complex<double> ** _Ey, complex<double> ** _Ez,
		complex<double> ** _Hx, complex<double> ** _Hy, complex<double> ** _Hz);

	const vector<vector<complex<double>>>& getEx() const;
	const vector<vector<complex<double>>>& getEy() const;

	void setDs(double _Ds);
	double getDs() const;

	//得到场的参数
	void getSourcePara(GraphTrans & _graphTransPara, 
		int & _N_width, int & _M_depth, double &_ds) const;

	double getWidth() const { return data[7]; }
	double getDepth() const { return data[8]; }

	void save(const std::string &fileName) const;

	void readData();
	errno_t readDataBinary();

	void setFileAddress(const string & file);

protected:
	GraphTrans graphTrans; // 旋转平移参数
	vector<double> data;

	// 场分量
	vector<vector<complex<double>>> Ex, Ey, Ez;
	vector<vector<complex<double>>> Hx, Hy, Hz;

	int N_width, M_depth;  //N_width = para[0] /ds

	double ds;
	string fileAddress;
	bool isSource;

};

#endif // !DRAW_H

