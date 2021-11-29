#pragma once
 
#ifndef PPVADLL_H
#define PPVADLL_H

#include <string>
#include "../Util/FieldBase.h"

class Mirror;
class _declspec(dllexport) PVVADll
{

public:
	PVVADll();
	~PVVADll();

	void setInField(FieldBase * _fin);
	void setOutField(FieldBase * _fout);
	void setThreadNum(int _threadNum) { ThreadNum = _threadNum; }
	void setAnalyticalModelFile(std::string _mirrorfile);
	void setSTLModelFile(std::string _stlfile);
	void setCATRMode(bool in) { CATRMode = in; }
	FieldBase getFieldout();
	void calculate(double fre);
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// 注册回调函数

private:
	FieldBase FieldBaseIn;
	FieldBase FieldBaseOut;
	std::string inputFieldFile;
	std::string stlMirrorFile;
	void(*returnFloat)(float, void*);
	void *user; // 回调函数的类指针
	Mirror *mirrorptr;
	bool AnalyticalMirror;
	bool CATRMode;
	int ThreadNum;

};


#endif // PPVADLL_H