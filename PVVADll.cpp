#include "PVVADll.h"
#include "Calculation\PVVA.h"
#include "VTK\Field.h"
#include "../ClonedMirrorsAndRestriction/Mirror.h"
#include "../ClonedMirrorsAndRestriction/STLMirror.h"
#include "../ClonedMirrorsAndRestriction/Paraboloid.h"
#include "../ClonedMirrorsAndRestriction/Restriction.h"
#include "../ClonedMirrorsAndRestriction/MirrorFactory.h"
#include "../ClonedMirrorsAndRestriction/Ellipsoid.h"

PVVADll::PVVADll()
{
	this->returnFloat = NULL;
	this->user = NULL;
	ThreadNum = 4;
	CATRMode = false;
}

PVVADll::~PVVADll()
{
	mirrorptr->~Mirror();
}

void PVVADll::setInField(FieldBase* _in) {
	FieldBaseIn.Ex = _in->Ex;
	FieldBaseIn.Ey = _in->Ey;
	FieldBaseIn.N_width = _in->N_width;
	cout << FieldBaseIn.N_width << endl;
	FieldBaseIn.M_depth = _in->M_depth;
	FieldBaseIn.graphTransField = _in->graphTransField;
	FieldBaseIn.ds_x = _in->ds_x;
	FieldBaseIn.ds_y = _in->ds_y;
}

void PVVADll::setOutField(FieldBase* _out) {
	FieldBaseOut.graphTransField = _out->graphTransField;
}

void PVVADll::setAnalyticalModelFile(std::string _mirrorfile) {
	
	Json::Reader reader;
	Json::Value js;
	ifstream file(_mirrorfile);
	reader.parse(file, js);
	file.close();
	mirrorptr = MirrorFactory::getMirrorByJson(js["Mirror"][0], ".");
	mirrorptr->updateData();
	AnalyticalMirror = true;
}

void PVVADll::setSTLModelFile(std::string _stlfile) {
	mirrorptr = new STLMirror;
	mirrorptr->setNameFile(_stlfile);
}

void PVVADll::calculate(double fre) {
	PVVA pvva;
	if (mirrorptr->getMirrorsType() == STLMIRROR) { CATRMode = false; }

	else if (mirrorptr->getRestriction(0)->getType() == Restriction::RES_POLYGON) {
		CATRMode = true;
	}
	else CATRMode = false;
	// 设置单位
	pvva.SetReturnFloat(returnFloat, user);
	pvva.setUnit(1);
	pvva.SetThreadNum(ThreadNum);
	// 设置频率
	pvva.setFre(fre);
	// 读入源并分配内存
	pvva.setSource(&FieldBaseIn);
	pvva.setOut(&FieldBaseOut);

	pvva.setMirror(mirrorptr);
	pvva.SetCATRMode(CATRMode);

	pvva.CalZ0Theta();
	pvva.ReflectCUDA();
	//pvva.Reflect();
	pvva.InterVal(); //插值

	pvva.Result();

	FieldBaseOut = pvva.getFieldBaseout();
	//pvva.getField(&inputField);
}

FieldBase PVVADll::getFieldout() {
	return FieldBaseOut;
}

void PVVADll::SetReturnFloat(void(*returnFloat)(float, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}
