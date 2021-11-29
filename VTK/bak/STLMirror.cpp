#include "STLMirror.h"
#include <vtkSTLReader.h>
#include <vtkTransformPolyDataFilter.h>

STLMirror::STLMirror()
{
}

STLMirror::~STLMirror()
{
}

void STLMirror::setNameFile(const std::string & file)
{
	fileName.clear();
	fileName = file;
	readData();
}

void STLMirror::readData()
{
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(fileName.c_str());
	reader->Update();
	polyDataSTL = reader->GetOutput();

}
