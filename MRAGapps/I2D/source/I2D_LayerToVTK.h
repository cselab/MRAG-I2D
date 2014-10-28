#pragma once

#include "I2D_Layer.h"

#include <vtkPoints.h> 
#include <vtkCell.h>
#include <vtkImageData.h>
#include <vtkImageNoiseSource.h>
#include <vtkFloatArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

template <int sizeX,int sizeY>
void dumpLayer2VTK(int istep, string sFileNamePattern, const I2D_Layer< sizeX,  sizeY, 1>& scalar_field)
{
	char buf[500];
	
	sprintf(buf, "%s_%07d.vti", sFileNamePattern.c_str(), istep);
	string sFileName(buf);
	
	vtkImageData * grid = vtkImageData::New();
	
	grid->SetExtent(0,sizeX-1,0,sizeY-1,0,0);
	grid->SetDimensions(sizeX, sizeY, 1);
	grid->SetScalarTypeToFloat();
	grid->SetNumberOfScalarComponents(1);
	grid->AllocateScalars();
	grid->SetSpacing(1./sizeX, 1./sizeX, 1);
	grid->SetOrigin(0,0,0);
	
	int iz = 0;
	for(int iy=0; iy<sizeY; iy++)
		for(int ix=0; ix<sizeX; ix++)
			grid->SetScalarComponentFromFloat(ix,iy,iz,0,scalar_field.read(ix,iy,0));
	
	vtkXMLImageDataWriter * writer = vtkXMLImageDataWriter::New();
	writer->SetFileName(sFileName.c_str());
	writer->SetInput(grid);
	writer->Write();
	
	writer->Delete();
	grid->Delete();
}

template <int sizeX,int sizeY>
void dumpLayer2VTK(int istep, string sFileNamePattern, const I2D_Layer< sizeX,  sizeY, 2>& scalar_field)
{
	char buf[500];
	
	sprintf(buf, "%s_%07d.vti", sFileNamePattern.c_str(), istep);
	string sFileName(buf);
	
	vtkImageData * grid = vtkImageData::New();
	
	grid->SetExtent(0,sizeX-1,0,sizeY-1,0,0);
	grid->SetDimensions(sizeX, sizeY, 1);
	grid->SetScalarTypeToFloat();
	grid->SetNumberOfScalarComponents(2);
	grid->AllocateScalars();
	grid->SetSpacing(1./sizeX, 1./sizeX, 1);
	grid->SetOrigin(0,0,0);
	
	int iz = 0;
	for(int iy=0; iy<sizeY; iy++)
		for(int ix=0; ix<sizeX; ix++)
		{
			grid->SetScalarComponentFromFloat(ix,iy,iz,0,scalar_field.read(ix,iy,0));
			grid->SetScalarComponentFromFloat(ix,iy,iz,1,scalar_field.read(ix,iy,1));
		}
	
	vtkXMLImageDataWriter * writer = vtkXMLImageDataWriter::New();
	writer->SetFileName(sFileName.c_str());
	writer->SetInput(grid);
	writer->Write();
	
	writer->Delete();
	grid->Delete();
}
