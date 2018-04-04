#!/usr/bin/python

import vtk
from gengeo import *


def STLMESH(filename="sample.stl"):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    bounds=data.GetBounds()
    minPoint = Vector3(bounds[0],bounds[2],bounds[4])
    maxPoint = Vector3(bounds[1],bounds[3],bounds[5])
    NUMBER_OF_BOUNDARIES = data.GetNumberOfCells()
    NUMBER_OF_BOUNDARIES_POINTS = data.GetNumberOfPoints()
    pyramidTris = TriPatchSet()
    for i in range(0,NUMBER_OF_BOUNDARIES):
        cell = data.GetCell(i)
        p1_id = cell.GetPointId(0)
        p2_id= cell.GetPointId(1)
        p3_id = cell.GetPointId(2)
        vpoint1=data.GetPoint(p1_id)
        vpoint2=data.GetPoint(p2_id)
        vpoint3=data.GetPoint(p3_id)
        p0 = Vector3(vpoint1[0],vpoint1[1],vpoint1[2])
        p1 = Vector3(vpoint2[0],vpoint2[1],vpoint2[2])
        p2 = Vector3(vpoint3[0],vpoint3[1],vpoint3[2])
        pyramidTris.addTriangle(p0,p1,p2,2)
    pyramid = MeshVolume(Mesh=pyramidTris)    
    return (pyramid,minPoint,maxPoint)        


def  GeneratePacking(meshSTL,minR,maxR,outputName,particle_id,bonds):
       (mesh,minPoint,maxPoint)=STLMESH(meshSTL)
       mntable = MNTable3D (
          minPoint,maxPoint,2.5*maxR,particle_id
       )
    
       packer = InsertGenerator3D (
          minR,maxR,10000,1000,1.0e-6,True
       )
    
       packer.generatePacking(mesh,mntable,0,1)
       if(bonds):
           mntable.generateBondsTagged(0,1.0e-6,0,1,1)
       mntable.write(outputName,2)
    
def  GeneratePackingSecond(meshSTL,minR,maxR,outputName,particle_id,bonds,OtherVTU):
       (mesh,minPoint,maxPoint)=STLMESH(meshSTL)
       mntable = MNTable3D (
          minPoint,maxPoint,2.5*maxR,particle_id
       )
       read=vtk.vtkXMLUnstructuredGridReader()
       read.SetFileName(OtherVTU)
       read.Update()
       
       
       for x in range(read.GetOutput().GetNumberOfPoints()):
           p=read.GetOutput().GetPoint(x)
           r=float(read.GetOutput().GetPointData().GetArray("radius").GetTuple1(x))
           tag=int(read.GetOutput().GetPointData().GetArray("particleTag").GetTuple1(x))
           xx=float(p[0])
           yy=float(p[1])
           zz=float(p[2])
           S=Sphere(Vector3(xx,yy,zz),float(r))
           mntable.insert(S,0)
           print(x)
           
           
           
       packer = InsertGenerator3D (
          minR,maxR,100,1000,1.0e-6,True
       )
    
       packer.generatePacking(mesh,mntable,0,1)
       #mntable.generateBondsTagged(0,1.0e-6,0,1,1)
       mntable.write(outputName,2)
    
    