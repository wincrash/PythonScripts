#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 12:45:37 2018

@author: ruslan
"""

import  vtk
import math
import numpy as np
from gengeo import *
import random
plotis=0.05
ilgis=0.01
R=0.001
R_KLIUTIS=0.001
KLIUCIU_KIEKIS=100

RMAX=R
if(RMAX<R_KLIUTIS):
    RMAX=R_KLIUTIS

virsus=0.01
kampas=-10
aukstis=0.02
atstumas=0.04

def RotateTranslate(data,kampas):
    tran=vtk.vtkTransform()
    tran.RotateZ(kampas)
    tranfilter=vtk.vtkTransformFilter()
    tranfilter.SetTransform(tran)
    tranfilter.SetInputData(data)
    tranfilter.Update()
    bounds=tranfilter.GetOutput().GetBounds()
    tran=vtk.vtkTransform()
    tran.Translate(-bounds[0],-bounds[2]+aukstis,-bounds[4])
    tranfilter1=vtk.vtkTransformFilter()
    tranfilter1.SetTransform(tran)
    tranfilter1.SetInputData(tranfilter.GetOutput())
    tranfilter1.Update()
    return tranfilter1.GetOutput()


def STLMESH(filename="sample.stl"):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    bounds=data.GetBounds()
    minPoint = Vector3(bounds[0]+RMAX,bounds[2]+RMAX,bounds[4]+RMAX)
    maxPoint = Vector3(bounds[1]-RMAX,bounds[3]-RMAX,bounds[5]-RMAX)
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
#plane=RotateTranslate(plane.GetOutput(),kampas)



def GenerateParticles():
   (mesh,minPoint,maxPoint)=STLMESH("particles.stl")
   mntable = MNTable3D (minPoint,maxPoint,2.0*RMAX,1)
   taskai=[]   
   
   for x in range(len(taskai)):
       S=Sphere(Vector3(float(taskai[x][0]),float(taskai[x][1]),float(taskai[x][2])),float(taskai[x][3]))
       mntable.insert(S,0)
       yra=1


   packer = InsertGenerator3D (
      R,R,100,1000,1.0e-9,True
   )
   
   packer.generatePacking(mesh,mntable,0,1)
   mntable.write("tempas.vtu",2)
   reader=vtk.vtkXMLUnstructuredGridReader()
   reader.SetFileName("tempas.vtu")
   reader.Update()
   poly=vtk.vtkPolyData()
   poly.SetPoints(reader.GetOutput().GetPoints())
   KIEKIS=poly.GetNumberOfPoints()
   rad_seg=vtk.vtkDoubleArray()
   rad_seg.SetName("UNIQUE_RADIUS")
   rad_seg.SetNumberOfComponents(1)
   rad_seg.SetNumberOfTuples(2)
   rad_seg.SetTuple1(0,R)
   rad_seg.SetTuple1(1,R_KLIUTIS)
    
   rad=vtk.vtkDoubleArray()
   rad.SetName("RADIUS")
   rad.SetNumberOfComponents(1)
   rad.SetNumberOfTuples(KIEKIS)
   vel=vtk.vtkDoubleArray()
   vel.SetName("VELOCITY")
   vel.SetNumberOfComponents(3)
   vel.SetNumberOfTuples(KIEKIS)
    
   part_type=vtk.vtkIntArray()
   part_type.SetName("PARTICLE_TYPE")
   part_type.SetNumberOfComponents(1)
   part_type.SetNumberOfTuples(KIEKIS)
   
   part_material=vtk.vtkIntArray()
   part_material.SetName("PARTICLE_MATERIAL")
   part_material.SetNumberOfComponents(1)
   part_material.SetNumberOfTuples(KIEKIS)    
    
   part_fix=vtk.vtkIntArray()
   part_fix.SetName("PARTICLE_FIX")
   part_fix.SetNumberOfComponents(1)
   part_fix.SetNumberOfTuples(KIEKIS)   
   
   
   
   reader.GetOutput().GetPointData().GetArray("radius").SetName("RADIUS")
   poly.GetPointData().SetScalars(reader.GetOutput().GetPointData().GetArray("RADIUS"))
   verts=vtk.vtkCellArray()
   for x in range(poly.GetNumberOfPoints()):
       vel.SetTuple3(x,0,0,0)
       
       part_material.SetTuple1(x,0)       
       tag=reader.GetOutput().GetPointData().GetArray("particleTag").GetTuple1(x)
  
       
       if(tag==0):
           part_type.SetTuple1(x,1)
           part_fix.SetTuple1(x,1)
           print(tag)
       else:
           part_type.SetTuple1(x,0)
           part_fix.SetTuple1(x,0)
       verts.InsertNextCell(1)
       verts.InsertCellPoint(x)
   poly.SetVerts(verts)
   poly.GetPointData().AddArray(part_fix)
   poly.GetPointData().AddArray(part_material)
   poly.GetPointData().AddArray(vel)
   poly.GetPointData().AddArray(part_type)
   poly.GetFieldData().AddArray(rad_seg)
   www=vtk.vtkXMLPolyDataWriter()
   www.SetInputData(poly)
   www.SetFileName("input.vtp")
   www.Write()


poly=vtk.vtkPolyData()
points=vtk.vtkPoints()
cells=vtk.vtkCellArray()


points.InsertNextPoint(0,0,0)#0
points.InsertNextPoint(plotis,0,0)#1
points.InsertNextPoint(plotis,0,ilgis)#2
points.InsertNextPoint(0,0,ilgis)#3
points.InsertNextPoint(0,virsus,0)#4
points.InsertNextPoint(plotis,virsus,0)#5
points.InsertNextPoint(plotis,virsus,ilgis)#6
points.InsertNextPoint(0,virsus,ilgis)#7
cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(1);
cells.InsertCellPoint(2);
cells.InsertCellPoint(3);

cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(1);
cells.InsertCellPoint(5);
cells.InsertCellPoint(4);

cells.InsertNextCell(4);
cells.InsertCellPoint(2);
cells.InsertCellPoint(3);
cells.InsertCellPoint(7);
cells.InsertCellPoint(6);

cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(3);
cells.InsertCellPoint(7);
cells.InsertCellPoint(4);


cells.InsertNextCell(4);
cells.InsertCellPoint(4);
cells.InsertCellPoint(5);
cells.InsertCellPoint(6);
cells.InsertCellPoint(7);



cells.InsertNextCell(4);
cells.InsertCellPoint(1);
cells.InsertCellPoint(2);
cells.InsertCellPoint(6);
cells.InsertCellPoint(5);


poly.SetPoints(points)
poly.SetPolys(cells)
box=RotateTranslate(poly,kampas)


test=vtk.vtkTriangleFilter()
test.SetInputData(box)
test.Update()
w=vtk.vtkSTLWriter()
w.SetInputData(test.GetOutput())
w.SetFileName("particles.stl")
w.Write()


GenerateParticles()


bounds=box.GetBounds()




print(bounds)









poly1=vtk.vtkPolyData()
points=vtk.vtkPoints()
cells=vtk.vtkCellArray()


points.InsertNextPoint(0,0,0)#0
points.InsertNextPoint(plotis,0,0)#1
points.InsertNextPoint(plotis,0,ilgis)#2
points.InsertNextPoint(0,0,ilgis)#3
cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(1);
cells.InsertCellPoint(2);
cells.InsertCellPoint(3);
poly1.SetPoints(points)
poly1.SetPolys(cells)
mesh=RotateTranslate(poly,kampas)


bb=mesh.GetBounds()
polymesh=vtk.vtkPolyData()
points=vtk.vtkPoints()
cells=vtk.vtkCellArray()


points.InsertNextPoint(0,0,0)#0
points.InsertNextPoint(plotis+atstumas,0,0)#1
points.InsertNextPoint(plotis+atstumas,0,ilgis)#2
points.InsertNextPoint(0,0,ilgis)#3
points.InsertNextPoint(0,bb[3]+RMAX,0)#4
points.InsertNextPoint(plotis+atstumas,bb[3]+RMAX,0)#5
points.InsertNextPoint(plotis+atstumas,bb[3]+RMAX,ilgis)#6
points.InsertNextPoint(0,bb[3]+RMAX,ilgis)#7
cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(1);
cells.InsertCellPoint(2);
cells.InsertCellPoint(3);

cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(1);
cells.InsertCellPoint(5);
cells.InsertCellPoint(4);

cells.InsertNextCell(4);
cells.InsertCellPoint(2);
cells.InsertCellPoint(3);
cells.InsertCellPoint(7);
cells.InsertCellPoint(6);

cells.InsertNextCell(4);
cells.InsertCellPoint(0);
cells.InsertCellPoint(3);
cells.InsertCellPoint(7);
cells.InsertCellPoint(4);


cells.InsertNextCell(4);
cells.InsertCellPoint(4);
cells.InsertCellPoint(5);
cells.InsertCellPoint(6);
cells.InsertCellPoint(7);



cells.InsertNextCell(4);
cells.InsertCellPoint(1);
cells.InsertCellPoint(2);
cells.InsertCellPoint(6);
cells.InsertCellPoint(5);



points.InsertNextPoint(box.GetPoint(0))#0
points.InsertNextPoint(box.GetPoint(1))#0
points.InsertNextPoint(box.GetPoint(2))#0
points.InsertNextPoint(box.GetPoint(3))#0
cells.InsertNextCell(4);
cells.InsertCellPoint(8);
cells.InsertCellPoint(9);
cells.InsertCellPoint(10);
cells.InsertCellPoint(11);


polymesh.SetPoints(points)
polymesh.SetPolys(cells)



test=vtk.vtkTriangleFilter()
test.SetInputData(polymesh)
test.Update()
vel=vtk.vtkDoubleArray()
vel.SetName("VELOCITY")
vel.SetNumberOfComponents(3)
KIEKIS=test.GetOutput().GetNumberOfPoints()
vel.SetNumberOfTuples(KIEKIS)
for x in range(KIEKIS):
    vel.SetTuple3(x,0,0,0)
test.GetOutput().GetPointData().AddArray(vel)
w=vtk.vtkDataSetWriter()
w.SetInputData(test.GetOutput())
w.SetFileName("mesh.vtk")
w.Write()
