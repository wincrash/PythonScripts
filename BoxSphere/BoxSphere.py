import  vtk
import math
import numpy as np
from gengeo import *
import random

boxSize=1
sphereSize=0.4
R=0.04

box=vtk.vtkCubeSource()
box.SetBounds(0,boxSize,0,boxSize,0,boxSize)
box.Update()
sphere=vtk.vtkSphereSource()
sphere.SetCenter(boxSize/2.0,boxSize/2.0,boxSize/2.0)
sphere.SetRadius(sphereSize)
sphere.Update()
tes=vtk.vtkTriangleFilter()
tes.SetInputData(sphere.GetOutput())
tes.Update()
w=vtk.vtkSTLWriter()
w.SetInputData(tes.GetOutput())
w.SetFileName("particles.stl")
w.Update()


test=vtk.vtkTriangleFilter()
test.SetInputData(box.GetOutput())
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




def GenerateParticles():
   (mesh,minPoint,maxPoint)=STLMESH("particles.stl")
   mntable = MNTable3D (minPoint,maxPoint,2.5*R,1)
   
   packer = InsertGenerator3D (           
      R,R,100,1000,1.0e-6,True
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
   rad_seg.SetNumberOfTuples(1)
   rad_seg.SetTuple1(0,R)   
    
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
   
   
GenerateParticles()