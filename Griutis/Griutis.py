#!/usr/bin/python

import  vtk
import math
import numpy as np
from gengeo import *
import random

plotis=0.05
ilgis=0.03
R=0.001
R_KLIUTIS=0.002
KLIUCIU_KIEKIS=100

RMAX=R
if(RMAX<R_KLIUTIS):
    RMAX=R_KLIUTIS



virsus=10*RMAX
kampas=-0
aukstis=0.01
atstumas=0.02
trans=[0,0,0]


def RotateTranslate(data):
    tran=vtk.vtkTransform()
    tran.RotateZ(kampas)
    tranfilter=vtk.vtkTransformFilter()
    tranfilter.SetTransform(tran)
    tranfilter.SetInputData(data)
    tranfilter.Update()
    bounds=tranfilter.GetOutput().GetBounds()
    tran=vtk.vtkTransform()
    tran.Translate(-bounds[0],-bounds[2]+aukstis,-bounds[4])
    trans[0]=-bounds[0]
    trans[1]=-bounds[2]+aukstis
    trans[2]=-bounds[4]
    tranfilter1=vtk.vtkTransformFilter()
    tranfilter1.SetTransform(tran)
    tranfilter1.SetInputData(tranfilter.GetOutput())
    tranfilter1.Update()
    return tranfilter1.GetOutput()

def RotateTranslate1(data,translate):
    tran=vtk.vtkTransform()
    tran.RotateZ(kampas)
    tranfilter=vtk.vtkTransformFilter()
    tranfilter.SetTransform(tran)
    tranfilter.SetInputData(data)
    tranfilter.Update()
    
    tran=vtk.vtkTransform()
    tran.Translate(translate[0],translate[1],translate[2])
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

def getKliutis(translate):
    x = np.linspace(0, plotis, math.ceil(plotis/(2.0*R_KLIUTIS)), endpoint=True)
    z = np.linspace(0, ilgis, math.ceil(ilgis/(2.0*R_KLIUTIS)), endpoint=True)
    ran_list=set()
    all_list=[]
    nn=0
    aaaaa=int(math.ceil(plotis/(2.0*R_KLIUTIS))-2)
    bbbbb=int(math.ceil(ilgis/(2.0*R_KLIUTIS))-2)
    for xx in range(2,aaaaa):
        for yy in range(2,bbbbb):
            all_list.append([xx,yy])
            nn=nn+1
    
    
    for xx in range(KLIUCIU_KIEKIS):
        ran1=random.randint(0,nn-1)        
        ran_list.add(ran1)
        
    new_list = []
    for xx in ran_list:
        new_list.append(all_list[xx])    
    
    
    #print(new_list)
    print("iss "+str(nn)+" gavom "+str(len(new_list)))
    poly=vtk.vtkPolyData()
    points=vtk.vtkPoints()
    for xx in new_list:
        points.InsertNextPoint(xx[0]*2.0*R_KLIUTIS,R_KLIUTIS,xx[1]*2.0*R_KLIUTIS)
    poly.SetPoints(points)
    
    
    poly=RotateTranslate1(poly,translate)
    taskai=[]
    for x in range(poly.GetNumberOfPoints()):
        p=poly.GetPoints().GetPoint(x)
        taskai.append([p[0],p[1],p[2],R_KLIUTIS])
    
    
    
    taskai=[]
    return taskai


def GenerateParticles(translate):
   (mesh,minPoint,maxPoint)=STLMESH("particles.stl")
   mntable = MNTable3D (minPoint,maxPoint,2.0*RMAX,1)
   taskai=getKliutis(translate)   
   
   for x in range(len(taskai)):
       S=Sphere(Vector3(float(taskai[x][0]),float(taskai[x][1]),float(taskai[x][2])),float(taskai[x][3]))
       mntable.insert(S,0)
       yra=1


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

poly.SetPoints(points)
poly.SetPolys(cells)


part=vtk.vtkCubeSource()
part.SetBounds(0,plotis,0,virsus,0,ilgis)
part.Update()


poly=RotateTranslate(poly)




poly1=RotateTranslate(part.GetOutput())
triangles=vtk.vtkTriangleFilter()
triangles.SetInputData(poly1)
triangles.Update()
w=vtk.vtkSTLWriter()
w.SetInputData(triangles.GetOutput())
w.SetFileName("particles.stl")
w.Write()






bounds=poly.GetBounds()
boxBounds=[]
boxBounds.append(bounds[0]-RMAX)
boxBounds.append(bounds[1]+atstumas+RMAX)
boxBounds.append(-RMAX)
boxBounds.append(bounds[3]+RMAX)
boxBounds.append(bounds[4]-RMAX)
boxBounds.append(bounds[5]+RMAX)
cube=vtk.vtkCubeSource()
cube.SetBounds(boxBounds)
cube.Update()

append=vtk.vtkAppendPolyData()
append.AddInputData (poly)
append.AddInputData (cube.GetOutput())
append.Update()

#mapper1 = vtk.vtkPolyDataMapper()
#mapper1.SetInputData(cube.GetOutput())
#actor1 = vtk.vtkActor()
#actor1.SetMapper(mapper1)
#actor1.GetProperty().SetOpacity(0.3)


GenerateParticles(trans)

w=vtk.vtkDataSetWriter()
w.SetInputData(append.GetOutput())
w.SetFileName("mesh.vtk")
w.Write()

read=vtk.vtkXMLPolyDataReader()
read.SetFileName("input.vtp")
ttt=vtk.vtkTransform()
ttt.Translate(RMAX,RMAX,RMAX)
ttfil=vtk.vtkTransformFilter()
ttfil.SetTransform(ttt)
ttfil.SetInputConnection(read.GetOutputPort())
w=vtk.vtkXMLPolyDataWriter()
w.SetFileName("input.vtp")
w.SetInputConnection(ttfil.GetOutputPort())
w.Write()

read=vtk.vtkDataSetReader()
read.SetFileName("mesh.vtk")
ttt=vtk.vtkTransform()
ttt.Translate(RMAX,RMAX,RMAX)
ttfil=vtk.vtkTransformFilter()
ttfil.SetTransform(ttt)
ttfil.SetInputConnection(read.GetOutputPort())
ttfil.Update()

vel=vtk.vtkDoubleArray()
vel.SetName("VELOCITY")
vel.SetNumberOfComponents(3)
KIEKIS=ttfil.GetOutput().GetNumberOfPoints()
vel.SetNumberOfTuples(KIEKIS)
for x in range(KIEKIS):
    vel.SetTuple3(x,0,0,0)
ttfil.GetOutput().GetPointData().AddArray(vel)

test=vtk.vtkTriangleFilter()
test.SetInputData(ttfil.GetOutput())
test.Update()
bbbb=test.GetOutput().GetBounds()
ttt=vtk.vtkTransform()
ttt.Translate(-bbbb[0],-bbbb[2],-bbbb[4])
ttfil=vtk.vtkTransformFilter()
ttfil.SetTransform(ttt)
ttfil.SetInputConnection(test.GetOutputPort())
ttfil.Update()



w=vtk.vtkDataSetWriter()
w.SetFileName("mesh.vtk")
w.SetInputData(ttfil.GetOutput())
w.Write()



