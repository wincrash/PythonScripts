#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 09:16:12 2018

@author: ruslan
"""
import GenGeoModule
import vtk
import numpy
import random

box="box.stl"
betonas="beton.stl"
rodas="rod.stl"

RMAX=0.0005
RMIN=RMAX
STEEL_RMAX=RMAX
ITER=1000
v=0.01
procN=0.1
F_N_RANGES=[[1,0.5],[2000000000,2000000000],[2000000000,2000000000]]
F_T_RANGES=[[F_N_RANGES[0][0]*procN,F_N_RANGES[0][1]*procN],[2000000000,2000000000],[2000000000,2000000000]]
ConcreateMaterials=4

#GenGeoModule.GeneratePacking(rodas,STEEL_RMAX,STEEL_RMAX,"rod.vtu",1,False);
#GenGeoModule.GeneratePackingSecond(betonas,RMIN*0.5,RMAX,"beton.vtu",2,False,"rod.vtu");


GenGeoModule.GeneratePacking(betonas,RMAX,RMAX,"tempas.vtu",1,True);
reader=vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("tempas.vtu")
reader.Update()
rodReader = vtk.vtkSTLReader()
rodReader.SetFileName(rodas)
rodReader.Update()
v1=vtk.vtkSelectEnclosedPoints()
v1.SetInputConnection(reader.GetOutputPort())
v1.SetSurfaceConnection(rodReader.GetOutputPort())
v1.CheckSurfaceOn ()
v1.InsideOutOff ()
v1.SetTolerance (1.0e-6)
v1.Update()
v1insideArray = v1.GetOutput().GetPointData().GetArray("SelectedPoints");
particleId=[]
for x in range(v1insideArray.GetNumberOfTuples()):
    particleId.append(v1insideArray.GetTuple1(x))
poly=vtk.vtkPolyData()
points=vtk.vtkPoints()
cells=vtk.vtkCellArray()

particle_id=vtk.vtkDoubleArray()
particle_id.SetName("ID")
particle_id.SetNumberOfComponents(1)
radius=vtk.vtkDoubleArray()
radius.SetName("RADIUS")
radius.SetNumberOfComponents(1)
bonds=vtk.vtkIntArray()
bonds.SetName("BONDS_ID")
bonds.SetNumberOfComponents(1)
data=reader.GetOutput()
for x in range(data.GetNumberOfPoints()):
    points.InsertNextPoint(data.GetPoint(x))
    radius.InsertNextTuple1(data.GetPointData().GetArray("radius").GetTuple1(x))
    particle_id.InsertNextTuple1(particleId[x])
for x in range(data.GetNumberOfCells()):
    cell=data.GetCell(x)
    cells.InsertNextCell(2);
    cells.InsertCellPoint(cell.GetPointId(0));
    cells.InsertCellPoint(cell.GetPointId(1));
    v1=particle_id.GetTuple1(cell.GetPointId(0))
    v2=particle_id.GetTuple1(cell.GetPointId(1))
    bonds.InsertNextTuple1(0)
    
    if(v1==0 and v2==0):
        bonds.SetTuple1(x,0)
    if((v1==1 and v2==0) or (v1==0 and v2==1)):
       bonds.SetTuple1(x,1)
    if(v1==1 and v2==1):
        bonds.SetTuple1(x,2)


poly.SetPoints(points)
poly.GetPointData().SetScalars(radius)
poly.GetPointData().AddArray(particle_id)
poly.GetCellData().SetScalars(bonds)
poly.SetLines(cells)
poly.GetPointData().SetActiveScalars("ID")


writer=vtk.vtkXMLPolyDataWriter()
writer.SetFileName("final.vtp")
writer.SetInputData(poly)
writer.Write()










#######
reader=vtk.vtkXMLPolyDataReader()
reader.SetFileName("final.vtp")
reader.Update()
bounds=reader.GetOutput().GetBounds()
KIEKIS=reader.GetOutput().GetNumberOfPoints()
NumberOfCells=reader.GetOutput().GetNumberOfCells()

poly=vtk.vtkPolyData()
pp=vtk.vtkPoints();
pp.SetNumberOfPoints(KIEKIS)

rad_seg=vtk.vtkDoubleArray()
rad_seg.SetName("UNIQUE_RADIUS")
rad_seg.SetNumberOfComponents(1)
rad_seg.SetNumberOfTuples(1)
rad_seg.SetTuple1(0,RMIN)

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



for x in range(KIEKIS):
    r=reader.GetOutput().GetPointData().GetArray("RADIUS").GetTuple1(x)
    pid=reader.GetOutput().GetPointData().GetArray("ID").GetTuple1(x)
    p=reader.GetOutput().GetPoint(x)
    pp.SetPoint(x,p)
    rad.SetTuple1(x,r)
    vel.SetTuple3(x,0,0,0)
    part_type.SetTuple1(x,0)
    part_fix.SetTuple1(x,0)
    if(pid!=1):
        part_material.SetTuple1(x,random.randint(1,ConcreateMaterials))
    else:
        part_material.SetTuple1(x,0)
    if(p[1]<(bounds[2]+1.0*RMAX) and pid==1):
        vel.SetTuple3(x,0,0,0)
        part_fix.SetTuple1(x,2)
    if(p[1]>(bounds[3]-1.5*RMAX)and pid==1):
        vel.SetTuple3(x,0,v,0)
        part_fix.SetTuple1(x,4)

    if(p[0]<(bounds[0]+1.0*RMAX) and pid==1):
        vel.SetTuple3(x,0,0,0)
        part_fix.SetTuple1(x,1)
    if(p[0]>(bounds[1]-1.5*RMAX)and pid==1):
        vel.SetTuple3(x,0,0,0)
        part_fix.SetTuple1(x,1)
        
    if(p[2]<(bounds[4]+1.0*RMAX) and pid==1):
        vel.SetTuple3(x,0,0,0)
        part_fix.SetTuple1(x,3)
    if(p[2]>(bounds[5]-1.5*RMAX)and pid==1):
        vel.SetTuple3(x,0,0,0)
        part_fix.SetTuple1(x,3)

           



cellsLines=vtk.vtkCellArray()
state=vtk.vtkIntArray()
state.SetName("STATE")
state.SetNumberOfComponents(1)
state.SetNumberOfTuples(NumberOfCells)   
force_N=vtk.vtkDoubleArray()
force_N.SetName("F_N_LIMIT")
force_N.SetNumberOfComponents(1)
force_N.SetNumberOfTuples(NumberOfCells)   

force_T=vtk.vtkDoubleArray()
force_T.SetName("F_T_LIMIT")
force_T.SetNumberOfComponents(1)
force_T.SetNumberOfTuples(NumberOfCells)   

for x in range(NumberOfCells):
    cell=reader.GetOutput().GetCell(x)
    cellsLines.InsertNextCell(2);
    cellsLines.InsertCellPoint(cell.GetPointId(0));
    cellsLines.InsertCellPoint(cell.GetPointId(1));
    state.SetTuple1(x,0)
    bond_id=int(reader.GetOutput().GetCellData().GetArray("BONDS_ID").GetTuple1(x))
    
    force_T.SetTuple1(x,0)
    force_N.SetTuple1(x,random.uniform(F_N_RANGES[bond_id][0],F_N_RANGES[bond_id][1]))
    force_T.SetTuple1(x,random.uniform(F_T_RANGES[bond_id][0],F_T_RANGES[bond_id][1]))





poly.SetPoints(pp)
poly.GetPointData().SetScalars(rad)
poly.GetPointData().SetVectors(vel)
poly.GetPointData().AddArray(part_type)
poly.GetPointData().AddArray(part_fix)
poly.GetPointData().AddArray(part_material)
poly.GetFieldData().AddArray(rad_seg)

poly.SetLines(cellsLines)
poly.GetCellData().SetScalars(state)
poly.GetCellData().AddArray(force_N)
poly.GetCellData().AddArray(force_T)


writer=vtk.vtkXMLPolyDataWriter()
writer.SetFileName("input.vtp")
writer.SetInputData(poly)
writer.Write()

poly.GetPointData().SetActiveScalars("PARTICLE_MATERIAL")
th=vtk.vtkThreshold()
th.SetInputData(poly)
th.ThresholdByLower (0.5)
th.Update()
print th.GetOutput()


writer=vtk.vtkDataSetWriter()
writer.SetFileName("STEEL.vtk")
writer.SetInputData(th.GetOutput())
writer.Write()



stl = vtk.vtkSTLReader()
stl.SetFileName(box)
stl.Update()
mat_id=0
materials=vtk.vtkIntArray()
materials.SetName("MATERIAL_ID")
materials.SetNumberOfComponents(1)
velocity=vtk.vtkDoubleArray()
velocity.SetName("VELOCITY")
velocity.SetNumberOfComponents(3)
for x in range(0,stl.GetOutput().GetNumberOfCells()):
	materials.InsertNextTuple1(mat_id)
for x in range(0,stl.GetOutput().GetNumberOfPoints()):
	velocity.InsertNextTuple3(0.0,0.0,0.0)
stl.GetOutput().GetPointData().AddArray(velocity)
stl.GetOutput().GetCellData().AddArray(materials)
tran=vtk.vtkTransform()
tran.Scale(1,1,1)
tranfilter=vtk.vtkTransformFilter()
tranfilter.SetTransform(tran)
tranfilter.SetInputData(stl.GetOutput())
tranfilter.Update()
wr=vtk.vtkDataSetWriter()
wr.SetFileName("mesh.vtk")
wr.SetInputData(tranfilter.GetOutput())
wr.Write()

