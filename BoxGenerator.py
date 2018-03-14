#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:21:33 2018

@author: ruslan
"""

import vtk
import numpy as np
import math
import random

outputFileName="input.vtk"
RMIN=0.001
RMAX=0.01
VMIN=-0.1
VMAX=0.1
R_COUNT=10
V_COUNT=1000

BOUNDS_X=1.0
BOUNDS_Y=1.0
BOUNDS_Z=1.0
CELL_SIZE=RMAX*2.0

NX=int(math.floor(BOUNDS_X/CELL_SIZE))
NY=int(math.floor(BOUNDS_Y/CELL_SIZE))
NZ=int(math.floor(BOUNDS_Z/CELL_SIZE))
print("NX=",NX,"NY=",NY,"NZ=",NZ)
KIEKIS=int(NX*NY*NZ)


R=np.linspace(RMIN, RMAX, num=R_COUNT)
V=np.linspace(VMIN, VMAX, num=V_COUNT)





rad_seg=vtk.vtkDoubleArray()
rad_seg.SetName("RADIUS")
rad_seg.SetNumberOfComponents(1)
rad_seg.SetNumberOfTuples(R_COUNT)
for x in range(R_COUNT):
    rad_seg.SetTuple1(x,R[x])




poly=vtk.vtkPolyData()
pp=vtk.vtkPoints();
pp.SetNumberOfPoints(KIEKIS)
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
cells=vtk.vtkCellArray()
id=0
for x in range(NX):
    for y in range(NY):
        for z in range(NZ):
            ran=random.randint(0,R_COUNT-1)
            part_type.SetTuple1(id,ran)
            rad.SetTuple1(id,R[ran])   
            vel.SetTuple3(id, V[random.randint(0,V_COUNT-1)],V[random.randint(0,V_COUNT-1)],V[random.randint(0,V_COUNT-1)])
            pp.SetPoint(id,x*CELL_SIZE+RMAX,y*CELL_SIZE+RMAX,z*CELL_SIZE+RMAX)
            cells.InsertNextCell(1);
            cells.InsertCellPoint(id);
            id=id+1
    
poly.SetPoints(pp)
poly.SetVerts(cells);
poly.GetPointData().SetScalars(rad)
poly.GetPointData().AddArray(part_type)
poly.GetPointData().SetVectors(vel)
poly.GetFieldData().AddArray(rad_seg)
writer=vtk.vtkDataSetWriter()
writer.SetInputData(poly)
writer.SetFileName(outputFileName)
writer.Write()
