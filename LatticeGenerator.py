#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 08:46:33 2018

@author: ruslan
"""

import vtk
import numpy as np
import math
import random
from scipy.spatial import Voronoi

outputFileName="input.vtp"
RMIN=0.1
RMAX=0.1
VMIN=-0.00001
VMAX=0.000001
R_COUNT=1
V_COUNT=1000

OFFSET_X=math.sqrt(3)*RMAX    
OFFSET_Y=RMAX  
OFFSET_Z=RMAX     





STATIC_VELOCITY=0.1
MOVE_NEGATIVE_X=-100
MOVE_POSITIVE_X=100
MOVE_NEGATIVE_Y=0.2
MOVE_POSITIVE_Y=1.5
MOVE_NEGATIVE_Z=-100
MOVE_POSITIVE_Z=100







LINE_FORCE_N_LIMIT_MIN=1000
LINE_FORCE_N_LIMIT_MAX=10000

LINE_FORCE_T_LIMIT_MIN=1000
LINE_FORCE_T_LIMIT_MAX=10000


BOUNDS_X=1
BOUNDS_Y=RMAX*4
BOUNDS_Z=2
CELL_SIZE=RMAX*2.0

NX=int(math.floor(BOUNDS_X/CELL_SIZE))
NY=int(math.floor(BOUNDS_Y/CELL_SIZE))
NZ=int(math.floor(BOUNDS_Z/CELL_SIZE))
print("NX=",NX,"NY=",NY,"NZ=",NZ)
KIEKIS=int(NX*NY*NZ)


R=np.linspace(RMIN, RMAX, num=R_COUNT)
V=np.linspace(VMIN, VMAX, num=V_COUNT)

rad_seg=vtk.vtkDoubleArray()
rad_seg.SetName("UNIQUE_RADIUS")
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


part_fix=vtk.vtkIntArray()
part_fix.SetName("PARTICLE_FIX")
part_fix.SetNumberOfComponents(1)
part_fix.SetNumberOfTuples(KIEKIS)


cellsLines=vtk.vtkCellArray()


points = np.zeros((KIEKIS,3))
id=0
for x in range(NX):
    for y in range(NY):
        for z in range(NZ):
            ran=random.randint(0,R_COUNT-1)
            part_type.SetTuple1(id,ran)
            rad.SetTuple1(id,R[ran])               
            
            #pp.SetPoint(id,x*CELL_SIZE+RMAX,y*CELL_SIZE+RMAX,z*CELL_SIZE+RMAX)
            xx=x*CELL_SIZE+RMAX
            yy=y*math.sqrt(3.0)*RMAX+RMAX
            zz=z*CELL_SIZE+RMAX

         #   xx=x*math.sqrt(3.0)*RMAX+RMAX
        #    yy=y*math.sqrt(3.0)*RMAX+RMAX
       #     zz=z*math.sqrt(3.0)*RMAX+RMAX
            if((y%2)==1):
                xx=x*math.sqrt(3.0)*RMAX
            #if((z%2)==1):
             #   yy=y*math.sqrt(3.0)*RMAX

            vx=V[random.randint(0,V_COUNT-1)]
            vy=V[random.randint(0,V_COUNT-1)]
            vz=V[random.randint(0,V_COUNT-1)]            
            part_fix.SetTuple1(id,0)
            
            if(xx<MOVE_NEGATIVE_X):
                part_fix.SetTuple1(id,1)
                vx=-STATIC_VELOCITY
            if(xx>MOVE_POSITIVE_X):
                part_fix.SetTuple1(id,1)
                vx=STATIC_VELOCITY            

            if(yy<MOVE_NEGATIVE_Y):
                part_fix.SetTuple1(id,1)
                vy=-STATIC_VELOCITY
            if(yy>MOVE_POSITIVE_Y):
                part_fix.SetTuple1(id,1)
                vy=STATIC_VELOCITY            


            if(zz<MOVE_NEGATIVE_Z):
                part_fix.SetTuple1(id,1)
                vz=-STATIC_VELOCITY
            if(zz>MOVE_POSITIVE_Z):
                part_fix.SetTuple1(id,1)
                vz=STATIC_VELOCITY            
                
            vel.SetTuple3(id, vx,vy,vz)
            xx=xx+OFFSET_X
            yy=yy+OFFSET_Y
            zz=zz+OFFSET_Z
            points[id]=[xx,yy,zz]
            pp.SetPoint(id,xx,yy,zz)
            
            
            
            
            id=id+1
vor = Voronoi(points)  
ll=vor.ridge_points


state=vtk.vtkIntArray()
state.SetName("STATE")
state.SetNumberOfComponents(1)
force_N=vtk.vtkDoubleArray()
force_N.SetName("F_N_LIMIT")
force_N.SetNumberOfComponents(1)

force_T=vtk.vtkDoubleArray()
force_T.SetName("F_T_LIMIT")
force_T.SetNumberOfComponents(1)

for l in ll:
    if l[0]>=0 and l[1]>=0:
        #print l[0],l[1]
        diffx=points[l[0]][0]-points[l[1]][0]
        diffy=points[l[0]][1]-points[l[1]][1]
        diffz=points[l[0]][2]-points[l[1]][2]
        distance=math.sqrt(diffx*diffx+diffy*diffy+diffz*diffz)
        skirtumas=distance-2*RMAX
        if(skirtumas<=RMAX*0.001):
            cellsLines.InsertNextCell(2);
            cellsLines.InsertCellPoint(l[0]);
            cellsLines.InsertCellPoint(l[1]);
            #print("celles ",l[0],l[1])
            state.InsertNextTuple1(0)
            force_N.InsertNextTuple1(random.uniform(LINE_FORCE_N_LIMIT_MIN, LINE_FORCE_N_LIMIT_MAX))
            force_T.InsertNextTuple1(random.uniform(LINE_FORCE_T_LIMIT_MIN, LINE_FORCE_T_LIMIT_MAX))
            
    
poly.SetPoints(pp)
poly.SetLines(cellsLines);
poly.GetCellData().SetScalars(state)
poly.GetCellData().AddArray(force_N)
poly.GetCellData().AddArray(force_T)

poly.GetPointData().SetScalars(rad)
poly.GetPointData().AddArray(part_type)
poly.GetPointData().SetVectors(vel)
poly.GetPointData().AddArray(part_fix)
poly.GetFieldData().AddArray(rad_seg)
writer=vtk.vtkXMLPolyDataWriter()
writer.SetInputData(poly)
writer.SetFileName(outputFileName)
writer.Write()











