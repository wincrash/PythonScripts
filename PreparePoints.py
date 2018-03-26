#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 09:01:33 2018

@author: ruslan
"""

import vtk
import numpy as np
import math
import random
from scipy.spatial import Voronoi

inputFileName="final.vtp"
outputFileName="output.vtp"
R_COUNT=1
V_MIN=-0.000000001
V_MAX=0.000000001
V_COUNT=1000
vel_val=np.linspace(V_MIN, V_MAX, num=V_COUNT)




#cells lines
CreateCells=True

LINE_FORCE_N_LIMIT_MIN=1000
LINE_FORCE_N_LIMIT_MAX=1500

LINE_FORCE_T_LIMIT_MIN=1000
LINE_FORCE_T_LIMIT_MAX=1500


CrackGeneration=True

PROC_XMIN=0.0
PROC_XMAX=0.2
PROC_YMIN=0.495
PROC_YMAX=0.505





FixedVelocity=True
DIRECTION=1
PROC=1.2
VEL_MOVE=0.1









reader=vtk.vtkXMLPolyDataReader()
reader.SetFileName(inputFileName)
reader.Update()
poly=reader.GetOutput()
radius=poly.GetPointData().GetArray("RADIUS")
radius_range=radius.GetRange()
print (radius_range)

KIEKIS=poly.GetNumberOfPoints()

R=np.linspace(radius_range[0], radius_range[1], num=R_COUNT)
rad_seg=vtk.vtkDoubleArray()
rad_seg.SetName("UNIQUE_RADIUS")
rad_seg.SetNumberOfComponents(1)
rad_seg.SetNumberOfTuples(R_COUNT)
for x in range(R_COUNT):
    rad_seg.SetTuple1(x,R[x])

#ideti radius skaidyma

   
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

verts=vtk.vtkCellArray()
points = np.zeros((KIEKIS,3))
bounds=poly.GetBounds()
print(bounds)

for x in range(KIEKIS):
    ppp=poly.GetPoint(x)
    #ppp[0]=ppp[0]-bounds[0]
    #ppp[1]=ppp[1]-bounds[2]
    #ppp[2]=ppp[2]-bounds[4]
    poly.GetPoints().SetPoint(x,ppp[0]-bounds[0],ppp[1]-bounds[2],ppp[2]-bounds[4])
    points[x]=[ppp[0]-bounds[0],ppp[1]-bounds[2],ppp[2]]
    part_type.SetTuple1(x,0)
    part_fix.SetTuple1(x,0)
    ran1=random.randint(0,V_COUNT-1)
    ran2=random.randint(0,V_COUNT-1)
    ran3=random.randint(0,V_COUNT-1)
    vel.SetTuple3(x,vel_val[ran1],vel_val[ran2],vel_val[ran3])
    
    if( FixedVelocity is True):
        tarpelis=radius_range[1]*PROC
        if(DIRECTION==0):            
            if(ppp[0]<(bounds[0]+tarpelis)):
                vel.SetTuple3(x,-VEL_MOVE,0,0)
                part_fix.SetTuple1(x,1)
            if(ppp[0]>(bounds[1]-tarpelis)):
                vel.SetTuple3(x,VEL_MOVE,0,0)    
                part_fix.SetTuple1(x,1)
        if(DIRECTION==1):            
            if(ppp[1]<(bounds[2]+tarpelis)):
                vel.SetTuple3(x,0,-VEL_MOVE,0)
                part_fix.SetTuple1(x,1)
            if(ppp[1]>(bounds[3]-tarpelis)):
                vel.SetTuple3(x,0,VEL_MOVE,0)    
                part_fix.SetTuple1(x,1)
        if(DIRECTION==2):            
            if(ppp[2]<(bounds[4]+tarpelis)):
                vel.SetTuple3(x,0,0,VEL_MOVE)
                part_fix.SetTuple1(x,1)
            if(ppp[2]>(bounds[5]-tarpelis)):
                vel.SetTuple3(x,0,0,VEL_MOVE)    
                part_fix.SetTuple1(x,1)
    
    if(CreateCells is False):
        verts.InsertNextCell(1);
        verts.InsertCellPoint(x);
            

if(CreateCells is False):
    poly.SetVerts(verts)
else:
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
    
    cellsLines=vtk.vtkCellArray()
   # bounds=poly.GetBounds()
    #print(bounds)
    for l in ll:
        if l[0]>=0 and l[1]>=0:
            #print l[0],l[1]
            diffx=points[l[0]][0]-points[l[1]][0]
            diffy=points[l[0]][1]-points[l[1]][1]
            diffz=points[l[0]][2]-points[l[1]][2]
            distance=math.sqrt(diffx*diffx+diffy*diffy+diffz*diffz)
            r1=poly.GetPointData().GetArray("RADIUS").GetTuple1(l[0])
            r2=poly.GetPointData().GetArray("RADIUS").GetTuple1(l[1])
            skirtumas=distance-(r1+r2)
            if(skirtumas<=(r1+r2)*0.001):
                cellsLines.InsertNextCell(2);
                cellsLines.InsertCellPoint(l[0]);
                cellsLines.InsertCellPoint(l[1]);
                busena=0
                if(CrackGeneration is True):
                    xminas=(bounds[1]-bounds[0])*PROC_XMIN
                    xmaxas=(bounds[1]-bounds[0])*PROC_XMAX
                    yminas=(bounds[3]-bounds[2])*PROC_YMIN
                    ymaxas=(bounds[3]-bounds[2])*PROC_YMAX
                    xvid=(points[l[0]][0]+points[l[1]][0])/2.0
                    yvid=(points[l[0]][1]+points[l[1]][1])/2.0
                    #print (xvid,yvid,xminas,xmaxas,yminas,ymaxas)
                    if( xvid>xminas and xvid<xmaxas and yvid>yminas and yvid<ymaxas):                                          
                        print("aaa")
                        busena=1                
                
                state.InsertNextTuple1(busena)
                force_N.InsertNextTuple1(random.uniform(LINE_FORCE_N_LIMIT_MIN, LINE_FORCE_N_LIMIT_MAX))
                force_T.InsertNextTuple1(random.uniform(LINE_FORCE_T_LIMIT_MIN, LINE_FORCE_T_LIMIT_MAX))
                
                
                
                
    poly.SetLines(cellsLines)
    poly.GetCellData().SetScalars(state)
    poly.GetCellData().AddArray(force_N)
    poly.GetCellData().AddArray(force_T)


poly.GetPointData().AddArray(part_type)
poly.GetPointData().SetVectors(vel)
poly.GetPointData().AddArray(part_fix)
poly.GetFieldData().AddArray(rad_seg)


writer=vtk.vtkXMLPolyDataWriter()
writer.SetInputData(poly)
writer.SetFileName(outputFileName)
writer.Write()
