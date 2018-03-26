#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 09:21:54 2018

@author: ruslan
"""

import vtk
import numpy as np
import math
import random

outputFileName="points.vtp"

NX=10
NY=100
R=0.016

poly=vtk.vtkPolyData()
points=vtk.vtkPoints()
radius=vtk.vtkDoubleArray()
radius.SetNumberOfComponents(1)
radius.SetName("RADIUS")
id=0
for i in range(NX):
    for k in range(NY):
        xx=i*2.0*R
        yy=k*2.0*R
        x=xx-yy*0.5      
        y=yy*math.sqrt(3.0)/2.0
        z=0
        radius.InsertNextTuple1(R)
        points.InsertNextPoint(x,y,z)
        id=id+1





poly.SetPoints(points)
poly.GetPointData().SetScalars(radius)



writer=vtk.vtkXMLPolyDataWriter()
writer.SetInputData(poly)
writer.SetFileName(outputFileName)
writer.Write()