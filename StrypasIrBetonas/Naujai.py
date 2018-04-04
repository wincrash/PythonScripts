#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 15:59:55 2018

@author: ruslan
"""

import vtk
import math
from subprocess import call
from gengeo import *
import random
import numpy as np
scale=1.0/1000.0
R=20
h=1.3
a=math.sqrt(3.0)*h
b=2
c=math.sqrt(3.0)*h
d=8


resolution=20

COUNT=3
boxSize=50
CalcBoxSizeProc=20





def GetForma(pradzia,points):
    startas=pradzia    
    points.append([R/2.0,startas])
    startas=startas+d/2
    points.append([R/2.0,startas])
    startas=startas+a
    points.append([h+R/2.0,startas])
    startas=startas+b
    points.append([h+R/2.0,startas])
    startas=startas+c
    points.append([R/2.0,startas])
    startas=startas+d/2
    points.append([R/2.0,startas])
    return startas,points


def SaveToSTL(data,name):
    tt=vtk.vtkTriangleFilter()
    tt.SetInputData(data)
    tt.Update()
    tran=vtk.vtkTransform()
    tran.Scale(scale,scale,scale)
    ttt=vtk.vtkTransformFilter()
    ttt.SetTransform(tran)
    ttt.SetInputData(tt.GetOutput())
    ttt.Update()
    w=vtk.vtkSTLWriter()
    w.SetFileName(name)
    w.SetInputData(ttt.GetOutput())
    w.Write()
    
    
def GetTriangles(data):
    tt=vtk.vtkTriangleFilter()
    tt.SetInputData(data)
    tt.Update()
    return tt.GetOutput()

def create3D():
    
    points=[]
    
    ilgis=0
    for x in range(COUNT):
    	ilgis,points=GetForma(ilgis,points)
    
    points.append([0,ilgis])
    points.append([0,0])
    points.append([R/2.0,0])
    rodData="translate([0,"+str(R/4.0)+","+str(R/4.0)+"])rotate ([-90,0,0])  rotate_extrude($fn=200) polygon( points="+str(points)+" );\n";
    rod=open("rod.scad","w")
    rod.write(rodData)
    rod.close()
    call(["openscad","-o","rod.stl","rod.scad"])
    read=vtk.vtkSTLReader()
    read.SetFileName("rod.stl")
    read.Update()
    tran=vtk.vtkTransform()
    bounds=read.GetOutput().GetBounds()
    tran.Translate(-read.GetOutput().GetBounds()[0]-(bounds[1]-bounds[0])/2.0,-read.GetOutput().GetBounds()[2]-(bounds[3]-bounds[2])/2.0,-read.GetOutput().GetBounds()[4]-(bounds[5]-bounds[4])/2.0)    
    tfilter=vtk.vtkTransformFilter()
    tfilter.SetInputData(read.GetOutput())
    tfilter.SetTransform(tran)
    tfilter.Update()
    
    w=vtk.vtkSTLWriter()
    w.SetFileName("rod.stl")
    w.SetInputData(tfilter.GetOutput())
    w.Update()
    
    plane1=vtk.vtkPlane()
    plane1.SetNormal(1,0,0)
    plane1.SetOrigin(0,0,0)
    plane2=vtk.vtkPlane()
    plane2.SetNormal(0,0,1)
    plane2.SetOrigin(0,0,0)      
        

    coll=vtk.vtkPlaneCollection()
    coll.AddItem(plane1)
    coll.AddItem(plane2)

    clip1=vtk.vtkClipClosedSurface()
    clip1.SetClippingPlanes(coll)
    clip1.GenerateFacesOn ()
    clip1.SetInputConnection(tfilter.GetOutputPort())
    clip1.Update()
    beton=vtk.vtkCubeSource()
    beton.SetBounds(-boxSize,boxSize,-ilgis/2.0,ilgis/2.0,-boxSize,boxSize);
    beton.Update()
    clip2=vtk.vtkClipClosedSurface()
    clip2.SetClippingPlanes(coll)
    clip2.GenerateFacesOn ()
    clip2.SetInputConnection(beton.GetOutputPort())
    clip2.Update()
    
    box=vtk.vtkCubeSource()
    b=clip2.GetOutput().GetBounds()
    bbbb=clip1.GetOutput().GetBounds()
    box.SetBounds(b[0],b[1],bbbb[2]-(bbbb[3]-bbbb[2])*CalcBoxSizeProc/1000.0,bbbb[3]+(bbbb[3]-bbbb[2])*CalcBoxSizeProc/1000.0,b[4],b[5])
    box.Update()
    SaveToSTL(box.GetOutput(),"box.stl")
    SaveToSTL(clip2.GetOutput(),"beton.stl")
    SaveToSTL(clip1.GetOutput(),"rod.stl")
    

    
        
        # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
     
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(clip1.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    ren.AddActor(actor)
    mapper1 = vtk.vtkPolyDataMapper()
    mapper1.SetInputConnection(clip2.GetOutputPort())
    actor1 = vtk.vtkActor()
    actor1.SetMapper(mapper1)
    ren.AddActor(actor1)
     
    # enable user interface interactor
    #iren.Initialize()
    #renWin.Render()
    #iren.Start()
    
create3D()
#GenerateParticles()
