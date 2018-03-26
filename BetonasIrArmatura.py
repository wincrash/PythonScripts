#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:11:57 2018

@author: ruslan
"""
import vtk


# create a rendering window and renderer

pointSource=vtk.vtkPointSource()
pointSource.SetCenter(0,0,0)
pointSource.SetRadius(10)
pointSource.SetNumberOfPoints(1000)

pointSource.Update()

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
 
# create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
# create source
source = vtk.vtkSphereSource()
source.SetCenter(0,0,0)
source.SetRadius(5.0)

source.Update()
normal=vtk.vtkPolyDataNormals()
normal.SetInputConnection(source.GetOutputPort())
normal.ComputeCellNormalsOn ()
normal.Update()
poly=normal.GetOutput()
rez=pointSource.GetOutput()

state=vtk.vtkDoubleArray()
state.SetName("busena")
state.SetNumberOfComponents(1)
state.SetNumberOfTuples(rez.GetNumberOfPoints())



selectEnclosedPoints=vtk.vtkSelectEnclosedPoints()
selectEnclosedPoints.SetInputConnection(pointSource.GetOutputPort())
selectEnclosedPoints.SetSurfaceConnection(source.GetOutputPort());
selectEnclosedPoints.Update()
state = selectEnclosedPoints.GetOutput().GetPointData().GetArray("SelectedPoints");
state.SetName("busena")
print(selectEnclosedPoints.GetOutput())
selectEnclosedPoints.GetOutput().GetPointData().SetScalars(state)

 
write=vtk.vtkDataSetWriter()
write.SetFileName("/tmp/aaa.vtk")
write.SetInputData(selectEnclosedPoints.GetOutput())
write.Write()


# mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(selectEnclosedPoints.GetOutputPort())
mapper.SetScalarRange(state.GetRange())
# actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetPointSize(10)
 
# assign actor to the renderer
ren.AddActor(actor)
 
# enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()