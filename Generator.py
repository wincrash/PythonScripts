#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 06:30:46 2018

@author: ruslan
"""

import vtk
import math
from subprocess import call
from gengeo import *
import random
import numpy as np

class Strypas:
    RMIN=0.8
    RMAX=0.8
    R_INER=10
    R=20
    E_steel_proc=3.14*R_INER*R_INER/(3.14*R*R)
    
    
    h=1.3
    a=math.sqrt(3)*h
    b=2
    c=math.sqrt(3)*h
    d=8
    COUNT=4
    boxSize=100
    AXIS1=[]
    AXIS2=[]
    ilgis=0
    
    contour2D=0
    contour3D=0
    betonas2D=0
    betonas3D=0
    box2D=0
    box3D=0
    BOX_RMAX_N=4
    particles=0
    
    v=0.02
    F_N_RANGES=[[0,0],[0.5,0.2],[1.0E+12,1.0E+12],[1.0E+12,1.0E+12]]    
    F_T_RANGES=[[0,0],[0.05,0.02],[1.0E+12,1.0E+12],[1.0E+12,1.0E+12]]    
    
    
    Max_material_count=5
    MATERIAL_PASIMETYMAS_PROC=0.1
    E=3.68E+07
    MIU=0.15
    RO=2700
    E_STEEL=2.03E+11*E_steel_proc
    MIU_STEEL=0.3
    RO_STEEL=8050

    def generateStrypas(self):
        for x in range(self.COUNT):
            self.AXIS1.append(self.R/2.0)
            self.AXIS2.append(self.ilgis)

            self.ilgis=self.ilgis+self.d/2.0
            self.AXIS1.append(self.R/2.0)
            self.AXIS2.append(self.ilgis)

            self.ilgis=self.ilgis+self.a
            self.AXIS1.append(self.h+self.R/2.0)
            self.AXIS2.append(self.ilgis)

            self.ilgis=self.ilgis+self.b
            self.AXIS1.append(self.h+self.R/2.0)
            self.AXIS2.append(self.ilgis)

            self.ilgis=self.ilgis+self.c
            self.AXIS1.append(self.R/2.0)
            self.AXIS2.append(self.ilgis)

            self.ilgis=self.ilgis+self.d/2
            self.AXIS1.append(self.R/2.0)
            self.AXIS2.append(self.ilgis)

        self.AXIS1.append(0)
        self.AXIS2.append(self.ilgis)

        self.AXIS1.append(0)
        self.AXIS2.append(0)

        self.AXIS1.append(self.R/2.0)
        self.AXIS2.append(0)
        points3D=[]
        points2D=[]
        for x in range(len(self.AXIS1)):
            points3D.append([self.AXIS1[x],self.AXIS2[x]])
            points2D.append([self.AXIS1[x],self.AXIS2[x]])

        
        rodData="rotate_extrude($fn=200) polygon( points="+str(points3D)+" );\n";
        rod=open("rod3D.scad","w")
        rod.write(rodData)
        rod.close()
        call(["openscad","-o","rod3D.stl","rod3D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("rod3D.stl")
        reader.Update()
        self.contour3D=reader.GetOutput()        
        points2D=[]
        for x in range(0,len(self.AXIS1)):
            
            points2D.append([self.AXIS1[x],self.AXIS2[x]])
            #/linear_extrude(height = 1)
        rodData="translate(["+str(-self.RMAX/2.0)+",0,0]) rotate([90,0,90])linear_extrude(height = "+str(self.RMAX)+") union(){polygon( points="+str(points2D)+" );\nrotate([0,180,0]) polygon( points="+str(points2D)+" );\n}";
        rod=open("rod2D.scad","w")
        rod.write(rodData)
        rod.close()
        call(["openscad","-o","rod2D.stl","rod2D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("rod2D.stl")
        reader.Update()
        self.contour2D=reader.GetOutput()
        bb=[self.boxSize,self.boxSize,self.ilgis-self.d/2]
        bbtranslate=[-bb[0]/2,-bb[1]/2,self.d/4]
        boxData="translate("+str(bbtranslate)+") "+ "cube(size = "+str(bb)+", center = false);\n"
        box=open("betonas3D.scad","w")
        box.write(boxData);
        box.close()
        call(["openscad","-o","betonas3D.stl","betonas3D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("betonas3D.stl")
        reader.Update()
        self.betonas3D=reader.GetOutput()        
        
        bb=[self.RMAX,self.boxSize,self.ilgis-self.d/2]
        bbtranslate=[-self.RMAX/2.0,-bb[1]/2,self.d/4]
        boxData="translate("+str(bbtranslate)+") "+ "cube(size = "+str(bb)+", center = false);\n"
        box=open("betonas2D.scad","w")
        box.write(boxData);
        box.close()
        call(["openscad","-o","betonas2D.stl","betonas2D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("betonas2D.stl")
        reader.Update()
        self.betonas3D=reader.GetOutput()
        
        didesnis=self.BOX_RMAX_N*self.RMAX
        bb=[self.RMAX,      self.boxSize+didesnis,    self.ilgis+didesnis]
        bbtranslate=[-self.RMAX/2.0,-bb[1]/2,-didesnis/2.0]
        boxData=" translate("+str(bbtranslate)+") "+ "cube(size = "+str(bb)+", center = false);\n"
        box=open("box2D.scad","w")
        box.write(boxData);
        box.close()
        call(["openscad","-o","box2D.stl","box2D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("box2D.stl")
        reader.Update()
        self.box2D=reader.GetOutput()
        
        didesnis=self.BOX_RMAX_N*self.RMAX
        bb=[self.boxSize+didesnis,      self.boxSize+didesnis,    self.ilgis+didesnis]
        bbtranslate=[-bb[1]/2,-bb[1]/2,-didesnis/2.0]
        boxData="translate("+str(bbtranslate)+") "+ "cube(size = "+str(bb)+", center = false);\n"
        box=open("box3D.scad","w")
        box.write(boxData);
        box.close()
        call(["openscad","-o","box3D.stl","box3D.scad"])
        reader=vtk.vtkSTLReader()
        reader.SetFileName("box3D.stl")
        reader.Update()
        self.box3D=reader.GetOutput()
        
   
    def STLMESH(self,filename="sample.stl"):
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

    def GenerateParticles(self):
        (mesh,minPoint,maxPoint)=self.STLMESH("box3D.stl")
        
        mntable = MNTable3D (minPoint,maxPoint,self.RMAX*2.0,1)
    
        packer = InsertGenerator3D (self.RMIN,self.RMAX,1000,1000,1.0e-6,True)
    
        packer.generatePacking(mesh,mntable,0,1)       
        
        mntable.generateBondsTagged(0,1.0e-6,0,1,1)
        mntable.write("tempas.vtu",2)
        reader=vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName("tempas.vtu")
        reader.Update()
        self.particles=reader.GetOutput()
        reader=vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName("tempas.vtu")
        reader.Update()
        rodReader = vtk.vtkSTLReader()
        rodReader.SetFileName("rod3D.stl")
        rodReader.Update()

        betonasReader = vtk.vtkSTLReader()
        betonasReader.SetFileName("betonas3D.stl")
        betonasReader.Update()
        v1=vtk.vtkSelectEnclosedPoints()
        v1.SetInputConnection(reader.GetOutputPort())
        v1.SetSurfaceConnection(rodReader.GetOutputPort())
        v1.CheckSurfaceOn ()
        v1.InsideOutOff ()
        v1.SetTolerance (1.0e-6)
        v1.Update()
        print(v1.GetOutput())
        v1insideArray = v1.GetOutput().GetPointData().GetArray("SelectedPoints");
        particleId=[]
        particleId1=[]
        print(v1insideArray.GetNumberOfTuples())
        for x in range(v1insideArray.GetNumberOfTuples()):
            particleId.append(v1insideArray.GetTuple1(x))
            
        v2=vtk.vtkSelectEnclosedPoints()
        v2.SetInputConnection(reader.GetOutputPort())
        v2.SetSurfaceConnection(betonasReader.GetOutputPort())
        v2.CheckSurfaceOn ()
        v2.InsideOutOff ()
        v2.SetTolerance (1.0e-6)
        v2.Update()
        v2insideArray = v2.GetOutput().GetPointData().GetArray("SelectedPoints");
        
        for x in range(v2insideArray.GetNumberOfTuples()):
            particleId1.append(v2insideArray.GetTuple1(x))
            #particleId.append(v2insideArray.GetTuple1(x))
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
            #val=v1insideArray.GetTuple1(x)
            ppp=data.GetPoint(x)
            points.InsertNextPoint(ppp[0]/1000.0,ppp[1]/1000.0,ppp[2]/1000.0)
            radius.InsertNextTuple1(data.GetPointData().GetArray("radius").GetTuple1(x)/1000.0)
            val=0
            if(particleId[x]==0 and particleId1[x]==0):
                val=0
            if(particleId1[x]!=0):
                val=1
            if(particleId[x]!=0):
                val=2
            particle_id.InsertNextTuple1(val)
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
                bonds.SetTuple1(x,0)
            if((v1==0 and v2==2) or (v1==2 and v2==0)):
                bonds.SetTuple1(x,0)
            if(v1==1 and v2==1):
                bonds.SetTuple1(x,1)
            if(v1==2 and v2==2):
                bonds.SetTuple1(x,2)
            if((v1==1 and v2==2) or (v1==2 and v2==1)):
                bonds.SetTuple1(x,3)
                
                
        poly.SetPoints(points)
        poly.GetPointData().SetScalars(radius)
        poly.GetPointData().AddArray(particle_id)
        poly.GetCellData().SetScalars(bonds)
        poly.SetLines(cells)
        poly.GetPointData().SetActiveScalars("ID")
        th=vtk.vtkThreshold()
        th.SetInputData(poly)
        th.ThresholdByUpper (0.5)    
        th.SetInputArrayToProcess(0,0,0,0,'ID')
        th.Update()
        writer=vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName("final.vtu")
        writer.SetInputConnection(th.GetOutputPort())    
        writer.Write()     
        
        boxSTL=vtk.vtkSTLReader()
        boxSTL.SetFileName("box3D.stl")
        boxSTL.Update()
        boxbounds=boxSTL.GetOutput().GetBounds()
        
        betonasSTL=vtk.vtkSTLReader()
        betonasSTL.SetFileName("betonas3D.stl")
        betonasSTL.Update()
        bounds=betonasSTL.GetOutput().GetBounds()
        print(bounds)
        reader=vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName("final.vtu")
        reader.Update()
        KIEKIS=reader.GetOutput().GetNumberOfPoints()
        NumberOfCells=reader.GetOutput().GetNumberOfCells()
        
        poly=vtk.vtkPolyData()
        pp=vtk.vtkPoints();
        pp.SetNumberOfPoints(KIEKIS)
        
        
        rad_seg=vtk.vtkDoubleArray()
        rad_seg.SetName("UNIQUE_RADIUS")
        rad_seg.SetNumberOfComponents(1)
        rad_seg.SetNumberOfTuples(1)
        rad_seg.SetTuple1(0,self.RMIN/1000.0)
        
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
            mat_id=0
            if(pid==2):
                mat_id=0
            if(pid==1):
                mat_id=1
                mat_id=random.randint(1,self.Max_material_count-1)
                    
            part_material.SetTuple1(x,mat_id)
            if(True):
                if(p[2]<bounds[4]/1000.0):
                    vel.SetTuple3(x,0,0,0)
                    part_fix.SetTuple1(x,1) 
                if(p[2]>bounds[5]/1000.0):
                    vel.SetTuple3(x,0,0,self.v)
                    part_fix.SetTuple1(x,1) 
                    
                    
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
            force_N.SetTuple1(x,random.uniform(self.F_N_RANGES[bond_id][0],self.F_N_RANGES[bond_id][1]))
            force_T.SetTuple1(x,random.uniform(self.F_T_RANGES[bond_id][0],self.F_T_RANGES[bond_id][1]))
    
            
            
            
        
        poly.SetPoints(pp)
        poly.GetPointData().SetScalars(rad)
        poly.GetPointData().AddArray(vel)
        poly.GetPointData().AddArray(part_type)
        poly.GetPointData().AddArray(part_fix)
        poly.GetPointData().AddArray(part_material)
        poly.GetFieldData().AddArray(rad_seg)
        
        poly.SetLines(cellsLines)
        poly.GetCellData().SetScalars(state)
        poly.GetCellData().AddArray(force_N)
        poly.GetCellData().AddArray(force_T)
        #
        tran=vtk.vtkTransform()
        tran.Scale(0.001,0.001,0.001)
        tranf=vtk.vtkTransformFilter()
        tranf.SetTransform(tran)
        tranf.SetInputConnection(boxSTL.GetOutputPort())
        tranf.Update()
        bbbbbb=tranf.GetOutput().GetBounds()
        tran=vtk.vtkTransform()
        tran.Translate(-bbbbbb[0],-bbbbbb[2],-bbbbbb[4])
        tranf1=vtk.vtkTransformFilter()
        tranf1.SetTransform(tran)
        tranf1.SetInputConnection(tranf.GetOutputPort())
        tranf1.Update()
        aa=vtk.vtkDataSetWriter()
        aa.SetFileName("mesh.vtk")
        aa.SetInputConnection(tranf1.GetOutputPort())
        aa.Write()
        tranf1=vtk.vtkTransformFilter()
        tranf1.SetTransform(tran)
        tranf1.SetInputData(poly)
        tranf1.Update()
        writer=vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("input.vtp")
        #writer.SetInputData(poly)
        writer.SetInputData(tranf1.GetOutput())
        writer.Write()       
        
        #minPoint,maxPoint
        tran=vtk.vtkTransform()
        tran.Scale(0.001,0.001,0.001)
        tran.Translate(-boxbounds[0],-boxbounds[2],-boxbounds[4])
        tranf=vtk.vtkTransformFilter()
        tranf.SetTransform(tran)
        tranf.SetInputConnection(betonasReader.GetOutputPort())
        tranf.Update()
        www=vtk.vtkDataSetWriter()
        www.SetInputConnection(tranf.GetOutputPort())
        www.SetFileName("betonas3D.vtk")
        www.Write()

        tranf=vtk.vtkTransformFilter()
        tranf.SetTransform(tran)
        tranf.SetInputConnection(rodReader.GetOutputPort())
        tranf.Update()
        www=vtk.vtkDataSetWriter()
        www.SetInputConnection(tranf.GetOutputPort())
        www.SetFileName("rod3D.vtk")
        www.Write()
        
        E=np.linspace(self.E-self.E*self.MATERIAL_PASIMETYMAS_PROC, self.E+self.E*self.MATERIAL_PASIMETYMAS_PROC, num=self.Max_material_count)
        MIU=np.linspace(self.MIU-self.MIU*self.MATERIAL_PASIMETYMAS_PROC, self.MIU+self.MIU*self.MATERIAL_PASIMETYMAS_PROC, num=self.Max_material_count)
        RO=np.linspace(self.RO-self.RO*self.MATERIAL_PASIMETYMAS_PROC, self.RO+self.RO*self.MATERIAL_PASIMETYMAS_PROC, num=self.Max_material_count)
        E[0]=self.E_STEEL
        MIU[0]=self.MIU_STEEL
        RO[0]=self.RO_STEEL
        eilute='"COUNT":'+str(self.Max_material_count)+',\n';
        eilute=eilute+'"E":'+str(E)+",\n"
        eilute=eilute+'"MIU":'+str(MIU)+",\n"
        eilute=eilute+'"RO":'+str(RO)+",\n"
        eilute=eilute+"\n\n"
        
        kk=0
        eil='"CONTACT_LIST":['
        mmm='"MIU_COULOMB":['
        for x in range(self.Max_material_count):
            kk=kk+1
            eil=eil+"["+str(x)+","+str(x)+"],"
            mmm=mmm+str(0.6)+","
        for x in range(self.Max_material_count):
            for y in range(self.Max_material_count):
                if(x!=y and x<y):
                    kk=kk+1
                    eil=eil+"["+str(x)+","+str(y)+"],"  
                    mmm=mmm+str(0.6)+","
    
        eil=eil[:-1]+"],\n"
        mmm=mmm[:-1]+"]\n"
        eilute=eilute+'"COUNT:'+str(kk)+'",\n';
        eilute=eilute+eil+"\n"+mmm+"\n"
        
        
        
        
        print eilute
        
    def Pjauti(self,data):
        plane1=vtk.vtkPlane()
        plane1.SetNormal(1,0,0)
        plane1.SetOrigin(0,0,0)
        plane2=vtk.vtkPlane()
        #plane2.SetNormal(-0.5,0.5,0)
        plane2.SetNormal(0,1,0)
        plane2.SetOrigin(0,0,0)        
        coll=vtk.vtkPlaneCollection()
        coll.AddItem(plane1)
        coll.AddItem(plane2)
        clip1=vtk.vtkClipClosedSurface()
        clip1.SetClippingPlanes(coll)
        clip1.SetInputData(data)
        clip1.Update()
        return clip1.GetOutput()
        
    def Ketvirtadalis(self):
        boxas=vtk.vtkSTLReader()
        boxas.SetFileName("box3D.stl")
        boxas.Update()
        betonas=vtk.vtkSTLReader()
        betonas.SetFileName("betonas3D.stl")
        betonas.Update()
        rod=vtk.vtkSTLReader()
        rod.SetFileName("rod3D.stl")
        rod.Update()
        
        
        w=vtk.vtkSTLWriter()
        zzz=self.Pjauti(betonas.GetOutput())
        w.SetInputData(zzz)
        w.SetFileName("betonas3D.stl")
        w.Write()
        bound1=zzz.GetBounds()
        w=vtk.vtkSTLWriter()
        zzz=self.Pjauti(rod.GetOutput())
        w.SetInputData(zzz)
        w.SetFileName("rod3D.stl")
        w.Write()
        bound2=zzz.GetBounds()
        
        bb=[]
        bb.append(bound1[0])
        bb.append(bound1[1])
        bb.append(bound1[2])
        bb.append(bound1[3])
        bb.append(bound1[4])
        bb.append(bound1[5])
        
        if(bb[0]>bound2[0]):
            bb[0]=bound2[0]
            
        if(bb[2]>bound2[2]):
            bb[2]=bound2[2]
            
        if(bb[4]>bound2[4]):
            bb[4]=bound2[4]
            
            
        if(bb[1]<bound2[1]):
            bb[1]=bound2[1]
            
        if(bb[3]<bound2[3]):
            bb[3]=bound2[3]
            
        if(bb[5]<bound2[5]):
            bb[5]=bound2[5]
            
            
            
            
        didesnis=self.BOX_RMAX_N*self.RMAX
        bb[0]=bb[0]-didesnis
        bb[2]=bb[2]-didesnis
        bb[4]=bb[4]-didesnis
        
        bb[1]=bb[1]+didesnis
        bb[3]=bb[3]+didesnis
        bb[5]=bb[5]+didesnis
            
        box=vtk.vtkCubeSource()
        box.SetBounds(bb)
        box.Update()
        tt=vtk.vtkTriangleFilter()
        tt.SetInputData(box.GetOutput())
        tt.Update()
        
        
        w=vtk.vtkSTLWriter()
        w.SetInputData(tt.GetOutput())
        w.SetFileName("box3D.stl")
        w.Write()
                
        
        
        #    Max_material_count=5
    #MATERIAL_PASIMETYMAS_PROC=10
    #E=5
    #MIU=5
    #RO=1000
        
        
        
        #"COUNT":3,
        #"E":[1.0E+8,3.563E+8,2.03E+8],
        #"G":[0,0,0],
        #"MIU":[0.3,0.15,0.3], 
        #"RO":[1000,2700,8050],



armatura=Strypas()
armatura.generateStrypas()
#armatura.Ketvirtadalis()
armatura.GenerateParticles()

