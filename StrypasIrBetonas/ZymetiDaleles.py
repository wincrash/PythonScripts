#!/usr/bin/python
from gengeo import *
import vtk
import random
import GenGeoModule

box="box.stl"
betonas="betonas.stl"
rodas="rod.stl"
RMIN=0.001
RMAX=0.001
ITER=1000
v=0.01
F_N_RANGES=[[0,0],[5,4],[60,50],[20,15]]
F_T_RANGES=F_N_RANGES


def DataSetPrepare():
    
    boxSTL=vtk.vtkSTLReader()
    boxSTL.SetFileName(box)
    boxSTL.Update()
    boxbounds=boxSTL.GetOutput().GetBounds()
    
    betonasSTL=vtk.vtkSTLReader()
    betonasSTL.SetFileName(betonas)
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
        part_material.SetTuple1(x,pid)
        if(False):
            if(p[2]<bounds[4]):
                vel.SetTuple3(x,0,0,0)
                part_fix.SetTuple1(x,1) 
            if(p[2]>bounds[5]):
                vel.SetTuple3(x,0,0,v)
                part_fix.SetTuple1(x,1)         
        else:
            if(p[1]<bounds[2]):
                vel.SetTuple3(x,0,0,0)
                part_fix.SetTuple1(x,1) 
            if(p[1]>bounds[3]):
                vel.SetTuple3(x,0,v,0)
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
        force_N.SetTuple1(x,random.uniform(F_N_RANGES[bond_id][0],F_N_RANGES[bond_id][1]))
        force_T.SetTuple1(x,random.uniform(F_T_RANGES[bond_id][0],F_T_RANGES[bond_id][1]))

        
        
        
    
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
    
    writer=vtk.vtkXMLPolyDataWriter()
    writer.SetFileName("input.vtp")
    #writer.SetInputData(poly)
    writer.SetInputData(poly)
    writer.Write()
    
    aa=vtk.vtkDataSetWriter()
    aa.SetFileName("mesh.vtk")
    aa.SetInputConnection(boxSTL.GetOutputPort())
    aa.Write()



def Sperate():
    reader=vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("tempas.vtu")
    reader.Update()
    rodReader = vtk.vtkSTLReader()
    rodReader.SetFileName(rodas)
    rodReader.Update()
    betonasReader = vtk.vtkSTLReader()
    betonasReader.SetFileName(betonas)
    betonasReader.Update()
    v1=vtk.vtkSelectEnclosedPoints()
    v1.SetInputConnection(reader.GetOutputPort())
    v1.SetSurfaceConnection(rodReader.GetOutputPort())
    v1.CheckSurfaceOn ()
    v1.InsideOutOff ()
    v1.SetTolerance (1.0e-6)
    v1.Update()
    v1insideArray = v1.GetOutput().GetPointData().GetArray("SelectedPoints");
    particleId=[]
    particleId1=[]
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
        points.InsertNextPoint(data.GetPoint(x))
        radius.InsertNextTuple1(data.GetPointData().GetArray("radius").GetTuple1(x))
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

def mkPacking (mesh=None,minPoint=None,maxPoint=None,Rmin=0.03,Rmax=0.035,maxFailureIterations=2000):

   mntable = MNTable3D (
      minPoint,maxPoint,2.5*Rmax,1
   )

   packer = InsertGenerator3D (
      Rmin,Rmax,maxFailureIterations,1000,1.0e-6,True
   )

   packer.generatePacking(mesh,mntable,0,1)
   mntable.generateBondsTagged(0,1.0e-6,0,1,1)
   mntable.write("tempas.vtu",2)
   
   

if __name__=="__main__":
#   print "Gengeo script generator: parameters is stlmesh minradius maxradius, maxFailureIterations=2000"
   #(mesh,minPoint,maxPoint)=STLMESH("box.stl")
   #mkPacking(mesh=mesh,minPoint=minPoint,maxPoint=maxPoint,Rmin=RMIN,Rmax=RMAX,maxFailureIterations=ITER)
   GenGeoModule.GeneratePacking(box,RMIN,RMAX,)
   Sperate()
   DataSetPrepare()