#!/usr/bin/python
import sys
import h5py
import vtk
import os

def GetVectorArray(f,step,name):
    dset = f[step][name][:]
    kiekis=dset.shape[0]
    array=vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(3)
    array.SetNumberOfTuples(kiekis)
    for x in range(0,kiekis):
        array.SetTuple3(x,dset[x][0],dset[x][1],dset[x][2])
    return array

def GetScalarArray(f,step,name):
    dset = f[step][name][:]
    kiekis=dset.shape[0]
    array=vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples(kiekis)
    for x in range(0,kiekis):
        array.SetTuple1(x,dset[x])
    return array

def CreateDataSet(filename,stepID):
    f=h5py.File(filename)
    step=f.keys()[int(stepID)]
    outputName=step+".vtp"
    print("Creating file "+outputName)
    dset = f[step]["POSITIONS"][:]
    cube    = vtk.vtkPolyData()
    points  = vtk.vtkPoints()
    polys   = vtk.vtkCellArray()
    kiekis=dset.shape[0]
    points.SetNumberOfPoints(kiekis);
    for x in range(0,kiekis):
        points.SetPoint(x,dset[x][0],dset[x][1],dset[x][2])
        polys.InsertNextCell(1);
        polys.InsertCellPoint(x);    
    if("UNIQUE_RADIUS" in f[step].keys()):
        array=GetScalarArray(f,step,"UNIQUE_RADIUS")
        cube.GetFieldData().SetScalars(array)
    if("RADIUS" in f[step].keys()):
        array=GetScalarArray(f,step,"RADIUS")
        cube.GetPointData().SetScalars(array)
    if("PARTICLE_TYPE" in f[step].keys()):
        array=GetScalarArray(f,step,"PARTICLE_TYPE")
        cube.GetPointData().AddArray(array)
    if("MASS" in f[step].keys()):
        array=GetScalarArray(f,step,"MASS")
        cube.GetPointData().AddArray(array)
    if("VELOCITY" in f[step].keys()):
        array=GetVectorArray(f,step,"VELOCITY")
        cube.GetPointData().SetVectors(array)
    if("FORCE" in f[step].keys()):
        array=GetVectorArray(f,step,"FORCE")
        cube.GetPointData().AddArray(array)
    if("TORQUE" in f[step].keys()):
        array=GetVectorArray(f,step,"TORQUE")
        cube.GetPointData().AddArray(array)
    if("ANGULAR_VELOCITY" in f[step].keys()):
        array=GetVectorArray(f,step,"ANGULAR_VELOCITY")
        cube.GetPointData().AddArray(array)
    if("ANGULAR_ACCELERATION" in f[step].keys()):
        array=GetVectorArray(f,step,"ANGULAR_ACCELERATION")
        cube.GetPointData().AddArray(array)
    if("ACCELERATION" in f[step].keys()):
        array=GetVectorArray(f,step,"ACCELERATION")
        cube.GetPointData().AddArray(array)
    cube.SetPoints(points)
    cube.SetVerts(polys)    
    writer=vtk.vtkXMLPolyDataWriter()    
    writer.SetInputData(cube)
    writer.SetFileName(outputName)
    writer.Write()
    f.close()
    return outputName




def CreateDataSetMesh(filename,stepID):
    f=h5py.File(filename)
    step=f.keys()[int(stepID)]
    outputName="MESH_"+step+".vtp"
    print("Creating file "+outputName)
    dset = f[step]["BOUNDARY_POINTS"][:]
    tria=f[step]["BOUNDARY_IDS"][:]
    cube    = vtk.vtkPolyData()
    points  = vtk.vtkPoints()
    polys   = vtk.vtkCellArray()
    kiekis=dset.shape[0]
    points.SetNumberOfPoints(kiekis);
    for x in range(0,kiekis):
        points.SetPoint(x,dset[x][0],dset[x][1],dset[x][2])
        
    for x in range(0,tria.shape[0]):
        polys.InsertNextCell(3);
        polys.InsertCellPoint(int(tria[x][0]));
        polys.InsertCellPoint(int(tria[x][1]));
        polys.InsertCellPoint(int(tria[x][2]));
        
    if("BOUNDARY_VELOCITY" in f[step].keys()):
        array=GetVectorArray(f,step,"BOUNDARY_VELOCITY")
        cube.GetPointData().AddArray(array)
    if("BOUNDARY_FORCE" in f[step].keys()):
        array=GetVectorArray(f,step,"BOUNDARY_FORCE")
        cube.GetCellData().AddArray(array)
        
    cube.SetPoints(points)
    cube.SetPolys(polys)    
    writer=vtk.vtkXMLPolyDataWriter()    
    writer.SetInputData(cube)
    writer.SetFileName(outputName)
    writer.Write()
    f.close()
    return outputName


def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
def PrintSteps(filename):
    f=h5py.File(filename)
    for x in range(0,len(f.keys())):
        print(x,str(f.keys()[x]))
    step=f.keys()[int(0)]
    print(f[step].keys())
    f.close()
        

if len(sys.argv) >= 3:
    print(sys.argv)
    if( "list" in str(sys.argv[2])):
        PrintSteps(sys.argv[1])
    if( "all" in str(sys.argv[2])):
        F = open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+".pvd",'w') 
        F.write("<VTKFile type=\"Collection\" version=\"0.1\" " +  "byte_order=\"LittleEndian\">\n"+  "<Collection>\n")
        f=h5py.File(sys.argv[1])
        kiekis=len(f.keys())
        f.close()
        for x in range(kiekis):
            filename=CreateDataSet(sys.argv[1],x)
            F.write("<DataSet part=\"0\"  timestep=\""+str(x)+ "\" file=\"" +filename+ "\"/>\n")
        F.write("\n</Collection>\n</VTKFile>\n")
        F.close()
    if( "single" in str(sys.argv[2])):
        F = open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+".pvd",'w') 
        F.write("<VTKFile type=\"Collection\" version=\"0.1\" " +  "byte_order=\"LittleEndian\">\n"+  "<Collection>\n")
        filename=CreateDataSet(sys.argv[1],int(sys.argv[3]))
        F.write("<DataSet part=\"0\"  timestep=\""+str(0)+ "\" file=\"" +filename+ "\"/>\n")
        F.write("\n</Collection>\n</VTKFile>\n")
        F.close()
    if( "range" in str(sys.argv[2])):
        F = open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+".pvd",'w') 
        F.write("<VTKFile type=\"Collection\" version=\"0.1\" " +  "byte_order=\"LittleEndian\">\n"+  "<Collection>\n")
        for x in range(int(sys.argv[3]),int(sys.argv[4])+1):
            filename=CreateDataSet(sys.argv[1],x)
            F.write("<DataSet part=\"0\"  timestep=\""+str(x)+ "\" file=\"" +filename+ "\"/>\n")
        F.write("\n</Collection>\n</VTKFile>\n")
        F.close()
    if( "mesh" in str(sys.argv[2])):
        F = open(os.path.splitext(os.path.basename(sys.argv[1]))[0]+"_MESH.pvd",'w')         
        F.write("<VTKFile type=\"Collection\" version=\"0.1\" " +  "byte_order=\"LittleEndian\">\n"+  "<Collection>\n")
        f=h5py.File(sys.argv[1])
        kiekis=len(f.keys())
        f.close()
        for x in range(kiekis):
            filename=CreateDataSetMesh(sys.argv[1],x)
            F.write("<DataSet part=\"0\"  timestep=\""+str(x)+ "\" file=\"" +filename+ "\"/>\n")
        F.write("\n</Collection>\n</VTKFile>\n")
        F.close()
        
else:
    print ("Galimi du parametrai:")
    print ("1 Parametras = HDF5 failas")
    print ("2 Parameters = galimi variatnai: ")
    print ("list - atspaudins stepu pavadinimas")
    print ("all - sukurs visus VTK failus")
    print ("mesh - sukurs mesh failus")
    print ("single - ir papildomai grupes nr")
    print ("range - ir papildomai grupes nr nuo iki")
    




