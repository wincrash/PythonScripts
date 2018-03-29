#!/usr/bin/python
from subprocess import call
import math
import vtk

R=30
h=1.3
a=math.sqrt(3)*h
b=2
c=math.sqrt(3)*h
d=8

COUNT=2
boxSize=30
CalcBoxSizeProc=5






from gengeo import *
import vtk
import sys

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

def mkPacking (mesh=None,minPoint=None,maxPoint=None,Rmin=0.03,Rmax=0.035,maxFailureIterations=2000,taskai=[],outputas="output.vtp"):
   mntable = MNTable3D (
      minPoint,maxPoint,2.5*Rmax,1
   )
   yra=0
   for x in range(len(taskai)):
       S=Sphere(Vector3(float(taskai[x][0]),float(taskai[x][1]),float(taskai[x][2])),float(taskai[x][3]))
       mntable.insert(S,0)
       yra=1


   packer = InsertGenerator3D (
      Rmin,Rmax,maxFailureIterations,1000,1.0e-6,True
   )
   
   packer.generatePacking(mesh,mntable,0,1)

   #mntable.generateBonds(0,1.e-5,1)
   #mntable.generateBonds(0,0.00001,0);
   if(yra==1):
       mntable.generateBondsTagged(0,1.0e-6,0,0,0)
       mntable.generateBondsTagged(0,1.0e-6,1,1,1)
       mntable.generateBondsTagged(2,10,2,0,1)
   mntable.write("tempas.vtu",2)
   reader=vtk.vtkXMLUnstructuredGridReader()
   reader.SetFileName("tempas.vtu")
   reader.Update()
   poly=vtk.vtkPolyData()
   poly.SetPoints(reader.GetOutput().GetPoints())
   reader.GetOutput().GetPointData().GetArray("radius").SetName("RADIUS")
   poly.GetPointData().SetScalars(reader.GetOutput().GetPointData().GetArray("RADIUS"))
   poly.GetPointData().AddArray(reader.GetOutput().GetPointData().GetArray("particleTag"))
   cells=vtk.vtkCellArray()
   for x in range(poly.GetNumberOfPoints()):
       cells.InsertNextCell(1);
       cells.InsertCellPoint(x);
   poly.SetVerts(cells)
   writer=vtk.vtkXMLPolyDataWriter()
   writer.SetInputData(poly)
   writer.SetFileName(outputas)
   writer.Write()



#def mkPacking (mesh=None,minPoint=None,maxPoint=None,Rmin=0.03,Rmax=0.035,maxFailureIterations=2000,taskai,outputas):

def GenerateParticles(stlRod,stlBeton,outputas,ROD_R_MIN,ROD_R_MAX,BETON_R_MIN,BETON_R_MAX,ITER):    
    points=[]
    (mesh,minPoint,maxPoint)=STLMESH(stlRod)
    mkPacking(mesh=mesh,minPoint=minPoint,maxPoint=maxPoint,Rmin=ROD_R_MIN,Rmax=ROD_R_MAX,maxFailureIterations=ITER,taskai=points,outputas="rod.vtp")
    reader=vtk.vtkXMLPolyDataReader()
    reader.SetFileName("rod.vtp")
    reader.Update()
    kiekis=reader.GetOutput().GetNumberOfPoints()
    print("kiekis",kiekis)
    for x in range(kiekis):
        p=reader.GetOutput().GetPoints().GetPoint(x)
        r=reader.GetOutput().GetPointData().GetArray("RADIUS").GetTuple1(x)
        points.append([p[0],p[1],p[2],r])
    
    (mesh,minPoint,maxPoint)=STLMESH(stlBeton)
    mkPacking(mesh=mesh,minPoint=minPoint,maxPoint=maxPoint,Rmin=BETON_R_MIN,Rmax=BETON_R_MAX,maxFailureIterations=ITER,taskai=points,outputas=outputas)










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




def create3D():
    
    points=[]
    
    ilgis=0
    for x in range(COUNT):
    	ilgis,points=GetForma(ilgis,points)
    
    points.append([0,ilgis])
    points.append([0,0])
    points.append([R/2.0,0])
    rodData="translate([0,"+str(R/4)+","+str(R/4)+"])rotate ([-90,0,0])  rotate_extrude($fn=200) polygon( points="+str(points)+" );\n";
    rod=open("rod.scad","w")
    rod.write(rodData)
    rod.close()
    
    cubesize=[boxSize,ilgis-d/2,boxSize]
    translate=[-cubesize[0]/2,d/4,-cubesize[1]/2]
    boxData="translate("+str(translate)+") "+ "cube(size = "+str(cubesize)+", center = false);\n"
    box=open("betonas.scad","w")
    #box.write("difference (){"+boxData+"\n"+rodData+"}\n");
    box.write(boxData);
    box.close()
        
    
    Calccubesize=[boxSize+CalcBoxSizeProc,boxSize+CalcBoxSizeProc,ilgis+CalcBoxSizeProc]
    Calctranslate=[-cubesize[0]/2-CalcBoxSizeProc/2,-cubesize[1]/2-CalcBoxSizeProc/2,-CalcBoxSizeProc/2]
    CalcboxData="translate("+str(Calctranslate)+") "+ "cube(size = "+str(Calccubesize)+", center = false);\n"
    Calcbox=open("box.scad","w")
    Calcbox.write(CalcboxData)
    Calcbox.close()
    
    
    
    call(["openscad","-o","box.stl","box.scad"])
    call(["openscad","-o","betonas.stl","betonas.scad"])
    call(["openscad","-o","rod.stl","rod.scad"])
    
    
 
def create2D(aukstis):
    points=[]
    
    ilgis=0
    for x in range(COUNT):
    	ilgis,points=GetForma(ilgis,points)
    
    points.append([0,ilgis])
    points.append([0,0])
    points.append([R/2.0,0])
    rodData="linear_extrude(height = "+str(aukstis)+", center = true, convexity = 10, twist = -fanrot, slices = 20, scale = 1.0) projection(cut = true) rotate([90,0,0]) rotate_extrude($fn=200) polygon( points="+str(points)+" );\n";
    rod=open("rod.scad","w")
    rod.write(rodData)
    rod.close()
    
    cubesize=[boxSize,aukstis,ilgis-d/2]
    translate=[-cubesize[0]/2,-cubesize[1]/2,d/4]
    boxData="rotate ([90,0,0])  translate("+str(translate)+") "+ "cube(size = "+str(cubesize)+", center = false);\n"
    box=open("betonas.scad","w")
    #box.write("difference (){"+boxData+"\n"+rodData+"}\n");
    box.write(boxData);
    box.close()
        
    
    Calccubesize=[boxSize+CalcBoxSizeProc,aukstis+CalcBoxSizeProc,ilgis+CalcBoxSizeProc]
    Calctranslate=[-cubesize[0]/2-CalcBoxSizeProc/2,-cubesize[1]/2-CalcBoxSizeProc/2,-CalcBoxSizeProc/2]
    CalcboxData="rotate ([90,0,0]) translate("+str(Calctranslate)+") "+ "cube(size = "+str(Calccubesize)+", center = false);\n"
    Calcbox=open("box.scad","w")
    Calcbox.write(CalcboxData)
    Calcbox.close()
    
    
    
    call(["openscad","-o","box.stl","box.scad"])
    call(["openscad","-o","betonas.stl","betonas.scad"])
    call(["openscad","-o","rod.stl","rod.scad"])
    
    
    
       
#create2D(0.8)
#DIMENSION=2


create3D()
DIMENSION=3






def translateAndRotate(data,bounds,scal,outputfilename,DIMENSION):
    
    if(DIMENSION==3):
        tran=vtk.vtkTransform()
        tran.RotateX(-90) 
        tfilter=vtk.vtkTransformFilter()
        tfilter.SetInputData(data)
        tfilter.SetTransform(tran)
        tfilter.Update()
        tran1=vtk.vtkTransform()
        tran1.Translate(-bounds[0],-bounds[2],-bounds[4])   
        tfilter1=vtk.vtkTransformFilter()
        tfilter1.SetInputData(tfilter.GetOutput())
        tfilter1.SetTransform(tran1)
        tfilter1.Update()
        
    
        
        aa=vtk.vtkSTLWriter()
        aa.SetFileName(outputfilename)
        aa.SetInputConnection(tfilter1.GetOutputPort())
        aa.Write()








def translateAndScale(data,bounds,scal,outputfilename,DIMENSION):
    
    
    tran=vtk.vtkTransform()
    tran.Translate(-bounds[0],-bounds[2],-bounds[4])    
    tfilter=vtk.vtkTransformFilter()
    tfilter.SetInputData(data)
    tfilter.SetTransform(tran)
    tfilter.Update()
    tran1=vtk.vtkTransform()
    tran1.Scale(scal,scal,scal)  
    tfilter1=vtk.vtkTransformFilter()
    tfilter1.SetInputData(tfilter.GetOutput())
    tfilter1.SetTransform(tran1)
    tfilter1.Update()
    

    
    aa=vtk.vtkSTLWriter()
    aa.SetFileName(outputfilename)
    aa.SetInputConnection(tfilter1.GetOutputPort())
    aa.Write()
    

boxas=vtk.vtkSTLReader()
boxas.SetFileName("box.stl")
boxas.Update()
bounds=boxas.GetOutput().GetBounds()

rodas=vtk.vtkSTLReader()
rodas.SetFileName("rod.stl")
rodas.Update() 
beton=vtk.vtkSTLReader()
beton.SetFileName("betonas.stl")
beton.Update() 
translateAndScale(boxas.GetOutput(),bounds,0.001,"box.stl",DIMENSION)
translateAndScale(rodas.GetOutput(),bounds,0.001,"rod.stl",DIMENSION)
translateAndScale(beton.GetOutput(),bounds,0.001,"betonas.stl",DIMENSION)
#boxas=vtk.vtkSTLReader()
#boxas.SetFileName("box.stl")
#boxas.Update()
#bounds=boxas.GetOutput().GetBounds()
#rodas=vtk.vtkSTLReader()
#rodas.SetFileName("rod.stl")
#rodas.Update() 
#beton=vtk.vtkSTLReader()
#beton.SetFileName("betonas.stl")
#beton.Update() 
#translateAndRotate(boxas.GetOutput(),bounds,0.001,"box.stl",DIMENSION)
#translateAndRotate(rodas.GetOutput(),bounds,0.001,"rod.stl",DIMENSION)
#translateAndRotate(beton.GetOutput(),bounds,0.001,"betonas.stl",DIMENSION)
    
    

#Calccubesize=[boxSize+CalcBoxSizeProc,ilgis+CalcBoxSizeProc,H2d+CalcBoxSizeProc]
#Calctranslate=[-cubesize[0]/2-CalcBoxSizeProc/2,-CalcBoxSizeProc/2,-H2d/2-CalcBoxSizeProc/2]
#CalcboxData="translate("+str(Calctranslate)+") "+ "cube(size = "+str(Calccubesize)+", center = false);\n"
#Calcbox2D=open("box2D.scad","w")
#Calcbox2D.write(CalcboxData)
#Calcbox2D.close()


#rod2d=open("rod2D.scad","w")
#rod2d.write("translate("+str([0,ilgis,0])+")"+" linear_extrude(height ="+ str(H2d)+", center = true, convexity = 10, twist = -fanrot, slices = 20, scale = 1.0) projection(cut = true) rotate([90,0,0]) "+rodData);
#rod2d.close()

#cubesize=[boxSize,ilgis-d/2,H2d]
#translate=[-cubesize[0]/2,d/4,-H2d/2]
#boxData="translate("+str(translate)+") "+ "cube(size = "+str(cubesize)+", center = false);\n"
#box2D=open("betonas2D.scad","w")
#box2D.write("difference (){"+boxData+"\n"+"translate("+str([0,ilgis,0])+")"+" linear_extrude(height ="+ str(H2d)+", center = true, convexity = 10, twist = -fanrot, slices = 20, scale = 1.0) projection(cut = true) rotate([90,0,0]) "+rodData+"}\n");
#box2D.close()
    

#call(["openscad","-o","rod2D.stl","rod2D.scad"])
#call(["openscad","-o","box2D.stl","box2D.scad"])
#call(["openscad","-o","betonas2D.stl","betonas2D.scad"])


#R_MIN=H2d/2.5
#R_MAX=H2d/2.5
#ITER=100

#GenerateParticles("rod2D.stl","rod2D.vtp",R_MIN,R_MAX,ITER)
#GenerateParticles("betonas2D.stl","betonas2D.vtp",R_MIN,R_MAX,ITER)
