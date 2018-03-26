
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

def mkPacking (mesh=None,minPoint=None,maxPoint=None,Rmin=0.03,Rmax=0.035,maxFailureIterations=2000):

   mntable = MNTable3D (
      minPoint,maxPoint,2.5*Rmax,1
   )
   
   Sphere S(Vector3(px+dx,py+dy,pz+dz),r);
	  bool fit=vol->isIn(S) && ntable->checkInsertable(S,gid);
	  if(fit){
	    S.setTag(tag);
	    ntable->insertChecked(S,gid);

   packer = InsertGenerator3D (
      Rmin,Rmax,maxFailureIterations,1000,1.0e-6,True
   )

   packer.generatePacking(mesh,mntable,0,1)

   #mntable.generateBonds(0,1.e-5,1)

   mntable.write("tempas.vtu",2)
   reader=vtk.vtkXMLUnstructuredGridReader()
   reader.SetFileName("tempas.vtu")
   reader.Update()
   poly=vtk.vtkPolyData()
   poly.SetPoints(reader.GetOutput().GetPoints())
   reader.GetOutput().GetPointData().GetArray("radius").SetName("RADIUS")
   poly.GetPointData().SetScalars(reader.GetOutput().GetPointData().GetArray("RADIUS"))
   writer=vtk.vtkXMLPolyDataWriter()
   writer.SetInputData(poly)
   writer.SetFileName("points.vtp")
   writer.Write()
   
   

if __name__=="__main__":
   print "Gengeo script generator: parameters is stlmesh minradius maxradius, maxFailureIterations=2000"
   (mesh,minPoint,maxPoint)=STLMESH(sys.argv[1])
   mkPacking(mesh=mesh,minPoint=minPoint,maxPoint=maxPoint,Rmin=float(sys.argv[2]),Rmax=float(sys.argv[3]),maxFailureIterations=int(sys.argv[4]))