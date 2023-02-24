import meshio
from dolfin import Mesh,MeshFunction,MeshEditor,Point,Measure,assemble
from dolfin import XDMFFile, HDF5File, MPI
# from dolfin_adjoint import *
import os
import numpy as np

def gmsh2vtk_hex(P):
    return [P[0], P[3], P[2], P[1], P[4], P[5], P[6], P[7]]
	
def vtk2dolfin_hex(P):
    return [P[0], P[1], P[3], P[2], P[4], P[5], P[7], P[6]]

def dolfin2vtk_hex(P):
    return [P[0], P[1], P[3], P[2], P[4], P[5], P[7], P[6]]

def vtk2gmsh_hex(P):
    return [P[0], P[3], P[2], P[1], P[4], P[5], P[6], P[7]]
	
def ufl_simplicial_order(P):
    return np.sort(np.array(P), axis=None).tolist()	

# filename = "STA26_27"
# loadPath = "./InputFiles/22-1215-Wavy_Sweep_v2/"
loadPath = "./InputFiles/22-1215-Wavy_Sweep/"
savePath = "./MeshForVSI/WavySweep06Amp_tet10"
filename = "ShearWavy_6.25MMDisp_Amp_0.6_tet"

def msh2xdmf(loadPath,savePath,filename):

    msh = meshio.read(os.path.join(loadPath, filename + '.inp')) # '.msh'

    #Find mesh size params
    NumVertex = msh.points.shape[0]
    NumCells = msh.cells[0].data.shape[0]

    #Initialize fenics mesh
    mesh = Mesh()

    editor = MeshEditor()
    editor.open(mesh, type="tetrahedron",tdim=3, gdim=3) #"hexahedron"
    editor.init_vertices(NumVertex)
    editor.init_cells(NumCells)
    #HexCellsDataID = msh.cell_data['gmsh:physical'][0]

    for i in range(NumVertex):
        editor.add_vertex(i, Point(msh.points[i,:]))
    for i in range(NumCells):
        try:
            #_cell_gmsh = msh.cells[0].data[i,:]
            #_cell_vtk = gmsh2vtk_hex(_cell_gmsh)
            #_cell_vtk = msh.cells[0].data[i,:]
            #_cell_dolfin = vtk2dolfin_hex(_cell_vtk)			
            _cell_dolfin = ufl_simplicial_order(msh.cells[0].data[i,:])
            editor.add_cell(i,_cell_dolfin)
        except RuntimeError:
            print("Error in cell index %i"%i)
            print(msh.cells[0].data[i,:4])
            raise
    try: 
        print('==============================')
        editor.close()
        print("Editor closed")
    except RuntimeError as err:
        print("Error in closing the editor")
        print(err)
        pass	

    try:
        print('==============================')
        mesh.order() 
        print('Mesh ordered: ', mesh.ordered())
        # Subdomains will be the cells 
        subdomains = MeshFunction('size_t',mesh,mesh.topology().dim())
        #mvc = MeshValueCollection('size_t',mesh,mesh.topology().dim())
        #subdomains.array()[:] = HexCellsDataID
    except RuntimeError:
        print("Error ordering the mesh and generating mesh function")
        raise

    #Write element id as a mesh function
    # write/export cellfunction

    try:
        print('==============================')
        # with XDMFFile(MPI.comm_world,savePath+'/'+filename+'.xdmf') as file:
        with XDMFFile(MPI.comm_world,os.path.join(savePath,filename+'.xdmf')) as file:
            file.write(mesh)
            file.write(subdomains)
        # with HDF5File(MPI.comm_world,savePath+'/'+filename+'_ID.h5','w') as file:
        with HDF5File(MPI.comm_world,os.path.join(savePath,filename+'_ID.h5'),'w') as file:
            file.write(mesh,'/mesh')
            file.write(subdomains,'/meshfunction')
            #file.write(mf,"my_mf")
            #file_IndexData = XDMFFile(MPI.comm_world,filename+'_ID.xdmf');
            #file_IndexData.write(mf)
            #XDMFFile('mesh.xdmf').write(mesh)
        print("Mesh saved")
    except RuntimeError:
        print("Error in saving file")
        raise

    #Test assembly
    ds_custom1 = Measure("ds", domain=mesh, subdomain_data=subdomains, subdomain_id=1)
    ds_custom2 = Measure("ds", domain=mesh, subdomain_data=subdomains, subdomain_id=2)
    ds_custom3 = Measure("ds", domain=mesh, subdomain_data=subdomains, subdomain_id=3)
    print(assemble(1*ds_custom1))
    print(assemble(1*ds_custom2))
    print(assemble(1*ds_custom3))

    # #Now try opening the file
    # mesh2 = Mesh()

    # try:
    #     print('==============================')
    #     with XDMFFile(MPI.comm_world,filename+'.xdmf') as file:
    #         file.read(mesh2)
    #         mvc = MeshValueCollection('size_t',mesh2,mesh2.topology().dim())
    #         file.read(mvc,'f')
    #         mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    #         print(mf.array()[:])
    #         #mf2 = MeshFunction('size_t',mesh2,mesh2.topology().dim())
    #         print('Mesh is readable')	
    # except RuntimeError as err:
    #     print("Error in reading the file")
    #     print(err)
msh2xdmf(loadPath,savePath,filename)