########################
# KUHN SIMPLEX CLASS
# Mesh generation utilities to create domain for embedded flow simulations
# - Collected several FRG functions/class into one file
# - Ref: bitbucket.org/frg
########################
class KuhnSimplex:
    def __init__(self,x,y,z, boundaryNames):
        '''
        :param x: float array
        :param y: float array
        :param z: float array
        :param boundaryNames:  [bottom z = zmin, top  z = zmax, left x = xmin,
                                      right x = xmax, front y = ymin, back y = ymax]
        '''
        self.x,  self.y,  self.z = x,y,z
        self.nx, self.ny, self.nz = len(x), len(y), len(z)
        self.nodes = self.create_nodes()
        self.eles = self.create_tet()
        self.boundaries = self.create_boundaries()
        self.boundaryNames = boundaryNames


    def _node_id(self,ix,iy,iz):
        nx, ny, nz = self.nx, self.ny, self.nz
        return ix + iy*nx + iz*nx*ny


    def create_nodes(self):

        import numpy as np

        x,y,z = self.x, self.y, self.z
        nx, ny, nz = self.nx, self.ny, self.nz

        nodes = np.empty((nx*ny*nz, 3))

        xx,yy,zz = np.meshgrid(x,y,z)

        nodes[:, 0] = np.reshape(np.reshape(xx, (-1, nz)), (1, -1), order='F')
        nodes[:, 1] = np.reshape(np.reshape(yy, (-1, nz)), (1, -1), order='F')
        nodes[:, 2] = np.reshape(np.reshape(zz, (-1, nz)), (1, -1), order='F')

        return nodes


    def create_tet(self):
        '''
        :param xx:
        :param yy:
        :param zz:
        cut each small cube into 6 tetrahedrons

        nodes = [n0,n1,n2,n3,n4,n5,n6,n7]

        n0: (0,0,0)
        n1: (1,0,0)
        n2: (0,1,0)
        n3: (1,1,0)
        n4: (0,0,1)
        n5: (1,0,1)
        n6: (0,1,1)
        n7: (1,1,1)

        T(1 2 3) n0(0, 0, 0), n1(1, 0, 0), n3(1, 1, 0), n7(1, 1, 1)
        T(1 3 2) n0(0, 0, 0), n1(1, 0, 0), n5(1, 0, 1), n7(1, 1, 1)
        T(2 1 3) n0(0, 0, 0), n2(0, 1, 0), n3(1, 1, 0), n7(1, 1, 1)
        T(3 1 2) n0(0, 0, 0), n2(0, 1, 0), n6(0, 1, 1), n7(1, 1, 1)
        T(3 2 1) n0(0, 0, 0), n4(0, 0, 1), n6(0, 1, 1), n7(1, 1, 1)
        T(2 3 1) n0(0, 0, 0), n4(0, 0, 1), n5(1, 0, 1), n7(1, 1, 1)

        :return:
        '''

        import numpy as np

        nx, ny, nz = self.nx, self.ny, self.nz

        eles = np.empty((6*(nx - 1) * (ny - 1) * (nz - 1), 4), dtype=int)
        for k in range(nz - 1):
            for j in range(ny - 1):
                for i in range(nx - 1):
                    cubeId = i + j*(nx-1) + k*(nx-1)*(ny-1)
                    nn = [self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i, j + 1, k), self._node_id(i + 1, j + 1, k),
                          self._node_id(i, j, k + 1), self._node_id(i + 1, j, k + 1), self._node_id(i, j + 1, k + 1), self._node_id(i + 1, j + 1, k + 1)]
                    eles[6 * cubeId:6 * (cubeId + 1), :] = [[nn[0], nn[1], nn[3], nn[7]],
                                                            [nn[0], nn[1], nn[5], nn[7]],
                                                            [nn[0], nn[2], nn[3], nn[7]],
                                                            [nn[0], nn[2], nn[6], nn[7]],
                                                            [nn[0], nn[4], nn[6], nn[7]],
                                                            [nn[0], nn[4], nn[5], nn[7]]]
        return eles


    def create_boundaries(self):
        nx, ny, nz = self.nx, self.ny, self.nz
        tri = [[],[],[],[],[],[]]
        # bottom z = zmin
        for j in range(ny - 1):
            for i in range(nx - 1):
                k = 0
                tri[0].append([self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i + 1, j + 1, k)])
                tri[0].append([self._node_id(i, j, k), self._node_id(i + 1, j + 1, k), self._node_id(i, j + 1, k)])


        # top  z = zmax
        for j in range(ny - 1):
            for i in range(nx - 1):
                k = nz - 1
                tri[1].append([self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i + 1, j + 1, k)])
                tri[1].append([self._node_id(i, j, k), self._node_id(i + 1, j + 1, k), self._node_id(i, j + 1, k)])

        # left x = xmin
        for k in range(nz - 1):
            for j in range(ny - 1):
                i = 0
                tri[2].append([self._node_id(i, j, k), self._node_id(i , j + 1, k),     self._node_id(i, j + 1, k + 1)])
                tri[2].append([self._node_id(i, j, k), self._node_id(i , j + 1, k + 1), self._node_id(i, j , k + 1)])
        # right x = xmax
        for k in range(nz - 1):
            for j in range(ny - 1):
                i = nx - 1
                tri[3].append([self._node_id(i, j, k), self._node_id(i, j + 1, k), self._node_id(i, j + 1, k + 1)])
                tri[3].append([self._node_id(i, j, k), self._node_id(i, j + 1, k + 1), self._node_id(i, j, k + 1)])
        # front y = ymin
        for k in range(nz - 1):
            for i in range(nx - 1):
                j = 0
                tri[4].append([self._node_id(i, j, k), self._node_id(i, j , k + 1), self._node_id(i + 1, j, k + 1)])
                tri[4].append([self._node_id(i, j, k), self._node_id(i + 1, j, k + 1), self._node_id(i + 1, j, k)])
        # back y = ymax
        for k in range(nz - 1):
            for i in range(nx - 1):
                j = ny - 1
                tri[5].append([self._node_id(i, j, k), self._node_id(i, j, k + 1), self._node_id(i + 1, j, k + 1)])
                tri[5].append([self._node_id(i, j, k), self._node_id(i + 1, j, k + 1), self._node_id(i + 1, j, k)])

        return tri


    def write_topfile(self, outputfile = 'domain.top', volFunc = lambda x: 0):
        nodes = self.nodes
        eles = self.eles
        boundaries = self.boundaries
        boundaryNames = self.boundaryNames
        self.write_tet(nodes, eles, boundaryNames, boundaries, outputfile, volFunc)


    def plot_mesh(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        eles  = self.eles
        nodes = self.nodes
        for e in eles:
            for i in range(4):
                for j in range(i+1,4):
                    x = [nodes[e[i],0],nodes[e[j],0]]
                    y = [nodes[e[i],1],nodes[e[j],1]]
                    z = [nodes[e[i],2],nodes[e[j],2]]
                    ax.plot(x,y,z)
        plt.show()



    def symmetry(self, xx, type):
        '''
        :param xx: array {x0, x1, x2}
        :param type: string, 'left' or 'right'
        :return: mirror the array to left or right
                if type = 'left' xx = {2x0 - x2, 2x0 - x1, x0, x1, x2}
                if type = 'right' xx = {x0, x1, x2, 2x2 - x1, 2x2 - x0}
        '''

        n = len(xx)

        if(type == 'left'):
            xx_temp = []
            xc = xx[0]
            for i in range(n - 1, 0, -1):
                xx_temp.append(2*xc - xx[i])
            xx = xx_temp + xx
        if(type == 'right'):
            xc = xx[-1]
            for i in range(n - 1):
                xx.append(2*xc - xx[n - 1 - i - 1])
        return xx



    def write_tet (self, nodes, elems, boundaryNames, boundaries, mshfile = 'domain.top', volFunc = lambda x: 0, flagPrintInfo = True):
        '''
        This function writes top file, mshfile, node coordinates are  in the list nodes,
        and element node numbers are in elems(the node number is 0-based index instead of top file'
        1-based index)
        :param mshfile: top file name
        :param nodes: a list of node coordinates
        :param elems: a list of elems node number
        :param boundaries: a list of several lists, each sublist is a list of boundary triangle node numbers
        :param volFunc: function returns 0, 1, 2, 3... input point, return the point's volume label
        '''

        import numpy as np

        if flagPrintInfo == True:
            print ('>> >> Saving mesh to TOP file: "{}"'.format (mshfile))

        with open(mshfile, 'w') as file:
            nNodes, nElems = len(nodes), len(elems)

            file.write('Nodes FluidNodes\n')
            for nN in range(nNodes):
                file.write("{:d}    {:.12f}  {:.12f}  {:.12f}\n".format (nN + 1, nodes[nN][0],nodes[nN][1],nodes[nN][2]))

            iElems, volLabel = 0, 0
            while(iElems < nElems):
                file.write("Elements FluidMesh_{:d} using FluidNodes\n".format (volLabel))
                for nE in range(nElems):
                    xc = (np.array(nodes[elems[nE][0]]) + np.array(nodes[elems[nE][1]]) + np.array(nodes[elems[nE][2]]) + np.array(nodes[elems[nE][3]]))/4.0
                    if volFunc(xc) == volLabel:
                        file.write("{:d}    {:d}    {:d}    {:d}    {:d}    {:d}\n".format (nE + 1, 5, elems[nE][0] + 1, elems[nE][1] + 1, elems[nE][2] + 1, elems[nE][3] + 1))
                        iElems = iElems + 1
                volLabel = volLabel + 1
            nBounds = len(boundaries)
            for i in range(nBounds):
                if(boundaries == 'None'):
                    continue
                file.write("Elements {:s} using FluidNodes\n".format (boundaryNames[i] + "_" + str(i+volLabel)))
                boundary = boundaries[i]
                nTris = len(boundary)
                for nT in range(nTris):
                    nE += 1
                    file.write("{:d}    {:d}    {:d}    {:d}    {:d}\n".format (nE, 4, boundary[nT][0] + 1, boundary[nT][1] + 1, boundary[nT][2] + 1))

        if flagPrintInfo == True:
            print (">> >> Wrote {} nodes and {} elements to TOP file".format (len(nodes), len(elems)))



########################
# MESHGEN UTIL FUNCTIONS
# Mesh generation utilities
# - createRectDomain: Create a simple rectangular domain for embedded flow simulations using Kuhn Simplex Class
# - rotateMesh: Rotate a mesh specified angles around an origin point
# - runXP2EXO: Use shell commands to run XP2EXO on mesh files - Note requires FRG software to be installed
# - runFEM2TOP: Use shell commands to run FEM2TOP on mesh files - Note requires FRG software to be installed
########################
def createRectDomain(xMin, xMax, yMin, yMax, zMin, zMax,
                     dx, dy, dz, 
                     boundaryNames=['InletFixed_Inlet_Surface_zmin','InletFixed_Inlet_Surface_zmax','InletFixed_Inlet_Surface_xmin', 'InletFixed_Inlet_Surface_xmax','InletFixed_Inlet_Surface_ymin','InletFixed_Inlet_Surface_ymax'], 
                     outputFile="domain.top",
                     flagPrintInfo=True):
    """
    Create a simple rectangular domain for embedded flow simulations
    Inputs:
        xMin...zMax:    boundary locations for rectangular volume
        dx...dz:        desired grid spacing (will approximate based on closest number of points to create this)
        boundaryNames:  names of each surface in mesh top file - note Kuhn simplex function above writes surfaces in order min/max of z, x, y
    """
    import numpy as np
    
    # Grid parameters
    xPoints = int ((xMax - xMin) / dx) + 1
    yPoints = int ((yMax - yMin) / dy) + 1
    zPoints = int ((zMax - zMin) / dz) + 1

    # Create Grid
    x = np.linspace (xMin, xMax, xPoints)
    y = np.linspace (yMin, yMax, yPoints)
    z = np.linspace (zMin, zMax, zPoints)

    # Create Kuhn Simplex ordered domain
    simpleKuhnSimplex = KuhnSimplex (x, y, z, boundaryNames)
    if flagPrintInfo == True:
        print (">> >> Domain contains {} x {} x {} nodes".format (xPoints, yPoints, zPoints))
        print (">> >> Writing to top file")
        
    # Write output to Xpost (TOP) file
    simpleKuhnSimplex.write_topfile (outputFile)
    
    # Return nodes and elements
    return simpleKuhnSimplex.nodes, simpleKuhnSimplex.eles



def rotateMesh (meshData, rotAngle, rotOrigin=[0,0,0]):
    """ Rotate mesh points around x, y, z axis in sequence
        Inputs:
            meshData:   numpy 2D array (i.e. list of vectors for each point)
            rotAngle:   list of rotation angles in degrees
            rotOrigin:  list defining origin for center of rotation
    """
    import numpy as np

    # Create rotation matrix based on multiplication of x, y, z rotation matrices
    rotX = np.radians(rotAngle[0])
    rotY = np.radians(rotAngle[1])
    rotZ = np.radians(rotAngle[2])
    R = np.zeros ((3, 3))
    R[0, 0] = np.cos(rotY) * np.cos(rotZ)
    R[0, 1] = -np.cos(rotY) * np.sin(rotZ)
    R[0, 2] = np.sin(rotY)
    R[1, 0] = np.cos(rotX) * np.sin(rotZ) + np.cos(rotZ) * np.sin(rotX) * np.sin(rotY)
    R[1, 1] = np.cos(rotX) * np.cos(rotZ) - np.sin(rotX) * np.sin(rotY) * np.sin(rotZ)
    R[1, 2] = -np.cos(rotY) * np.sin(rotX)
    R[2, 0] = np.sin(rotX) * np.sin(rotZ) - np.cos(rotX) * np.cos(rotZ) * np.sin(rotY)
    R[2, 1] = np.cos(rotZ) * np.sin(rotX) + np.cos(rotX) * np.sin(rotY) * np.sin(rotZ)
    R[2, 2] = np.cos(rotX) * np.cos(rotY)

    # Apply rotations to each node
    meshDataTmp = np.copy (meshData)                                                   # create copy of array so we don't modify input
    for nodeTmp in range (len (meshDataTmp)):
        meshDataTmp[nodeTmp, :] -= np.array (rotOrigin)                                # translate points based on center of rotation
        meshDataTmp[nodeTmp, :] = np.matmul (R, meshDataTmp[nodeTmp, :])               # apply rotation matrix
        meshDataTmp[nodeTmp, :] += np.array (rotOrigin)                                # translate back to original origin

    # Return result
    return meshDataTmp



def runXP2EXO (filename):
    """ Run XP2EXO to convert file from AERO-F to Paraview format
        - Assumes user has compiled and set up XP2EXO on the machine, and binary is accessible via XP2EXO environment variable
        - Will source config.sh file in meshGen directory - use this to source appropriate files and setup XP2EXO variable
    """
    import subprocess
    
    # Print info
    print ("Running XP2EXO to convert file: {}".format (filename+".top"))
    
    # Create array of string commands
    cmd = ["source ./config.sh",
           "$XP2EXO {} {}".format (filename+".top", filename+".exo")]
           
    # Run commands joined by semi-colons to allow consecutive execution in bash
    try:
        subprocess.run (";".join(cmd), shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
        print ("ERROR running XP2EXO: Make sure FRG suite is installed, and config.sh file exists and has XP2EXO command defined. Exodus file not generated.")
    
    
    
def runFEM2TOP (outputDir, filename):
    """ Run AEROS to convert file from AERO-S format to TOP
        - Assumes user has compiled and set up AERO-S on the machine, and binary is accessible via AEROS environment variable
        - Will source config.sh file in meshGen directory - use this to source appropriate files and setup AEROS variable
    """
    import subprocess
    
    # Print info
    print ("Running AEROS to convert file: {}".format (filename+".fem"))
    
    # Create filename for conversion deck
    filenameConvert = "convert_fem2top.aeros"
    
    # Write AEROS input deck to do conversion
    f = open (outputDir+"/"+filenameConvert, "w")
    f.write ('CONTROL\n')
    f.write ('"{}"\n'.format (filename))
    f.write ('1\n')
    f.write ('"nodeset"\n')
    f.write ('"elemset"\n')
    f.write ('***************************************************************************\n')
    f.write ('INCLUDE "{}"\n'.format (filename+".fem"))
    f.write ('END\n')
    f.close ()
    
    # Create array of string commands
    cmd = ["source ./config.sh",
           "cd {}".format (outputDir),
           "$AEROS -t {}".format (filenameConvert)]
           
    # Run commands joined by semi-colons to allow consecutive execution in bash
    try:
        subprocess.run (";".join(cmd), shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
        print ("ERROR running FEM2TOP: Make sure FRG suite is installed, and config.sh file exists and has AEROS command defined. TOP file not generated.")
