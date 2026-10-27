#####################################################################################
# RUN MESHGEN
#####################################################################################
# Main script to run MeshGen and create meshes for AERO Suite / FSI analysis
# - Provides capability to convert STL files into appropriate format for AERO-F/S
# - Includes ability to use parachute classes to parametrically create parachute meshes
#   as well for use in FSI analysis
# - Users can add additional mesh modules similar to how parachute mesh is implemented
#####################################################################################
 

########################
# INPUTS - Main Settings
########################
runMeshModule = "parachuteDGB"              # None, parachuteDGB, parachuteRingsail, etc.
outputDir = "./output"                      # Output directory
outputEmbeddedPrefix = "meshEmbedded_"      # Prefix for embedded mesh files
outputDomainFilename = "domain.top"         # Filename for domain mesh


########################
# INPUTS - Mesh Definition
########################

# DGB Parachute
domainX                 = [-80.0, 80.0]
domainY                 = [-80.0, 80.0]
domainZ                 = [-20.0, 200.0]
domainSpacing           = 4.0
inputEmbedded = {"Payload" : {
                    "filename" : "input/aspire_mesh_cubit.stl",
                    "bcType" : "StickFixed",                            # StickFixed, StickMoving, SlipFixed, SlipMoving
                    "scaleFactor" : 0.001,                              # Scale factor to convert source file data to [m]
                    "rotation" : [0.0, 180.0, 0.0],                     # Applied in order after scaling in [deg]
                    "offset" : [0.0, 0.0, -6.67677],                    # Applied after rotation in [m]
                    "cg" : [0.0, 0.0, -(6.67677-2.8227)]               # Relative to scaled/rotated/translated object in [m]
                    }
                }
# inputEmbedded = {"Payload" : {
#                     "filename" : "input/msl_mesh_cubit.stl",
#                     "bcType" : "StickFixed",                            # StickFixed, StickMoving, SlipFixed, SlipMoving
#                     "scaleFactor" : 0.001,                              # Scale factor to convert source file data to [m]
#                     "rotation" : [90.0, 0.0, 0.0],                      # Applied in order after scaling in [deg]
#                     "offset" : [0.0, 0.0, -2.9134],                     # Applied after rotation in [m]
#                     "cg" : [0.0, 0.0, -(2.9134-1.278)]                   # Relative to scaled/rotated/translated object in [m]
#                     }
#                 }


# RAD Test
# domainX                 = [-200.0, 200.0]
# domainY                 = [-200.0, 200.0]
# domainZ                 = [-200.0, 200.0]
# domainSpacing           = 10.0
# inputEmbedded = {"RAD_body" : {
#                     "filename" : "input/skyfall/rad_body.stl",
#                     "bcType" : "StickFixed",                            # StickFixed, StickMoving, SlipFixed, SlipMoving
#                     "scaleFactor" : 0.001,                              # Scale factor to convert source file data to [m]
#                     "rotation" : [90.0, 0.0, 180.0],                      # Applied in order after scaling in [deg]
#                     "offset" : [0.0, 0.0, 0.08636],                     # Applied after rotation in [m]
#                     "cg" : [0.0, 0.0, 0.0]                  # Relative to scaled/rotated/translated object in [m]
#                     },
#                  "Nozzle_Exit" : {
#                     "filename" : "input/skyfall/nozzle_exit.stl",
#                     "bcType" : "OutletFixed",                            # StickFixed, StickMoving, SlipFixed, SlipMoving, OutletFixed, OutletMoving
#                     "scaleFactor" : 0.001,                              # Scale factor to convert source file data to [m]
#                     "rotation" : [90.0, 0.0, 180.0],                      # Applied in order after scaling in [deg]
#                     "offset" : [0.0, 0.0, 0.08636],                     # Applied after rotation in [m]
#                     "cg" : [0.0, 0.0, 0.0]                  # Relative to scaled/rotated/translated object in [m]
#                     }
#                 }
# inputEmbedded = {"RAD_body" : {
#                     "filename" : "input/skyfall/rad_nonfiring.stl",
#                     "bcType" : "StickFixed",                            # StickFixed, StickMoving, SlipFixed, SlipMoving
#                     "scaleFactor" : 0.001,                              # Scale factor to convert source file data to [m]
#                     "rotation" : [90.0, 0.0, 180.0],                      # Applied in order after scaling in [deg]
#                     "offset" : [0.0, 0.0, 0.08636],                     # Applied after rotation in [m]
#                     "cg" : [0.0, 0.0, 0.0]                  # Relative to scaled/rotated/translated object in [m]
#                     },
#                 }



########################
# MAIN FUNCTION
########################
def meshGenMain ():
    """
    Main analysis function to create fluid and structure meshes for AERO suite
    Author: Marcus Lobbia
    """
    import copy
    import os
    import numpy as np
    from meshGen_utils import createRectDomain, rotateMesh, runXP2EXO, runFEM2TOP
    import meshio
    
    # Other constants
    dimSurf = 2
    dimVol = 3

    # Print info to screen
    print (">> ============================================================")
    print (">> MESHGEN: Mesh generator for AERO Suite / FSI analysis")
    print (">> by Marcus Lobbia")
    print (">> ============================================================")

    # Create output dir if non-existant, or delete output inside if exists
    if not os.path.exists (outputDir):
        os.makedirs (outputDir)
    else:
        filelist = [f for f in os.listdir (outputDir)]
        for f in filelist:
            os.remove (os.path.join (outputDir, f)) 
    
    
    # FLUID DOMAIN MESH
    
    # Create fluid domain mesh
    print (">> Creating rectangular domain for embedded mesh simulation...")
    domainNodes, domainElements = createRectDomain   (xMin=domainX[0], xMax=domainX[1], 
                                                      yMin=domainY[0], yMax=domainY[1], 
                                                      zMin=domainZ[0], zMax=domainZ[1],
                                                      dx=domainSpacing, dy=domainSpacing, dz=domainSpacing,
                                                      outputFile=outputDir+"/"+outputDomainFilename)
    print (">> >> Writing exodus file...")
    meshDomain = meshio.Mesh (points=domainNodes, cells=[("tetra", domainElements)])        # XPost (TOP) file
    meshDomain.write (outputDir+"/"+outputDomainFilename.replace(".top", ".exo"))           # Exodus file
    meshDomain.write (outputDir+"/"+outputDomainFilename.replace(".top", ".dat"))           # Tecplot file
    
    
    # EMBEDDED SURFACE MESH
    
    # Create empty dicts to store embedded mesh data
    meshEmbedded = {}
    labelEmbedded = {}
    
    # Set initial surface ID
    surfaceID = 1
    
    # Read input data
    print (">> Creating embedded surface mesh data for AERO-F...")
    for name, input in inputEmbedded.items():
        print (">> >> {}: Reading input mesh...".format (name))
        meshEmbedded[name] = meshio.read (input["filename"])
        
        # Create label for current surface
        labelEmbedded[name] = "{}_{}_{}".format (input["bcType"], name, surfaceID)
        surfaceID = surfaceID + 1
        
        # Apply scaling, translation, and rotation to mesh
        print (">> >> {}: Applying scaling, translation, and rotation to mesh...".format (name))
        meshEmbedded[name].points *= input["scaleFactor"]
        meshEmbedded[name].points = rotateMesh (meshEmbedded[name].points, input["rotation"])
        meshEmbedded[name].points += input["offset"]      
        
        # Write embedded surface
        print (">> >> {}: Writing embeded mesh data for each body...".format (name))
        meshEmbedded[name].write (os.path.join (outputDir, outputEmbeddedPrefix + name + "_tecplot.dat"))                           # Tecplot file   
        meshEmbedded[name].write (os.path.join (outputDir, outputEmbeddedPrefix + name + ".stl"))                                   # STL file
        meshEmbedded[name].write (os.path.join (outputDir, outputEmbeddedPrefix + name + ".exo"))                                   # Exodus file

    # Write TOP file and convert to Exodus
    print (">> Writing TOP and Exodus output for all bodies...")
    writeMeshEmbeddedAeroF (meshEmbedded, labelEmbedded, os.path.join (outputDir, outputEmbeddedPrefix + "Fluid"))       # XPost (TOP) file
    runXP2EXO (os.path.join (outputDir, outputEmbeddedPrefix + "Fluid"))
    writeMeshEmbeddedAeroS (meshEmbedded, inputEmbedded, os.path.join (outputDir, outputEmbeddedPrefix + "Structure"))   # AERO-S FEM file
    runFEM2TOP (outputDir, outputEmbeddedPrefix + "Structure")
    runXP2EXO (os.path.join (outputDir, outputEmbeddedPrefix + "Structure"))

    
    # USER MODULES

    # Run user-defined additional mesh modules
    if runMeshModule == "parachuteDGB":
        from meshGen_parachuteDGB import ParachuteClass
        print (">> Creating DGB Parachute surfaces and meshes...")
        payloadName = "Payload"                                                 # Set name for payload (must be in inputEmbedded defined above)
        Parachute = ParachuteClass()                                            # Instantiate class
        Parachute.main(outputDir, payloadName, meshEmbedded, inputEmbedded)     # Create parachute
    if runMeshModule == "parachuteRingsail":
        from meshGen_parachuteRingsail import ParachuteClass
        print (">> Creating Ringsail Parachute surfaces and meshes...")
        payloadName = "Payload"                                                 # Set name for payload (must be in inputEmbedded defined above)
        Parachute = ParachuteClass()                                            # Instantiate class
        Parachute.main(outputDir, payloadName, meshEmbedded, inputEmbedded)     # Create parachute
        

    # Print message when complete
    print (">> Finished!")
    
        

########################
# OTHER FUNCTIONS
######################## 
def writeMeshEmbeddedAeroF (meshEmbedded, labelEmbedded, filename):
    """ Write ascii file in XPost (TOP) format with embedded mesh data for AEROF
        - Assumes mesh is meshio object
        - Iterates over each item in mesh (i.e. multiple surfaces or bodies)
    """
    import numpy as np
    
    # Open output file
    f = open (filename+".top", "w")
    
    # Setup vars for iterating over multiple bodies
    nodeMin = {}
    nodeMaxGlobal = 0
    elementMin = {}
    elementMaxGlobal = 0
    
    # Get node and element min/max indices for each body
    for name, mesh in meshEmbedded.items():

        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
        elements = np.array (mesh.cells_dict["triangle"])
        
        # Get min/max indices for current body
        nodeMin[name] = nodeMaxGlobal
        nodeMaxGlobal = nodeMin[name] + len (points)
        elementMin[name] = elementMaxGlobal
        elementMaxGlobal = elementMaxGlobal + len (elements)        

    # Write NODE data
    f.write ("Nodes EmbeddedNodes\n")
    for name, mesh in meshEmbedded.items():

        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
                
        # Write node data
        for n in range (len (points)):
            f.write ("{:d}      {:.16e}  {:.16e}  {:.16e}\n".format (nodeMin[name]+n+1, points[n, 0], points[n, 1], points[n, 2]))
        
    # Write ELEMENTS data
    for name, mesh in meshEmbedded.items():
        f.write ("Elements {} using EmbeddedNodes\n".format (labelEmbedded[name]))

        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        elements = np.array (mesh.cells_dict["triangle"])
        
        # ELEMENTS
        elementType = 4     # triangle = 4
        for e in range (len (elements)):
            f.write ("{:d}   {:d}       {:d}  {:d}  {:d}\n".format (elementMin[name]+e+1, elementType, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))
    
    # Close output file
    f.close()
    

def writeMeshEmbeddedAeroS (meshEmbedded, inputEmbedded, filename):
    """ Write ascii file in XPost (TOP) format with embedded mesh data for AEROS
        - Assumes rigid structure with single node at CG, and phantom elements
        - Assumes mesh is meshio object
        - Iterates over each item in mesh (i.e. multiple surfaces or bodies)
    """
    import numpy as np

    # Open output files
    f = open (filename+".fem", "w")
    f_surfTopo = open (filename+"_surfaceTopo.fem", "w")
    f_nodeIvel = open (filename+"_nodeIvel.txt", "w")

    # NODES and IVEL
    f.write ("NODES\n")
    f_nodeIvel.write ("IVEL\n")
    nodeMaxGlobal = 0
    nodeMin = {}

    # Write CGs as first nodes
    for name, mesh in meshEmbedded.items():
        f.write ("* CG for {}:\n".format (name))
        f.write ("{:<12d} {:>25.16e} {:>25.16e} {:>25.16e}\n".format (nodeMaxGlobal+1, inputEmbedded[name]["cg"][0], inputEmbedded[name]["cg"][1], inputEmbedded[name]["cg"][2]))
        f_nodeIvel.write ("* CG IVEL for {}:\n".format (name))
        for dir in range (3):
            f_nodeIvel.write ("{} {}   0.0\n".format (nodeMaxGlobal+1, dir+1))     # init velocity to 0.0 --> update to desired vals when putting into AERO-S input deck
        nodeMaxGlobal = nodeMaxGlobal + 1
    
    # Write surface nodes
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
        nodeMin[name] = nodeMaxGlobal
        nodeMaxGlobal = nodeMin[name] + len (points)

        # Write node data
        f.write ("* Surface Nodes for {}:\n".format (name))
        for n in range (len (points)):
            f.write ("{:<12d}      {:>25.16e} {:>25.16e} {:>25.16e}\n".format (nodeMin[name]+n+1, points[n, 0], points[n, 1], points[n, 2]))
            
        # Write IVEL data
        f_nodeIvel.write ("* Surface Nodes IVEL for {}:\n".format (name))
        for dir in range (3):
            f_nodeIvel.write ("{} thru {}     {} 0.0\n".format (nodeMin[name]+1, nodeMaxGlobal, dir+1))
            
    

    # TOPOLOGY and SURFACETOPO
    f.write ("TOPOLOGY\n")
    elementMaxGlobal = 0
    
    # BEAMS
    # Iterate over each body 
    meshCount = 0
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
        
        # Write rigid beams - elemtype=133
        f.write ("* Rigid Beams for {}:\n".format (name))
        for n in range (len (points)):
            #         element#  type       node1   node2     
            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 133, meshCount+1, nodeMin[name]+n+1))
            elementMaxGlobal = elementMaxGlobal + 1
        
        # Iterate for next body
        meshCount = meshCount + 1
    
    # Save elementMax # for rigid beams
    elementMaxRigidBeams = elementMaxGlobal

    # SURFACES
    # Iterate over each body 
    meshCount = 0
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        elements = np.array (mesh.cells_dict["triangle"])
        
        # Write surface
        f.write ("* Phantom Shells for {}:\n".format (name))
        f_surfTopo.write ("* Surface for {}:\n".format (name))
        f_surfTopo.write ("SURFACETOPO {}\n".format (meshCount+1))
        for e in range (len (elements)):
            
            # TOPOLOGY - triangle is type=8 (phantom membrane)
            #         element#  type       node1   node2   node3
            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 8, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))        
             
            # SURFACETOPO - triangle is type=3
            #                  element#  type       node1   node2   node3
            f_surfTopo.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 3, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))        
            
            # Increment global count for elements
            elementMaxGlobal = elementMaxGlobal + 1      
            
        # Iterate for next body
        meshCount = meshCount + 1
    
    # Close topoloy and surfacetopo files
    f.close()
    f_surfTopo.close()
    f_nodeIvel.close()
    
    # Save elementMax # for Phantom Shells
    elementMaxPhantomShells = elementMaxGlobal        


    # ATTRIBUTES
    f = open (filename+"_attributes.fem", "w")
    f.write ("ATTRIBUTES\n")
    f.write ("* Rigid Beams:    elements {}-{}, material={}\n".format (1, elementMaxRigidBeams, 1))
    f.write ("* Phantom Shells: elements {}-{}, material={}\n".format (elementMaxRigidBeams+1, elementMaxPhantomShells, -1))
    f.write ("{} {} {}\n".format (1, elementMaxRigidBeams, 1))
    f.write ("{} {} {}\n".format (elementMaxRigidBeams+1, elementMaxPhantomShells, -1))
    f.write ("END\n")
    f.close()
    
    
    
def writeMeshEmbeddedAeroS_volume (meshEmbedded, inputEmbedded, filename):
    """ Write ascii file in XPost (TOP) format with embedded mesh data for AEROS
        - Assumes rigid structure with single node at CG, and phantom membranes for surface
        - Attributes and mass defined based on volume of tets
        - Assumes mesh is meshio object
        - Iterates over each item in mesh (i.e. multiple surfaces or bodies)
    """
    import numpy as np

    # Open output files
    f = open (filename+".fem", "w")
    f_surfTopo = open (filename+"_surfaceTopo.fem", "w")
    f_nodeIvel = open (filename+"_nodeIvel.txt", "w")

    # NODES and IVEL
    f.write ("NODES\n")
    f_nodeIvel.write ("IVEL\n")
    nodeMaxGlobal = 0
    nodeMin = {}

    # Write CGs as first nodes
    for name, mesh in meshEmbedded.items():
        f.write ("* CG for {}:\n".format (name))
        f.write ("{:<12d} {:>25.16e} {:>25.16e} {:>25.16e}\n".format (nodeMaxGlobal+1, inputEmbedded[name]["cg"][0], inputEmbedded[name]["cg"][1], inputEmbedded[name]["cg"][2]))
        f_nodeIvel.write ("* CG IVEL for {}:\n".format (name))
        for dir in range (3):
            f_nodeIvel.write ("{} {}   0.0\n".format (nodeMaxGlobal+1, dir+1))     # init velocity to 0.0 --> update to desired vals when putting into AERO-S input deck
        nodeMaxGlobal = nodeMaxGlobal + 1
    
    # Write surface nodes
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
        nodeMin[name] = nodeMaxGlobal
        nodeMaxGlobal = nodeMin[name] + len (points)

        # Write node data
        f.write ("* Surface Nodes for {}:\n".format (name))
        for n in range (len (points)):
            f.write ("{:<12d}      {:>25.16e} {:>25.16e} {:>25.16e}\n".format (nodeMin[name]+n+1, points[n, 0], points[n, 1], points[n, 2]))
            
        # Write IVEL data
        f_nodeIvel.write ("* Surface Nodes IVEL for {}:\n".format (name))
        for dir in range (3):
            f_nodeIvel.write ("{} thru {}     {} 0.0\n".format (nodeMin[name]+1, nodeMaxGlobal, dir+1))
            
    
    # CREATE VOLS
    
    # Create volume tets and calc volume fraction (will be used later to determine material mass to achieve desired CG)
    # Iterate over each body
    meshCount = 0
    bodyDict = {}
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        elements = np.array (mesh.cells_dict["triangle"])
        eMax = len (elements)
        
        # Create list of tets
        bodyDict[meshCount] = {}
        bodyDict[meshCount]["tetNodes"] = np.empty ((eMax, 4))
        bodyDict[meshCount]["tetNodes"][:, 0] = meshCount+1
        bodyDict[meshCount]["tetNodes"][:, 1] = nodeMin[name] + elements[:, 0] + 1
        bodyDict[meshCount]["tetNodes"][:, 2] = nodeMin[name] + elements[:, 1] + 1
        bodyDict[meshCount]["tetNodes"][:, 3] = nodeMin[name] + elements[:, 2] + 1
        
        
        # Iterate over each element for current body
        for e in range (len (elements)):
            
            # Set node numbers for current tet assuming CG as first point, and surface element as other three points
            nodeList = [meshCount+1, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1]
            volumeTmp.append (nodeList)
            
            # Increment global count for elements
            elementMaxGlobal = elementMaxGlobal + 1      
        
        # Iterate for next body
        meshCount = meshCount + 1
    
    # Save elementMax # for rigid beams
    elementMaxRigidBeams = elementMaxGlobal
    
   
    # TOPOLOGY and SURFACETOPO
    f.write ("TOPOLOGY\n")
    elementMaxGlobal = 0
    
    # VOLUME TETS
    # Iterate over each body 
    meshCount = 0
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        points = np.array (mesh.points)
        
        # Write rigid tetrahedra assuming surface triangle + CG point
        f.write ("* Rigid tetrahedra for {}:\n".format (name))
        for e in range (len (elements)):
            
            # TOPOLOGY - volume tetrahedra is type=150 (3d 4-node rigid tetrahedron with 3 dof/node )
            #         element#  type       node1   node2   node3
            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 129, meshCount+1, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))        
             
            # Increment global count for elements
            elementMaxGlobal = elementMaxGlobal + 1      
        
        # Iterate for next body
        meshCount = meshCount + 1
    
    # Save elementMax # for rigid beams
    elementMaxRigidBeams = elementMaxGlobal

    # SURFACES
    # Iterate over each body 
    meshCount = 0
    for name, mesh in meshEmbedded.items():
        
        # Copy meshio data numpy arrays to speed up output - note for some reason meshio can be very slow when accessing the data for output on large meshes
        elements = np.array (mesh.cells_dict["triangle"])
        
        # Write surface
        f.write ("* Phantom Membranes for {}:\n".format (name))
        f_surfTopo.write ("* Surface for {}:\n".format (name))
        f_surfTopo.write ("SURFACETOPO {}\n".format (meshCount+1))
        for e in range (len (elements)):
            
            # TOPOLOGY - triangle is type=129 (phantom membrane = 3d 3-node plane stress triangular membrane with 3 dof/node)
            #         element#  type       node1   node2   node3
            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 129, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))        
             
            # SURFACETOPO - triangle is type=3
            #                  element#  type       node1   node2   node3
            f_surfTopo.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementMaxGlobal+1, 3, nodeMin[name]+elements[e, 0]+1, nodeMin[name]+elements[e, 1]+1, nodeMin[name]+elements[e, 2]+1))        
            
            # Increment global count for elements
            elementMaxGlobal = elementMaxGlobal + 1      
            
        # Iterate for next body
        meshCount = meshCount + 1
    
    # Close topoloy and surfacetopo files
    f.close()
    f_surfTopo.close()
    f_nodeIvel.close()
    
    # Save elementMax # for Phantom Shells
    elementMaxPhantomShells = elementMaxGlobal        


    # ATTRIBUTES
    f = open (filename+"_attributes.fem", "w")
    f.write ("ATTRIBUTES\n")
    f.write ("* Rigid Beams:    elements {}-{}, material={}\n".format (1, elementMaxRigidBeams, 1))
    f.write ("* Phantom Shells: elements {}-{}, material={}\n".format (elementMaxRigidBeams+1, elementMaxPhantomShells, -1))
    f.write ("{} {} {}\n".format (1, elementMaxRigidBeams, 1))
    f.write ("{} {} {}\n".format (elementMaxRigidBeams+1, elementMaxPhantomShells, -1))
    f.write ("END\n")
    f.close()
    
    

if __name__ == "__main__":
    meshGenMain ()
    