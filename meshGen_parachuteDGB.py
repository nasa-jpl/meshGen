#####################################################################################
# MESHGEN PARACHUTE DGB
#####################################################################################
# Mesh generation of DGB (Disk-Gap Band) Parachute for FSI simulations in AERO suite
# - Creates "as-built" configuration with flat disk and cylindrical band
# - Uses gmsh for meshing
# - Includes routine for calculating folding displacements
# - Parses Excel file with parachute parameters
#####################################################################################


# Imports
from meshGen_parachuteUtils import unit_vector, createRotationMatrix, calcTriangleNormal, calcQuadNormal, calcVectorRotationAngles
import gmsh


# Globals
dimLine = 1
dimSurf = 2
dimVol  = 3
flagEmbedLines          = False                     # If parachute mesh = True, embed suspension lines (i.e. create dressing surfaces)?
flagRecombine           = False                     # False = triangular surface elements, True = quad surface elements
flagWriteTecplot        = False                     # Write Tecplot files for parachute mesh? (can be slow)
flagWriteFrames         = True                      # Write c/eframes for AERO?
flagWriteSurfaceTopo    = True                      # Write surface topology files for AERO?
flagWriteDecomp         = True                      # Write decomposition files for AERO?
flagVentForeshorten     = True                      # shorten vent lines for foreshortening? [Note applied in analytical folding function]
flagCoarsenPayloadMesh  = True                      # coarsen payload mesh for structure?
foldType                = "triangle"                # triangle or parabolic for gore shape - note parabolic will put gores closer together near radials... may result in more self-contact
foldDist                = 0.2                       # [m] folding distance input - i.e. distance of band radials to z-axis when folded (smaller will be tighter folding)
apparateAngleX          = 0.0                       # [deg] rotate parachute system about riser start point by specified degrees (x-axis rotation)
apparateAngleY          = 0.0                       # [deg] rotate parachute system about riser start point by specified degrees (y-axis rotation)
flagTransfiniteMesh     = False                      # Do transfinite mesh (structured triangles)? Or more unstructured appraoch?
EPS                     = 1.0e-14                   # A very small number


# Class definition
class ParachuteClass:
    """ Parachute mesh constructor class for FSI simulations
    """
    def __init__(self):

        self.excelInputFilename = "FSI_materialProps_v11_for_distribution.xlsx"
        self.excelInputSheetname = "ASPIRE_SR02"
        self.param = {}                 # parameters defining parachute
        self.mesh = {}                  # meshio representation
        self.data = {}                  # internal data respresentation of mesh - used due to issues with meshio output
        self.meshType = {}              # fluid or structure
        self.attributes = {}            # structure model attributes dict
        self.minPointsCircle = 20       # min points in a circle/ellipse
        self.minPointsEdge = 5          # min points along an edge - approximate as actual spacing set at nodes which average two edge spacings
        self.BCtype = "StickMoving"     # StickMoving for AERO-F surface BC input
        self.goresBandPerSubdomain = 8  # number of gores per subdomain in AERO-S decomposition file output
        self.goresDiskPerSubdomain = 4  # number of gores per subdomain in AERO-S decomposition file output
        self.numOtherSubdomain = 8      # number of other (non-gore) subdomains in AERO-S (i.e. rigging, payload, etc.)
        self.nBridleSupportBeams = 3    # number of support beams to help transfer bridle load to payload (number per bridle attach point)
        self.payloadMeshSize = 0.6      # [m] Approximate mesh size reparameterization of payload for structure mesh
        self.verbosity = 0              # verbosity level for gmsh --> 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug

    

    def main (self, outputDir, payloadName, meshEmbedded, inputEmbedded):
        """ Run through tasks to create parachute model/mesh
        """
        import os
        from runMeshGen import outputEmbeddedPrefix
        from meshGen_utils import runXP2EXO, runFEM2TOP

        # Set filename for stl file
        payloadStlFilename = os.path.join (outputDir, outputEmbeddedPrefix + payloadName + ".stl")

        # Parse Excel file for inputs
        self.parseInputExcel (outputDir+"/"+"out_ParachuteSettingsData.txt")

        # Create mesh
        self.createMesh ("fluid", flagEmbedLines)    
        self.addMeshPayload ("fluid", meshEmbedded[payloadName], inputEmbedded[payloadName]["cg"])
        self.createMesh ("structure", flagEmbedLines, flagRecombine)
        self.addMeshPayload ("structure", meshEmbedded[payloadName], inputEmbedded[payloadName]["cg"], payloadStlFilename)
        self.calcElementFrames ("structure")
        
        # Write data output - fluid mesh
        print (">> Writing Parachute Fluid Mesh Data...")
        self.writeMeshAeroF (outputDir+"/"+"parachute_meshFluid_combined", flagCombineSurfs=True)
        self.writeMeshAeroF (outputDir+"/"+"parachute_meshFluid_individualGores", flagCombineSurfs=False)
        runXP2EXO (outputDir+"/"+"parachute_meshFluid_combined")
        runXP2EXO (outputDir+"/"+"parachute_meshFluid_individualGores")
        
        # Write data output - Tecplot
        if flagWriteTecplot == True:
            print (">> Writing Tecplot data...")
            self.writeMeshTecplot (outputDir+"/"+"parachute_meshFluid", meshType="fluid", flagWriteFrames=False)
            self.writeMeshTecplot (outputDir+"/"+"parachute_meshStructure", meshType="structure", flagWriteFrames=True)
        
        # Write data output - structure mesh (payload + parachute)
        print (">> Writing Parachute Structure Mesh Data...")
        self.writeMeshAeroS (outputDir+"/"+"parachute_meshStructure", flagPayloadOnly=False)
        runFEM2TOP (outputDir, "parachute_meshStructure")
        runXP2EXO (outputDir+"/"+"parachute_meshStructure")
        
        # Write data output - structure mesh (payload only)
        print (">> Writing PayloadOnly Structure Mesh Data...")        
        self.writeMeshAeroS (outputDir+"/"+"payloadOnly_meshStructure", flagPayloadOnly=True)
        runFEM2TOP (outputDir, "payloadOnly_meshStructure")
        runXP2EXO (outputDir+"/"+"payloadOnly_meshStructure")        

        

    def createMesh (self, meshType, flagEmbedLines=True, flagRecombine=False):
        """ Create parachute mesh
        """
        import numpy as np
        import os
        
        # General settings
        var = self.param                 # Save params to variable      

        # Init pygmsh models to collect all elements
        gmsh.initialize()
        gmsh.model.add(meshType)

        # Set verbosity level
        gmsh.option.setNumber("General.Verbosity", self.verbosity)

        # BRIDLE
        name = "Bridle"
        p0 = [var["Bridle"]["l_proj"], 0.0, 0.0]
        p1 = [0.0, 0.0, var["Bridle"]["h_proj"]]
        if meshType == "structure"  or  (meshType == "fluid" and flagEmbedLines == True):
            self.createMeshRadials (name=name, ds=0.05, pStart=p0, pEnd=p1, radius=var["Bridle"]["r"], factorSpacing=0.05, iMax=var["Bridle"]["n_lines"], meshType=meshType)
        self.attributes[name] = {"matAttribute": 1, "elementType": 191, "elementName": "NLBAR"}
        
        # RISER
        name = "Riser"
        p0 = p1
        p1 = [0.0, 0.0, p1[2]+var["Riser"]["L"]]
        if meshType == "structure"  or  (meshType == "fluid" and flagEmbedLines == True):
            self.createMeshRadials (name=name, ds=0.05, pStart=p0, pEnd=p1, radius=var["Riser"]["r"], factorSpacing=0.05, iMax=1, meshType=meshType)
        self.attributes[name] = {"matAttribute": 2, "elementType": 191, "elementName": "NLBAR"}
        pRiserStart = p0
                                 
        # SUSPENSION LINES
        name = "SuspensionLines"
        dx = var["Band"]["d"]/2.0
        dz = np.sqrt ((var["SuspensionLines"]["L"])**2.0 - dx**2.0)    
        p0 = p1
        p1 = [dx, 0.0, p1[2]+dz]
        if meshType == "structure"  or  (meshType == "fluid" and flagEmbedLines == True):
            self.createMeshRadials (name=name, ds=0.05, pStart=p0, pEnd=p1, radius=var["SuspensionLines"]["r"], factorSpacing=0.05, iMax=var["SuspensionLines"]["n_lines"], meshType=meshType)
        self.attributes[name] = {"matAttribute": 3, "elementType": 191, "elementName": "NLBAR"}
        pSuspensionLineApex = p0
                                 
        # BAND - Radials
        name = "BandRadials"
        p0 = p1
        p1 = [p1[0], p1[1], p1[2]+var["Band"]["L"]]
        self.attributes[name] = {"matAttribute": 4, "elementType": 191, "elementName": "NLBAR"}
        pBandStart = p0
        pBandEnd = p1
        
        # BAND - Circumferential Reinforcements
        name = "BandLE"
        rLoc0 = var["Band"]["d"] / 2.0
        zLoc0 = var["Bridle"]["h_proj"] + var["Riser"]["L"] + np.sqrt ((var["SuspensionLines"]["L"])**2.0 - rLoc0**2.0)   
        self.attributes[name] = {"matAttribute": 6, "elementType": 191, "elementName": "NLBAR"}
        name = "BandTE"
        rLoc1 = rLoc0
        zLoc1 = zLoc0 + var["Band"]["L"]
        self.attributes[name] = {"matAttribute": 7, "elementType": 191, "elementName": "NLBAR"}
        rBand = [rLoc0, rLoc1]      # save coords - gore mesh will be created at end
        zBand = [zLoc0, zLoc1]      # save coords - gore mesh will be created at end     
             
        # GAP LINES
        name = "GapLines"
        p0 = p1
        p1 = [p1[0], p1[1], p1[2]+var["Gap"]["L"]]
        if meshType == "structure"  or  (meshType == "fluid" and flagEmbedLines == True):
            self.createMeshRadials (name=name, ds=0.05, pStart=p0, pEnd=p1, radius=var["SuspensionLines"]["r"], factorSpacing=0.05, iMax=var["SuspensionLines"]["n_lines"], meshType=meshType)
        self.attributes[name] = {"matAttribute": 3, "elementType": 191, "elementName": "NLBAR"}

        # DISK - Radials
        name = "DiskRadials"
        p0 = p1
        p1 = [var["Vent"]["d"]/2, 0.0, p1[2]]
        self.attributes[name] = {"matAttribute": 5, "elementType": 191, "elementName": "NLBAR"}
        pDiskStart = p0
        pDiskEnd = p1
                                 
        # DISK - Circumferential Reinforcements
        name = "DiskLE"
        rLoc0 = p0[0]
        zLoc0 = p1[2]
        self.attributes[name] = {"matAttribute": 8, "elementType": 191, "elementName": "NLBAR"}
        name = "DiskTE"
        rLoc1 = p1[0]
        zLoc1 = zLoc0
        self.attributes[name] = {"matAttribute": 9, "elementType": 191, "elementName": "NLBAR"}  
        rDisk = [rLoc0, rLoc1]      # save coords - gore mesh will be created at end
        zDisk = [zLoc0, zLoc1]      # save coords - gore mesh will be created at end

        # VENT LINES
        name = "VentLines"
        p0 = p1
        p1 = [0.0, 0.0, p1[2]]
        if meshType == "structure"  or  (meshType == "fluid" and flagEmbedLines == True):
            self.createMeshRadials (name=name, ds=0.05, pStart=p0, pEnd=p1, radius=var["SuspensionLines"]["r"], factorSpacing=0.05, iMax=var["SuspensionLines"]["n_lines"], meshType=meshType)
        self.attributes[name] = {"matAttribute": 10, "elementType": 191, "elementName": "NLBAR"} 
        pVentApex = p1

        # # CREATE PANELS
        
        # Set gore element type
        if flagRecombine == True:
            elementType = 128
            elementName = "QUAD4"
        else:
            elementType = 129
            elementName = "TRI3"
        
        # BAND - Gores
        name = "BandGores"
        iMax = var["SuspensionLines"]["n_lines"]
        for i in range (iMax):
            self.attributes[name+"_{}".format(i)] = {"matAttribute": 50, "elementType": elementType, "elementName": elementName}
        self.createMeshPanels (name=name, radialName="BandRadials", circumName=["BandLE", "BandTE"], ds=0.05, rLoc=rBand, zLoc=zBand, iMax=var["SuspensionLines"]["n_lines"], meshType=meshType)
   
        # DISK - Gores
        name = "DiskGores"
        iMax = var["SuspensionLines"]["n_lines"]
        for i in range (iMax):
            self.attributes[name+"_{}".format(i)] = {"matAttribute": 51, "elementType": elementType, "elementName": elementName, "flagAlternateGores": True}
        self.createMeshPanels (name=name, radialName="DiskRadials", circumName=["DiskLE", "DiskTE"], ds=0.05, rLoc=rDisk, zLoc=zDisk, iMax=var["SuspensionLines"]["n_lines"], meshType=meshType)

        
        # MESH GEOMETRY
        
        # Apply min number of points to curved lines
        gmsh.option.setNumber ("Mesh.MeshSizeFromCurvature", self.minPointsCircle)
        
        # Apply transfinite method to all entities
        if flagTransfiniteMesh == True:
            gmsh.model.mesh.setTransfiniteAutomatic(recombine=flagRecombine)
        
        # Else set meshing algorithm:
        #  1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
        else:
            entityList = gmsh.model.getEntities (dimSurf)               # Get list of entity (dim, tag) pairs for 2D surfaces
            for entity in entityList:
                gmsh.model.mesh.setAlgorithm (entity[0], entity[1], 6)
        
        # Mesh structure and output
        print (">> >> >> Meshing {}...".format (meshType))
        gmsh.model.mesh.generate(dimSurf)
        
        # Remove duplicate nodes
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()
        
        # If quad mesh, start up gmsh GUI to inspect - note need to work more output/etc. for future
        if flagRecombine == True:
            gmsh.model.mesh.recombine()
            #gmsh.fltk.run()
            
        
        # CONVERT TO INTERNAL DATA    
        
        # Create meshGen data structure to hold data
        print (">> >> >> Converting mesh to internal meshGen format...")
        self.createMeshGenData (meshType)
        
        
        # CALCULATE FOLDING DISPLACEMENTS AND ROTATIONS
        if meshType == "structure":
            print (">> >> >> Calculating mesh displacements for folding...")
            self.createFoldParachuteDisplacements (foldDist, pRiserStart, pSuspensionLineApex, pBandStart, pBandEnd, pDiskStart, pDiskEnd, pVentApex)
            print (">> >> >> Calculating mesh rotations for folding...")
            self.calcFoldParachuteRotations()
               
        # Finalize mesh
        gmsh.finalize()
     


    def createFoldParachuteDisplacements (self, foldDist, pRiserStart, pSuspensionLineApex, pBandStart, pBandEnd, pDiskStart, pDiskEnd, pVentApex):
        """ Create DISPLACEMENTS file for structure sim to fold parachute
            - Note: assumes parachute mesh has already been created and meshGen data structure exists
            - Inputs
                foldDist:                                                                       min distance from z axis for band radials (smaller will be tighter folding)
                pSuspensionLineApex, pBandStart, pBandEnd, pDiskStart, pDiskEnd, pVentApex:     points defining x,y,z coordinates for folding regions of interest
        """
        import numpy as np
        from math import pi
         
         
        # Function defs for coordinate conversion
        def cartesianToSpherical (nodeCartesian):
            """ Convert cartesian coord (x,y,z) to spherical (r, theta, phi)
            """
            import numpy as np
            x, y, z = nodeCartesian[0], nodeCartesian[1], nodeCartesian[2]
            r = np.sqrt (x**2 + y**2 + z**2)                            # radius
            if abs (r) <= EPS:
                return [0, 0, 0]                                        # if point at origin, return zeros
            theta = np.acos (z / r)                                     # inclination
            phi = np.atan2 (y, x)  #np.sign (y) * np.acos (x / np.sqrt (x**2 + y**2))     # azimuth
            return [r, theta, phi]
            
        def sphericalToCartesian (nodeSpherical):
            """ Convert spherical coord (r, theta, phi) to cartesian (x, y, z)
            """
            import numpy as np
            r, theta, phi = nodeSpherical[0], nodeSpherical[1], nodeSpherical[2]
            x = r * np.sin (theta) * np.cos (phi)
            y = r * np.sin (theta) * np.sin (phi)
            z = r * np.cos (theta)
            return [x, y, z]

        
        # SETUP FOLDING PARAMS
        
        # Tolerance for checking regions - this should be << than min expected spacing of lines/etc., but > EPS to allow some precision effects in meshing
        regionTol = 1.0e-6
        
        # Copy params to local var
        var = self.param
        
        # Calculate as-fabricated suspension line angle and rDist from z-axis, and offset in r and z for new folded band start
        length = var["SuspensionLines"]["L"]
        rDistOrig = var["Disk"]["d"] / 2.0
        thetaOrig = np.asin (rDistOrig / length)
        theta_lines = np.asin (foldDist / length)
        dr_lines = foldDist - rDistOrig
        dz_lines = length * (np.cos (theta_lines) - np.cos (thetaOrig))
        
        # Calculate angle for folding to align to vent apex
        length = var["Disk"]["d"] / 2.0
        theta_disk = pi - np.asin (foldDist / length)           # note disk is folded in negative z dir, so will be between 90 and 180 deg
        dz_disk = abs (length * np.cos (theta_disk))            # dz_disk will be positive as offset from current disk location to folding apex
        
        # Get number of gores/lines and spacing phi around parachute
        iMax = self.param["SuspensionLines"]["n_lines"]     
        phi_spacing = 2 * pi / iMax
        
        # Get vent foreshortening fraction for vent line shrinkage
        lengthVent = var["Vent"]["d"] / 2                       # length of radius for ventTE (non-foreshortened)
        lengthVentLine = var["VentLines"]["L"]                  # length of foreshortened vent line
        hVent = lengthVent * np.sin (theta_disk)                # height (dist from z-axis) of folded ventTE
        theta_vent = pi - np.asin (hVent / lengthVentLine)      # note vent is folded in negative z dir, so will be between 90 and 180 deg
        
        
        # Init tmp dict to store data, get node data, and init displacement output
        data = self.data["structure"]
        node = data["node"]
        nodeMax = len (node)
        displac = np.zeros ((nodeMax, 3))
        
                
        # CALC FOLDING DISPLACEMENTS
        
        # Iterate over each node and figure out displacements based on regions
        for n in range (nodeMax):
            

            # SKIP: PRE-FOLDING REGION
            if node[n, 2] < pSuspensionLineApex[2]:
                nodeOrig = node[n, :]
                nodeFold = nodeOrig         # skip folding

                
            # SUSPENSION LINE REGION
            elif node[n, 2] < (pBandStart[2] - regionTol):
            
                # Calc new node location based on folding angle
                nodeOrig = node[n, :] - pSuspensionLineApex                         # Get node location in local frame (origin at apex)
                [r, theta, phi] = cartesianToSpherical (nodeOrig)                   # Convert to spherical coordinates
                nodeFold = np.array (sphericalToCartesian ([r, theta_lines, phi]))  # Get cartesian coord for folded node
                
            
            # BAND / GAP REGION
            elif node[n, 2] < (pDiskStart[2] - regionTol):
                
                # Set params of interest for region
                rTmp = var["Band"]["d"] / 2.0
                
                # Get node location in local frame (origin at current z location but on z-axis)
                nodeOrig = node[n, :] - np.array([0, 0, node[n, 2]])                # Get node location in local frame (origin at z-axis, current slice)
                [r, theta, phi] = cartesianToSpherical (nodeOrig)                   # Convert to spherical coordinates - note theta should be 90 deg
                
                # If phi indicates gap line or band gore edge, set adjustment for gore folding to zero
                if abs (phi % phi_spacing) <= regionTol  or  abs (abs (phi % phi_spacing) - phi_spacing) <= regionTol:
                    rAdjFold = 0.0
                
                # Else calculate radius adjustment for folded gore shape
                else:
                    
                    # Get bounding line phi angles for current node, and calc nodes on edges and gore width
                    phi1 = phi_spacing * np.floor (phi / phi_spacing)                       # phi corresponding to edge 1
                    phi2 = phi_spacing * np.ceil (phi / phi_spacing)                        # phi corresponding to edge 2
                    node1 = sphericalToCartesian ([rTmp, theta, phi1])                      # node on edge 1
                    node2 = sphericalToCartesian ([rTmp, theta, phi2])                      # node on edge 2
                    widthOrig = np.linalg.norm (np.array(node2) - np.array(node1))          # width of gore
                    widthLocFrac = np.linalg.norm (nodeOrig - np.array(node1)) / widthOrig  # fraction of width of location of node
                    
                    # Get radial adjustment for location on gore based on folding (assume triangular or parabolic)
                    node1 = sphericalToCartesian ([rTmp + dr_lines, theta, phi1])           # node on edge 1 - folded
                    node2 = sphericalToCartesian ([rTmp + dr_lines, theta, phi2])           # node on edge 2 - folded
                    widthFold = np.linalg.norm (np.array(node2) - np.array(node1))          # width of gore - folded
                    W = widthOrig
                    b = widthFold
                    x = widthLocFrac * b
                    if foldType == "triangle":
                        rAdjFold = np.sqrt((W/2)**2 - (b/2)**2) * (1 - 2/b * abs (x - b/2))                                 # calc height along fold assuming parametric equation for isosceles triange fold of original width
                    else:
                        rAdjFold = 4/(b**2) * (W/2 - b**2 / (8*W) * np.log (4*W/b) - (b**2) / (16*W)) * (x * (b - x))       # else assume parabolic curve - use higher-order logarithmic approximation for high curvature arc
                    
                # Apply offsets for folding and get new node location
                nodeFold = sphericalToCartesian ([r + dr_lines + rAdjFold, theta, phi])     # location of folded node in z slice
                nodeFold[2] = nodeFold[2] + dz_lines                                        # offset z slice based on folded suspension lines
                                
            
            # DISK / VENT REGION
            else:

                # Get node location in local frame (origin at current z location but on z-axis)
                nodeOrig = node[n, :] - np.array([0, 0, node[n, 2]])                        # Get node location in local frame (origin at z-axis, current slice)
                [r, theta, phi] = cartesianToSpherical (nodeOrig)                           # Convert to spherical coordinates - note theta should be 90 deg
                
                # If phi indicates gap line or band gore edge, set adjustment for gore folding to zero
                if abs (phi % phi_spacing) <= regionTol  or  abs (abs (phi % phi_spacing) - phi_spacing) <= regionTol:
                    rAdjFold = 0.0
                    
                # Else calculate radius adjustment for folded gore shape
                else:
                    
                    # Get bounding line phi angles for current node, and calc nodes on edges and gore width
                    phi1 = phi_spacing * np.floor (phi / phi_spacing)                       # phi corresponding to edge 1
                    phi2 = phi_spacing * np.ceil (phi / phi_spacing)                        # phi corresponding to edge 2
                    node1 = sphericalToCartesian ([r, theta, phi1])                         # node on edge 1 - note approx as it doesn't account for flat LE between gores... but should be close
                    node2 = sphericalToCartesian ([r, theta, phi2])                         # node on edge 2 - note approx as it doesn't account for flat LE between gores... but should be close
                    widthOrig = np.linalg.norm (np.array(node2) - np.array(node1))          # width of gore
                    widthLocFrac = np.linalg.norm (nodeOrig - np.array(node1)) / widthOrig  # fraction of width of location of node            
                    
                    # Get radial adjustment for location on gore based on folding (assume triangular or parabolic)
                    node1 = sphericalToCartesian ([r, theta_disk, phi1])                    # node on edge 1 - folded - note approx as it doesn't account for flat LE between gores... but should be close
                    node2 = sphericalToCartesian ([r, theta_disk, phi2])                    # node on edge 2 - folded - note approx as it doesn't account for flat LE between gores... but should be close
                    widthFold = np.linalg.norm (np.array(node2) - np.array(node1))          # width of gore - folded
                    W = widthOrig
                    b = widthFold
                    x = widthLocFrac * b
                    if foldType == "triangle":
                        rAdjFold = np.sqrt((W/2)**2 - (b/2)**2) * (1 - 2/b * abs (x - b/2))                                 # calc height along fold assuming parametric equation for isosceles triange fold of original width
                    else:
                        rAdjFold = 4/(b**2) * (W/2 - b**2 / (8*W) * np.log (4*W/b) - (b**2) / (16*W)) * (x * (b - x))       # else assume parabolic curve - use higher-order logarithmic approximation for high curvature arc
                
                # Get adjustments for folding gore shape if necessary (i.e. if not edge or line)
                if r <= EPS  or  rAdjFold == 0:
                    thetaAdjFold = 0.0
                else:
                    thetaAdjFold = np.atan (rAdjFold / r)                                   # get angle to raise gore point based on folding height of point
                
                # Fold disk and apply offsets for folding to get new node location
                nodeFold = sphericalToCartesian ([r, theta_disk - thetaAdjFold, phi])       # location of folded node assuming gore edge or line
                nodeFold[2] = nodeFold[2] + dz_lines + dz_disk                              # offset z slice based on disk folding moving vent apex aft
                
                # Apply vent foreshortening to shorten vent lines
                if flagVentForeshorten == True:
                    hTmp = np.sqrt (nodeFold[0]**2.0 + nodeFold[1]**2.0)                    # get current point radial distance from z-axis
                    if hTmp < hVent - regionTol:
                        hInv = hVent - hTmp                                                 # distance from ventTE in radial direction to z-axis
                        dzAdj = hInv / np.tan (theta_disk) - hInv / np.tan (theta_vent)     # dz adjustment for shorter vent lines
                        nodeFold[2] = nodeFold[2] + dzAdj                                   # offset z based on foreshortening
                        
            
            # SAVE DISPLACEMENT
        
            # Save displacement from original to folded node
            displac[n, :] = nodeFold - nodeOrig


            # APPLY APPARATE ANGLE
            if node[n, 2] >= pRiserStart[2]:
                
                # Get folded node point
                nodeOrig = node[n] + displac[n]
           
                # Rotate point
                from meshGen_utils import rotateMesh
                nodeApparate = rotateMesh (nodeOrig.reshape(1,3), [apparateAngleX, apparateAngleY, 0], rotOrigin=pRiserStart)

                # Add additional displacement
                displac[n, :] = displac[n, :] + (nodeApparate.reshape(3) - nodeOrig)


        # Update class data with results
        self.data["structure"]["displac"] = displac
            
        
        

    def calcFoldParachuteRotations (self):
        """ Calculate node rotation angles from original to displaced mesh - required for IDISP file
        """
        import numpy as np
        
        # Init tmp dict to store data, get node data, and init displacement output
        data = self.data["structure"]
        node = data["node"]
        nodeMax = len (node)
        displac = data["displac"]
        rotation = np.zeros ((nodeMax, 3))
        weight = np.zeros ((nodeMax))
        nodeDisplac = node + displac
        
        # Iterate over each set/element to calculate node rotation angles
        for setName, setData in data["cell_set"].items():
            print (">> >> >> >> Calculating rotations for: {}".format (setName))
            for element in setData["elements"]:
                
                # If line element, calc vectors of orig and displaced mesh
                if setData["type"] == "line":
                    
                    # Calculate vectors defining orientation of orig and displaced line, and get axis-angle rotation vector for rotation of element
                    v1 = np.array(node[element[1]]) - np.array(node[element[0]])
                    v2 = np.array(nodeDisplac[element[1]]) - np.array(nodeDisplac[element[0]])
                    rotVec = calcVectorRotationAngles (v1, v2, type="rotVec")
                   
                # Else if triangular element, calculate vectors using two sides
                elif setData["type"] == "triangle":
                    
                    # Calculate normal vectors of orig and displaced triangle, and get axis-angle rotation vector for rotation of element
                    v1 = calcTriangleNormal (node[element[0]], node[element[1]], node[element[2]])
                    v2 = calcTriangleNormal (nodeDisplac[element[0]], nodeDisplac[element[1]], nodeDisplac[element[2]])
                    rotVec = calcVectorRotationAngles (v1, v2, type="rotVec")
                   
                # Else if triangular element, calculate vectors using two sides
                elif setData["type"] == "quad":
                    
                    # Calculate normal vectors of orig and displaced triangle, and get axis-angle rotation vector for rotation of element
                    v1 = calcQuadNormal (node[element[0]], node[element[1]], node[element[2]], node[element[3]])
                    v2 = calcQuadNormal (nodeDisplac[element[0]], nodeDisplac[element[1]], nodeDisplac[element[2]], nodeDisplac[element[3]])
                    rotVec = calcVectorRotationAngles (v1, v2, type="rotVec")
                    
                # Save rotation angles and increment weight (counter) to later calc average angles for each node
                for n in element:
                    rotation[n, :] += np.array (rotVec)
                    weight[n] += 1
        
        # Average all element rotations applied to each node to calc final node rotation
        nFilt = weight[:] > 0
        rotation[nFilt, :] = rotation[nFilt, :] / weight[nFilt][:, None]
        
        # Save final rotations to data class
        self.data["structure"]["rotation"] = rotation
                   
  
  
  

    def addMeshPayload (self, meshType, payloadMesh, payloadCG, payloadStlFilename=None):
        """ Add payload into mesh
            - Note: assumes parachute mesh has already been created and meshGen data structure exists
            - Inputs
                payloadMesh:        meshio format of triangulated surface mesh [i.e. class with .points and .elements["triangle"] data
                payloadCG:          array of [cg_x, cg_y, cg_z]
        """
        import numpy as np
        from math import pi
        
        # Error checking
        if "node" not in self.data[meshType]:
            raise Exception ("Error in addMeshPayload: 'node' dict not defined yet - make sure parachute mesh is created first.")
            
        # Init tmp dict to store data
        data = self.data[meshType]
      
        # Copy payload meshio data to numpy arrays - assumes meshio has been used to read STL file earlier
        name = "Payload"
        points = np.array (payloadMesh.points)
        elements = np.array (payloadMesh.cells_dict["triangle"])
                             
        # FLUID 
        if meshType == "fluid":

            # Set offsetCG and get max number of nodes in payload mesh
            offsetCG = 0
            nodeMaxPayload = len (points)
            print (">> >> >> >> Adding Payload: nodeMax = {}".format (nodeMaxPayload))
        
            # Insert payload nodes to front of meshGen dataset
            data["node"] = np.append (points, data["node"], axis=0)

            # Offset all parachute cell_set element node ids to accound for payload nodes
            for setName in data["cell_set"].keys():
                data["cell_set"][setName]["elements"] = np.array (data["cell_set"][setName]["elements"])            # convert fluid elements from numpy format if not already
                data["cell_set"][setName]["elements"] += nodeMaxPayload                                             # offset elements
                data["cell_set"][setName]["elements"] = data["cell_set"][setName]["elements"].tolist()              # need to convert fluid elements from numpy back to list format
                
            # Create set for payload surface
            setName = name
            data["cell_set"][setName] = {}
            data["cell_set"][setName]["type"] = "triangle"
            data["cell_set"][setName]["elements"] = elements + offsetCG     # Add payload surface elements - offset if necessary to account for CG as first node

        
        # STRUCTURE
        if meshType == "structure":
            
            # Coarsen payload mesh if desired (can help improve AERO-S stability time step)
            if flagCoarsenPayloadMesh == True:
                print (">> >> >> >> Creating coarse payload mesh for structure")
                
                # Use gmsh to reclassify mesh and convert to parameterized surface, and remesh to coarser resolution
                import gmsh
                gmsh.initialize()
                gmsh.option.setNumber("General.Verbosity", self.verbosity)                          # Set verbosity level
                gmsh.model.add(meshType)
                gmsh.merge (payloadStlFilename)                                                     # read STL
                gmsh.model.mesh.classifySurfaces (pi/8, boundary=True, forReparametrization=True)   # reparameterize STL into new surface(s) to remesh
                gmsh.model.mesh.createGeometry()                                                    # create geometry out of reparameterized surface
                gmsh.option.setNumber ("Mesh.MeshSizeMin", self.payloadMeshSize)                    # set mesh size
                gmsh.option.setNumber ("Mesh.MeshSizeMax", self.payloadMeshSize)                    # set mesh size
                gmsh.model.mesh.generate (dimSurf)                                                  # create new mesh
                
                # Remove duplicate nodes
                gmsh.model.mesh.removeDuplicateNodes()
                gmsh.model.mesh.removeDuplicateElements()
                
                # Get nodes for mesh
                #   note gmsh returns lists nodeTags=[1, 2, 3, ...], coord=[x1, y1, z1, x2, y2, z2, ...], and parametricCoord (not used)
                nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
                points = np.array (coord.reshape ((int(len(coord)/3), 3)))                          # reshape list to be 2D array with columns for xyz
                nodeMaxPayload = len (points)
                
                # Get elements for mesh
                elementTypes, elementTags, elements = gmsh.model.mesh.getElements (dimSurf)         # get all elements for remeshed payload surface
                elements = np.hstack (elements) - 1                                                 # horz stack array to remove array inside array and make zero offset for python
                elements = np.array (elements)                                                      # convert to numpy array
                elements = elements.astype(int)                                                     # convert to type int
                elements = elements.reshape ((int(len(elements)/(dimSurf+1)), dimSurf+1))           # reshape list to be 2D array with columns for node ids
               # gmsh.fltk.run()
                
            
            # Set offsetCG, add to mesh, and get max number of nodes in payload mesh
            offsetCG = 1
            points = np.append (np.array ([payloadCG]), points, axis=0)   
            nodeMaxPayload = len (points)
            print (">> >> >> >> Adding Payload: nodeMax = {}".format (nodeMaxPayload))
                
            # Save IDs of bridle attach points for structure mesh - assume they will be at z=0 slice per parachute mesh creation assumptions above
            nodeID_bridleAttach = []
            for n in range (len (data["node"])):
                if abs (data["node"][n, 2]) < EPS:
                    nodeID_bridleAttach.append (n + nodeMaxPayload)             # zero-offset node ID in original parachute mesh + offset for payload nodes
            bridleAttachMax = len (nodeID_bridleAttach)                         # number of attach points
            
            # Insert payload nodes to front of meshGen dataset
            data["node"] = np.append (points, data["node"], axis=0)
            
            # Prepend zeros for node displacement/rotation data for payload nodes
            data["displac"] = np.append (np.zeros ((nodeMaxPayload, 3)), data["displac"], axis=0)
            data["rotation"] = np.append (np.zeros ((nodeMaxPayload, 3)), data["rotation"], axis=0)
            
            # Offset all parachute cell_set element node ids to accound for payload nodes
            for setName in data["cell_set"].keys():
                data["cell_set"][setName]["elements"] = np.array (data["cell_set"][setName]["elements"])                # convert fluid elements from numpy format if not already
                data["cell_set"][setName]["elements"] += nodeMaxPayload                                                 # offset elements
                
                
            # PHANTOM SURFACE
            
            # Create surface
            setName = name + "_shell"           # for strucutre mesh, add shell to denote shell vs. beams
            setName = name                      # else if fluid, just keep Payload name
            elementType = 129                       # 3d 3-node plane stress triangular membrane with 3 dof/node
            self.attributes[setName] = {"matAttribute": -1, "elementType": elementType, "elementName": "PHANTOM_SHELL"} 
            data["cell_set"][setName] = {}
            data["cell_set"][setName]["type"] = "triangle"
            
            # Add payload surface elements - offset if necessary to account for CG as first node
            data["cell_set"][setName]["elements"] = elements + offsetCG


            # RIGID BEAMS
            
            # Set parameters
            setName = name + "_rigidBeam"
            elementType = 106                       # 3d 2-node rigid beam of arbitrary length, including zero length
            self.attributes[setName] = {"matAttribute": 100, "elementType": elementType, "elementName": "RIGID_BEAMS"} 
            data["cell_set"][setName] = {}
            data["cell_set"][setName]["type"] = "line"
                  
            # Create empty array for elements - offset by -1 to skip first (CG) point, and add number of bridle attach points
            elementMax = len (points) - offsetCG + bridleAttachMax
            data["cell_set"][setName]["elements"] = np.zeros((elementMax, 2)).astype(int)
            
            # Iterate over each shell node and bridle attach nodes to create beams to CG - skip CG point
            i = 0
            for e in range (elementMax):
                if e < len (points) - offsetCG:
                    data["cell_set"][setName]["elements"][e, 0] = e + offsetCG              # current payload surface node
                    data["cell_set"][setName]["elements"][e, 1] = 0                         # CG = first node
                else:
                    data["cell_set"][setName]["elements"][e, 0] = nodeID_bridleAttach[i]    # current bridle attach node
                    data["cell_set"][setName]["elements"][e, 1] = 0                         # CG = first node
                    i = i + 1
                    
                    
            # ADD ADDITIONAL BEAMS TO DISTRIBUTE BRIDLE LOAD
                    
            # Add additional beams to connect triple bridle attach points to other nodes and minimize penalty constraint issues
            for i in range (bridleAttachMax):
                
                # Attach bridle points to each other
                if i == 0:
                    new_element = np.array ([[nodeID_bridleAttach[-1], nodeID_bridleAttach[i]]])
                else:
                    new_element = np.array ([[nodeID_bridleAttach[i-1], nodeID_bridleAttach[i]]])
                data["cell_set"][setName]["elements"] = np.append (data["cell_set"][setName]["elements"], new_element, axis=0)

                # Find nodes of closet points to create additional support beams
                nSupports = self.nBridleSupportBeams                                # number of beams to connect each attach point to help distribute load
                pointAttach = data["node"][nodeID_bridleAttach[i]]                  # node coord of current attach point
                vectors = points - pointAttach                                      # get vectors from attach point to each payload node
                vectorsMag = np.linalg.norm (vectors, axis=1)                       # get magnitude of each vector
                nodeSupports = np.argsort(vectorsMag)[:nSupports]                   # get indices of closest nodes to current attach point

                # Add additional support beams to connect trip bridle attach points to closest points
                for n in nodeSupports:
                    new_element = np.array ([[nodeID_bridleAttach[i], n]])
                    data["cell_set"][setName]["elements"] = np.append (data["cell_set"][setName]["elements"], new_element, axis=0)    
                    
                # # Create support beams between bridle attach points and each surface node of payload to help distribute load (skip CG point as already beam already created)
                # # NOTE this version is *all* nodes so creates 4x the number of rigid beams... will influence the mass calculation
                # for n in range (1, len (points)):
                    # new_element = np.array ([[nodeID_bridleAttach[i], n]])
                    # data["cell_set"][setName]["elements"] = np.append (data["cell_set"][setName]["elements"], new_element, axis=0) 
                    


        # Save data back to class
        self.data[meshType] = data

                
        

    def parseInputExcel (self, outputFilename):
        """ Parse excel worksheet with parachute parameters
        """
        try:
            import openpyxl
        except:
            print ("ERROR: Required python dependency not available. Please install using the following command:")
            print ("   pip install {}".format ("openpyxl"))
            exit (0)
            
        # Settings
        maxHeaderRows = 20          # max number of rows to search for header and start of data
        maxHeaderCols = 20          # max number of cols to search for header and start of data
        maxDataRows = 1000          # max number of rows with data
        flagWriteSettings = True    # Write settings to text file?
        inputFilename = self.excelInputFilename
        inputSheetname = self.excelInputSheetname
                   
        # Open xlsx file and set worksheet
        wb = openpyxl.load_workbook (filename=inputFilename, data_only=True)
        ws = wb[inputSheetname]
        
        # Find row and columns for headers
        flagFoundItem = flagFoundVariable = flagFoundValue = False
        for colTmp in range (1, maxHeaderCols+1):
            for rowTmp in range (1, maxHeaderRows+1):
                valTmp = ws.cell (row=rowTmp, column=colTmp).value
                if flagFoundItem == False:
                    if (str (valTmp)).find("Item") >= 0:
                        rowHeader = rowTmp
                        colItem = colTmp
                        flagFoundItem = True
                if flagFoundVariable == False:
                    if (str (valTmp)).find("Variable") >= 0:
                        colVariable = colTmp
                        flagFoundVariable = True
                if flagFoundValue == False:
                    if (str (valTmp)).find("Value") >= 0:
                        colValue = colTmp
                        flagFoundValue = True
        if flagFoundItem != True  or  flagFoundVariable != True  or  flagFoundValue != True:
            raise Exception ("Error in parseInputExcel: header data not found...")
        
        # Read data into dictionaries for items
        for rowTmp in range (rowHeader, maxDataRows+1):
            
            # Save current row data to tmp variables
            itemTmp = ws.cell (row=rowTmp, column=colItem).value
            variableTmp = ws.cell (row=rowTmp, column=colVariable).value
            valueTmp = ws.cell (row=rowTmp, column=colValue).value
            
            # Pull data from rows with item name defined
            if itemTmp != None:
            
                # if current item does not exist, add to parachute parameter dict and set equal to blank dict to store data
                if itemTmp not in self.param:
                    self.param[itemTmp] = {}
                
                # If current row has variable name defined, save data to parameter dict
                if variableTmp != None:
                    self.param[itemTmp][variableTmp] = valueTmp           
        
        # Write setting dictionaries to text file to provide check if needed
        if flagWriteSettings == True:
            import pprint
            f = open (outputFilename, "w")
            pprint.pprint (self.param, stream=f, sort_dicts=False)
            f.write ("\n\n")
            f.close ()

       


    def createMeshRadials (self, name, ds, pStart, pEnd, factorSpacing, radius, iMax, meshType):
        """ Create Bridle lines, Suspension lines, Radials, etc. - i.e. lines aligned with z-axis
                name = name of component
                ds = mesh spacing [m]
                pStart = [x,y,z] start of line
                pEnd = [x,y,z] end of line
                factorSpacing = non-dimensional fractin of line radius to space out lines near apex and avoid overlapping surfaces
                radius = radius of dressing (for embedded lines in fluid simulation)
                iMax = number of lines
                meshType = fluid or structure
        """
        import numpy as np
        from math import pi
                
        # Create pointer to body parameters; create empty lists for meshes
        var = self.param
        print (">> >> >> Creating Geometry: {}".format (name))

        
        # STRUCTURE
        if meshType == "structure":
                          
            # Set parameters for line
            r0 = np.sqrt (pStart[0]**2.0 + pStart[1]**2.0)
            r1 = np.sqrt (pEnd[0]**2.0 + pEnd[1]**2.0)
            z0 = pStart[2]
            z1 = pEnd[2]
            
            # Iterate over segments
            points = []
            lines = []
            for i in range (iMax):            
                theta = 2.0 * pi * (i / iMax)
                p0 = [r0*np.cos(theta), r0*np.sin(theta), z0]
                p1 = [r1*np.cos(theta), r1*np.sin(theta), z1]
                dsTmp = min (ds, (np.linalg.norm (np.array(p1) - np.array(p0))) / self.minPointsEdge)
                points.append (gmsh.model.occ.addPoint (*p0, meshSize=dsTmp))
                points.append (gmsh.model.occ.addPoint (*p1, meshSize=dsTmp))
                lines.append  (gmsh.model.occ.addLine (points[2*i], points[2*i+1]))

            # Assign physical group to lines and synchronize
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup (dimLine, lines, name="{}".format(name))

            
        
        # FLUID
        if meshType == "fluid":
                    
            # Set parameters for line
            r0 = np.sqrt (pStart[0]**2.0 + pStart[1]**2.0)
            r1 = np.sqrt (pEnd[0]**2.0 + pEnd[1]**2.0)
            z0 = pStart[2]
            z1 = pEnd[2]  
            
            # Adjust ends of line to offset if needed and avoid overlapping another line
            if iMax > 1:
            
                # Adjust z to avoid overlapping line with another body/surface due to angle of line
                dr = r1 - r0
                dz = z1 - z0
                if abs (dr) >= EPS  and  abs (dz) >= EPS:       # only adjust if not horizontal or vertical line
                    theta = np.arctan (dr / dz)                 # angle of line
                    dz = abs (radius * np.sin (theta))          # amount of overlap in z-dir
                    dr = dz * np.tan (theta)                    # r pos where line will not overlap
                    r0 = r0 + dr
                    r1 = r1 - dr
                    z0 = z0 + dz
                    z1 = z1 - dz
                    
                # Adjust z based on cutback for line dressing at apex and not overlapping surfaces - only for multiple line cases
                theta = 2.0 * pi / iMax                                                 # angle between radial location of each circle
                ringRadius = (1.0 + factorSpacing) * radius / np.sin (theta / 2.0)      # calculate radial location for non-overlapping circle size
                if ringRadius > r0 and ringRadius > r1:
                    raise Exception ("Error in createMeshRadials / Fluid: cannot create dressing due to overlapping surfaces... Exiting")
                elif ringRadius > r0:
                    z0 = z0 + (z1 - z0) * ringRadius / r1       # calc z location based on angle of line
                    r0 = ringRadius
                elif ringRadius > r1:
                    z1 = z1 + (z0 - z1) * ringRadius / r0       # calc z location based on angle of line
                    r1 = ringRadius
                    
            # Create disk at start of line, extrude into cylinder
            surfaceGroup = []
            volumeGroup = []
            for i in range (iMax):
                
                surface = []
                volume = []
                # Define start/end points of current line center axis
                theta = 2.0 * pi * (i / iMax)
                p0 = [r0*np.cos(theta), r0*np.sin(theta), z0]
                p1 = [r1*np.cos(theta), r1*np.sin(theta), z1]
                
                # Set spacing on 
                dsTmp = min (ds, (2.0*pi*radius)/4)
                
                # Create face of parachute line -- NOTE THIS WORKS!!!
                surface.append (gmsh.model.occ.addDisk (*p0, rx=radius, ry=radius))

                # Rotate disk to align with line if not parallel or perpendicular to z-axis - only expect this for >1 line cases
                if iMax > 1:
                
                    # Calculate line vectors and rotation axis
                    v = np.array(p1) - np.array(p0)         # vector along line
                    v1 = np.array([v[0], v[1], 0.0])        # radial component
                    v2 = np.array([0.0,  0.0,  v[2]])       # axial component
                    if np.linalg.norm (v2) < EPS:
                        v2 = np.array ([0.0, 0.0, 1.0])     # reset v to z-axis if line is perpendicular
                    v3 = np.cross (v1, v2)                  # axis of rotation
                    
                    # Set rotation angle
                    if abs (z0 - z1) < EPS:
                        rotAngle = pi/2                                                         # set angle to 90deg for z-perpendicular lines
                    else:
                        rotAngle = -np.arctan2 (np.linalg.norm (v1), np.linalg.norm (v2))       # else caclulate from line vector
                        
                    # Only apply non-zero rotation angle
                    if abs (rotAngle) > EPS:
                        gmsh.model.occ.rotate ([(dimSurf, surface[0])], *p0, *v3, rotAngle)
                
                # Extrude line
                length = np.linalg.norm (np.array(p1) - np.array(p0))
                numLayers = max (int (length/ds)+1, self.minPointsEdge)
                dx = p1[0] - p0[0]
                dy = p1[1] - p0[1]
                dz = p1[2] - p0[2]
                extruded = gmsh.model.occ.extrude ([(dimSurf, surface[0])], dx, dy, dz, [numLayers])    # returns dim,tag pairs for top surf, volume, and side surfaces (note bottom surface is original used as input to extrude)
                for (dim, tag) in extruded:
                    if dim == dimSurf:
                        surface.append (tag)
                    elif dim == dimVol:
                        volume.append (tag)           
                gmsh.model.occ.synchronize()
                
                # Save surfaces and volumes to label as physical group below
                surfaceGroup = surfaceGroup + surface
                volumeGroup = volumeGroup + volume
            
            # Assign physical groups to model
            gmsh.model.addPhysicalGroup (dimVol, volumeGroup, name="{}".format(name))
            gmsh.model.addPhysicalGroup (dimSurf, surfaceGroup, name="{}".format(name))

                
                ### TO DO LATER
                # Create cone for cutback portion if multiple lines
                # if iMax > 1:
                    # if abs (z0 - pStart[2]) > EPS:
                        # pBase = p0
                        # pTip = pStart
                    # else:
                        # pBase = p1
                        # pTip = pEnd
                    # model.characteristic_length_max = dsTmp
                    # volume.append (model.add_cone (pBase, (np.array(pTip)-np.array(pBase)), radius, 0))
        


        



    def createMeshPanels (self, name, radialName, circumName, ds, rLoc, zLoc, iMax, meshType):
        """ Create circumferential reiforcement lines mesh - i.e. segmented lines that are circumferential around z-axis in single z-plane
                name = name of component
                radialName = name of set of radials (sides of panel)
                circumName = name of set of circumferentials (LE/TE of panel
                ds = nominal mesh spacing [m]
                rLoc = 2-element list of radial locations (relative to z-axis) of line segment endpoints
                zLoc = 2-element list of z locations of line
                iMax = number of segments (e.g. gores)
                meshType = fluid or structure
        """
        import numpy as np
        from math import pi
                
        print (">> >> >> Creating Geometry: {}".format (name))
                
        # Iterate over segments
        linesLE = []
        linesTE = []
        linesRadials = []
        for i in range (iMax):      
            theta0 = 2.0 * pi * (i / iMax)
            theta1 = 2.0 * pi * ((i + 1) / iMax)
            
            # Calculate points of polygon surface
            p0 = [rLoc[0]*np.cos(theta0), rLoc[0]*np.sin(theta0), zLoc[0]]
            p1 = [rLoc[0]*np.cos(theta1), rLoc[0]*np.sin(theta1), zLoc[0]]
            p2 = [rLoc[1]*np.cos(theta1), rLoc[1]*np.sin(theta1), zLoc[1]]
            p3 = [rLoc[1]*np.cos(theta0), rLoc[1]*np.sin(theta0), zLoc[1]]
            
            # Calculate minimum spacing
            ds0 = min (ds, (np.linalg.norm (np.array(p1) - np.array(p0)) / self.minPointsEdge))
            ds1 = min (ds, (np.linalg.norm (np.array(p2) - np.array(p1)) / self.minPointsEdge))
            ds2 = min (ds, (np.linalg.norm (np.array(p3) - np.array(p2)) / self.minPointsEdge))
            ds3 = min (ds, (np.linalg.norm (np.array(p0) - np.array(p3)) / self.minPointsEdge))
            
            # Create points and lines for mesh
            points = []
            points.append (gmsh.model.occ.addPoint (*p0, meshSize=0.5 * (ds3 + ds0)))
            points.append (gmsh.model.occ.addPoint (*p1, meshSize=0.5 * (ds0 + ds1)))
            points.append (gmsh.model.occ.addPoint (*p2, meshSize=0.5 * (ds1 + ds2)))
            points.append (gmsh.model.occ.addPoint (*p3, meshSize=0.5 * (ds2 + ds3)))
            gmsh.model.occ.synchronize()            
            lines = []
            lines.append (gmsh.model.occ.addLine (points[0], points[1]))
            lines.append (gmsh.model.occ.addLine (points[1], points[2]))
            lines.append (gmsh.model.occ.addLine (points[2], points[3]))
            lines.append (gmsh.model.occ.addLine (points[3], points[0]))
            gmsh.model.occ.synchronize()
                         
            # Create curve loop and surface
            curves = []
            curves.append (gmsh.model.occ.addCurveLoop (lines))
            gmsh.model.occ.synchronize()
            surface = gmsh.model.occ.addPlaneSurface (curves)
            gmsh.model.occ.synchronize()
                
            # Append lines to lists
            linesLE.append      (lines[0])
            linesRadials.append (lines[1])
            linesTE.append      (lines[2])
            linesRadials.append (lines[3])
                
            # Assign physical group to surfaces
            gmsh.model.addPhysicalGroup (dimSurf, [surface], name="{}_{}".format(name, i))
            
            # Save theta [rad] bounds of current gore to attributes
            self.attributes[name+"_{}".format(i)]["thetaBounds"] = [theta0, theta1]
                
        # Assign physical group to lines for structure mesh
        if meshType == "structure":
            if radialName != None:
                gmsh.model.addPhysicalGroup (dimLine, linesRadials, name=radialName)
            if circumName != None:
                gmsh.model.addPhysicalGroup (dimLine, linesLE, name=circumName[0])
                gmsh.model.addPhysicalGroup (dimLine, linesTE, name=circumName[1])
                        
        # Synchronize model
        gmsh.model.occ.synchronize()




    def calcElementFrames (self, meshType):
        """ Calculate locale coordinate frames for structure material props
            - Will create s = <s1, s2, s3> which represents primary, secondary, and tertiary axis for material laws
              in default coordinate system (i.e. z-dir is flow direction)
            - Assumes right-hand convention
            - For 1D elements: s1 in line-direction, s2 and s3 are calculated based on normals to main coordinate system
            - For 2D elements: s1 is primary direction, s2 is rotated 90deg in plane, and s3 is cross product of s1 and s2
            - Note for panels, material directions may be rotated based on attributes (i.e. block vs. bias construction)
        """
        import numpy as np
        
        # Make sure this is only called on structure mesh
        if meshType != "structure":
            raise Exception ("Error in calcElementFrames: can only be called on 'structure' meshType.")
        
        # Gore naming label - used to allow alternating directions if flagged
        goreLabel = "Gores_"
        
        # Copy data to local arrays
        var = self.param
        data = self.data["structure"]
        node = np.array (data["node"])
 
        # Init arrays for points and vectors
        p0 = np.zeros ((3))
        p1 = np.zeros ((3))
        p2 = np.zeros ((3))
        p3 = np.zeros ((3))
        v1 = np.zeros ((3))
        v2 = np.zeros ((3))
        pCenter = np.zeros ((3))
        
        # Iterate over cell sets
        print (">> >> >> Calculating material property coordinate frames...")
        for setName, setData in data["cell_set"].items():
            print (">> >> >> >> Calculating frames... setName = {}".format (setName))
            
            # Save elements to local array, and create array to store coordinate frames
            setElement = np.array (setData["elements"])
            setElementMax = len (setElement)
            frameList = []
            
            # Get name of var block based on current set name
            for varName, varData in var.items():
                if setName.find(varName) >= 0:
                    name = varName
                    break
                    
            # PAYLOAD ELEMENTS
            if setName.find ("Payload") >= 0:
                
                # Iterate over elements
                for e in range (setElementMax):

                    # For payload elements, just set frames to x, y, z axes as they are not needed
                    frameList.append({"s0": np.array ([1, 0, 0]),
                                      "s1": np.array ([0, 1, 0]),
                                      "s2": np.array ([0, 0, 1])})
                
          
            # PARACHUTE 1D ELEMENTS
            elif setData["type"] == "line":
            
                # Iterate over elements
                for e in range (setElementMax):
                    
                    # Get x,y,z data for each node of element
                    for xyz in range (3):
                        p0[xyz] = node[setElement[e, 0], xyz]
                        p1[xyz] = node[setElement[e, 1], xyz]
                        
                    # Calculate dx, dy, dz for line
                    dx = p1[0] - p0[0]
                    dy = p1[1] - p0[1]
                    dz = p1[2] - p0[2]
                    
                    # Set frame axes s0, s1, s2
                    frameTmp = {}
                    frameTmp["s0"] = unit_vector (np.array([dx, dy, dz]))                               # set s0 to line direction
                    if abs (dx) < EPS and abs (dy) < EPS:
                        frameTmp["s1"] = unit_vector (np.array([1, 0, 0]))                              # if along z-axis, set s1 to x axis
                    else:
                        frameTmp["s1"] = unit_vector (np.cross (frameTmp["s0"], np.array([0, 0, 1])))   # else set s1 to cross product of line and z-axis
                    frameTmp["s2"] = unit_vector (np.cross (frameTmp["s0"], frameTmp["s1"]))            # set s2 as cross product of s0 and s1
                    frameList.append (frameTmp)                                                         # append frame to list
               
            
            # PARACHUTE 2D ELEMENTS
            elif setData["type"] == "triangle"  or  setData["type"] == "quad":
                
                # Determine radial angle for current gore (i.e. angle projected from z-axis along center of gore)
                theta = np.average (self.attributes[setName]["thetaBounds"])
                
                # Add theta offset if needed if gore panel and alternate attribute is true
                if setName.find (goreLabel) >= 0 and "flagAlternateGores" in self.attributes[setName]:
                    id = int (setName.split(goreLabel)[1])                                              # get gore ID
                    if id % 2 != 0 and self.attributes[setName]["flagAlternateGores"] == True:
                        theta = theta - np.pi / 2.0                                                     # if odd id, set theta offset to -90deg to reverse S0 and S1
                
                # Adjust theta for weave angle (e.g. block vs. bias construction)
                theta = theta + np.radians (var[name]["theta"])
                    
                # Iterate over elements
                for e in range (setElementMax):
                    
                    if setData["type"] == "triangle":
                        for xyz in range (3):                               # Get x,y,z data for each node of element
                            p0[xyz] = node[setElement[e, 0], xyz]
                            p1[xyz] = node[setElement[e, 1], xyz]
                            p2[xyz] = node[setElement[e, 2], xyz]
                            pCenter[xyz] = (p0[xyz] + p1[xyz] + p2[xyz]) / 3.0
                        for xyz in range (3):                               # Calculate vectors on edges of triangle
                            v1[xyz] = p1[xyz] - p0[xyz]
                            v2[xyz] = p2[xyz] - p0[xyz]
                            
                    elif setData["type"] == "quad":
                        for xyz in range (3):                               # Get x,y,z data for each node of element
                            p0[xyz] = node[setElement[e, 0], xyz]
                            p1[xyz] = node[setElement[e, 1], xyz]
                            p2[xyz] = node[setElement[e, 2], xyz]
                            p3[xyz] = node[setElement[e, 3], xyz]
                            pCenter[xyz] = (p0[xyz] + p1[xyz] + p2[xyz] + p3[xyz]) / 4.0
                        for xyz in range (3):                               # Calculate vectors on edges of triangle
                            v1[xyz] = p2[xyz] - p0[xyz]
                            v2[xyz] = p3[xyz] - p1[xyz]    
                        
                    # Calc surface normal and orient so it points outward (based on vector from origin)
                    vNormal = np.array (unit_vector (np.cross (v1, v2)))
                    if np.dot (vNormal, unit_vector (np.array (pCenter))) < 0.0:
                        vNormal = -vNormal
                    
                    # Calc local vectors for direction based on theta (which defines warp/weft) and surface normal
                    frameTmp = {}
                    frameTmp["s0"] = unit_vector (np.array ([np.cos(theta), np.sin(theta), 0]))                     # set s0 to align with gore + weave offset angle
                    frameTmp["s1"] = unit_vector (np.array ([np.cos(theta+np.pi/2), np.sin(theta+np.pi/2), 0]))     # set s1 to be rotated 90 deg from s0
                    frameTmp["s2"] = unit_vector (np.array ([0, 0, 1]) )                                            # set s3 to be z-axis (i.e. normal to local 2D frame)
                    
                    # Create rotation matrix transform local frame to global frame
                    R = createRotationMatrix (frameTmp["s2"], vNormal)
                    for nameFrame, vector in frameTmp.items():
                        frameTmp[nameFrame] = np.dot (R, vector)                                                    # Convert local frame vectors to global frame
                    frameList.append (frameTmp)                                                                     # append frame to list            
                
                
            # If unrecognized element type raise exception
            else:
                raise Exception ("calcElementFrames: unrecognized element type")
            
            # Copy frame data back to cell set
            self.data["structure"]["cell_set"][setName]["frames"] = frameList
 
            
            
        
        
    def writeMeshAeroF (self, filename, flagCombineSurfs=False):
        """ Write ascii file in XPost (TOP) format with fluid mesh data for AERO-F
            - Assumes mesh is meshGen data object
            - Assumes parachute mesh tree - i.e. multiple surfaces stored in a dictionary, with "fluid" = embedded surfaces
        """
        import numpy as np
        
        
        # COMBINE SURFACES OR COPY DATA
        flagSeparateBandGores = False
        
        # If outputting all surfaces, copy data to local variable
        if flagCombineSurfs == False:
            data = self.data["fluid"]
        
        # Else copy gores to just have one surface for disk and one for band
        elif flagSeparateBandGores == True:
            
            # Init data dict
            data = {}
            data["node"] = self.data["fluid"]["node"]
            data["cell_set"] = {}
            data["cell_set"]["DiskGores"] = {"elements": [], "type": "triangle"}
            data["cell_set"]["BandGores"] = {"elements": [], "type": "triangle"}
            
            # Copy disk and band gores to one cell_set for each, otherwise just copy cell set
            for setName, setData in self.data["fluid"]["cell_set"].items():
                if (setName.find ("DiskGores")) >= 0:
                    data["cell_set"]["DiskGores"]["elements"] = data["cell_set"]["DiskGores"]["elements"] + setData["elements"]
                elif (setName.find ("BandGores")) >= 0: 
                    data["cell_set"]["BandGores"]["elements"] = data["cell_set"]["BandGores"]["elements"] + setData["elements"]
                else:
                    data["cell_set"][setName] = setData

        # Else copy gores to just have one surface for both disk and band
        else:
            
            # Init data dict
            data = {}
            data["node"] = self.data["fluid"]["node"]
            data["cell_set"] = {}
            data["cell_set"]["Canopy"] = {"elements": [], "type": "triangle"}
            
            # Copy disk and band gores to one cell_set, otherwise just copy cell set
            for setName, setData in self.data["fluid"]["cell_set"].items():
                if (setName.find ("Gores")) >= 0:
                    data["cell_set"]["Canopy"]["elements"] = data["cell_set"]["Canopy"]["elements"] + setData["elements"]
                else:
                    data["cell_set"][setName] = setData
                    
                          
        # REORDER SET DATA
        
        # Get ordered list of keys with Payload ones first
        setOrderList = []
        for setName in data["cell_set"].keys():
            if setName.find ("Payload") >= 0:
                setOrderList.append (setName)               # add Payload dict keys first
        for setName in data["cell_set"].keys():
            if setName.find ("Payload") < 0:
                setOrderList.append (setName)               # then add non-Payload dict keys 
                
        
        # WRITE DATA
        
        # Open output file
        print (">> >> Writing Fluid Mesh Data... {}".format (filename+".top"))
        f = open (filename+".top", "w")   
        
        # NODES
        # Write nodes for all surfaces
        points = np.array (data["node"])
        nodeMax = len (points)
        print (">> >> >> nodeMax = {}".format (nodeMax))            

        # Write nodes
        f.write ("Nodes EmbeddedNodes\n")   
        for n in range (nodeMax):
            f.write ("{:<12d}      {:>25.16e} {:>25.16e} {:>25.16e}\n".format (n+1, points[n, 0], points[n, 1], points[n, 2]))
        
        # ELEMENTS
        # Write elements for all surfaces
        surfCount = 1
        elementCount = 1
        
        # Iterate over cell sets - skip line cell sets such as radials
        for setName in setOrderList:
            setData = data["cell_set"][setName]
            if setData["type"] == "triangle":
        
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                print (">> >> >> Writing Fluid Mesh: setName = {:24s}  elementMax = {:<12d}".format (setName, setElementMax))
                
                # Write elements
                elementType = 4     # triangle = 4
                f.write ("Elements {}_{}_{} using EmbeddedNodes\n".format (self.BCtype, setName, surfCount))
                for e in range (setElementMax):
                    f.write ("{:<12d}   {:<6d}       {:>12d} {:>12d} {:>12d}\n".format (elementCount, elementType, setElement[e, 0]+1, setElement[e, 1]+1, setElement[e, 2]+1))
                    elementCount = elementCount + 1

                # Increment counters
                surfCount = surfCount + 1
            
        # Close output file
        f.close()




    def writeMeshAeroS (self, filename, flagPayloadOnly=False):
        """ Write ascii file in FEM format with structure mesh data for AERO-S
            - Assumes mesh is meshGen data object
            - Assumes parachute mesh tree - i.e. multiple surfaces stored in a dictionary
            - Option to write just payload mesh only (nodes_topology and attributes files only)
            - Inputs:
                filename                        prefix to output filenames
                flagPayloadOnly                 True: write payload mesh only, False: write all entities
            - Writes AERO-S required blocks: 
                filename.fem:                   NODES, TOPOLOGY
                filename_attributes.fem:        ATTRIBUTES
                filename_frames.fem:            EFRAMES, CFRAMES
                filename_surfacetopo.fem        SURFACETOPO
            - Writes AERO-S decomposition file for parallel analysis:
                filename.dec
        """
        import numpy as np
        
        # Save mesh data to local dict
        data = self.data["structure"]
        
        # Other settings
        flagWriteEFrames = False
        
        #####################
        # NODES_TOPOLOGY FILE        
        # Open output file
        print (">> >> >> Writing Structure Mesh Data... NODES TOPOLOGY...")
        f = open (filename+".fem", "w")

        # NODES
        
        # Write nodes for all surfaces 
        points = np.array (data["node"])
        nodeMax = len (points)         

        # Write nodes
        f.write ("NODES\n") 
        for n in range (nodeMax):
            f.write ("{:<12d}      {:>25.16e} {:>25.16e} {:>25.16e}\n".format (n+1, points[n, 0], points[n, 1], points[n, 2]))
        
        # Init dict to store element IDs - note these will be unit offset as required by AERO suite
        setElementID = {}
        
        # Iterate over cell sets
        elementCount = 1
        for setName, setData in data["cell_set"].items():
            if flagPayloadOnly == False  or  setName.find ("Payload") >= 0:
        
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                setType = setData["type"]
                
                # Init list to store element IDs for current set
                setElementID[setName] = []

                # Get element name and type
                elementName = self.attributes[setName]["elementName"]
                elementType = self.attributes[setName]["elementType"]
                
                # Write Header
                f.write ("*\n")
                f.write ("* =============================================================================\n")
                f.write ("* mesh set\n")
                f.write ("* -----\n")
                f.write ("* name:         {}\n".format (setName))
                f.write ("* elementType:  {} {}\n".format (elementName, elementType))
                f.write ("* elementMax:   {}\n".format (setElementMax))
                f.write ("* meshType:     {}\n".format (setType))
                f.write ("*\n")
                
                # TOPOLOGY
                # Write elements - make sure to associate with global node number
                f.write ("TOPOLOGY\n")       
                if setType == "line":
                    for e in range (setElementMax):
                        #         element#  type       node1   node2     
                        f.write ("{:<16d}   {:<6d}     {:>16d} {:>16d}\n".format (elementCount+e, elementType, setElement[e, 0]+1, setElement[e, 1]+1))   
                        setElementID[setName].append (elementCount+e)        # save element ID
                elif setType == "triangle":
                    for e in range (setElementMax):
                        #         element#  type       node1   node2   node3
                        f.write ("{:<16d}   {:<6d}     {:>16d} {:>16d} {:>16d}\n".format (elementCount+e, elementType, setElement[e, 0]+1, setElement[e, 1]+1, setElement[e, 2]+1))
                        setElementID[setName].append (elementCount+e)        # save element ID
                elif setType == "quad":
                    for e in range (setElementMax):
                        #         element#  type       node1   node2   node3
                        f.write ("{:<16d}   {:<6d}     {:>16d} {:>16d} {:>16d} {:>16d}\n".format (elementCount+e, elementType, setElement[e, 0]+1, setElement[e, 1]+1, setElement[e, 2]+1, setElement[e, 3]+1))
                        setElementID[setName].append (elementCount+e)        # save element ID
                print (">> >> >> >> Creating Topology setName={}, numElements={}".format(setName, len(setElementID[setName])))
                
                # Iterate elementCount and id
                elementCount = elementCount + setElementMax
            
        # Close output file
        f.close()


        #################
        # ATTRIBUTES / MATUSAGE FILE        
        # Open output file
        print (">> >> >> Writing Structure Mesh Data... ATTRIBUTES MATUSAGE...")
        f = open (filename+"_attributes.fem", "w")

        # Iterate over cell sets
        for setName, setData in data["cell_set"].items():
            if flagPayloadOnly == False  or  setName.find ("Payload") >= 0:
        
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                setType = setData["type"]

                # Get element name and type
                elementName = self.attributes[setName]["elementName"]
                elementType = self.attributes[setName]["elementType"]            
                
                # Get material attribute number - note this will need to map to AERO-S input deck
                setMatAttribute = self.attributes[setName]["matAttribute"]

                # Write Header
                f.write ("*\n")
                f.write ("* =============================================================================\n")
                f.write ("* attributes set\n")
                f.write ("* -----\n")
                f.write ("* name:         {}\n".format (setName))
                f.write ("* elementType:  {} {}\n".format (elementName, elementType))
                f.write ("* elementMax:   {}\n".format (setElementMax))
                f.write ("* meshType:     {}\n".format (setType))
                f.write ("* matAttribute: {}\n".format (setMatAttribute))
                f.write ("*\n")
                
                                    
                # PAYLOAD: Write payload rigid beams and phantom shells - MATUSAGE not needed as no MATLAW specified for these
                if setName.find ("Payload") >= 0:
                    f.write ("ATTRIBUTES\n")
                    f.write ("{} {}     {}\n".format (setElementID[setName][0], setElementID[setName][-1], setMatAttribute))    # start element, end element, materialID
                    
                # OTHER: Write other rigid beams and shell elements - MATUSAGE requlired if using MATLAW in input deck
                else:
                    # Write data for beam elements 
                    if setType == "line":
                        f.write ("ATTRIBUTES\n")
                        f.write ("{} {}     {}\n".format (setElementID[setName][0], setElementID[setName][-1], setMatAttribute))        # start element, end element, materialID
                        f.write ("MATUSAGE\n")
                        f.write ("{} {}     {}\n".format (setElementID[setName][0], setElementID[setName][-1], setMatAttribute))        # start element, end element, matlawID
                        
                    # Write data for shell elements
                    elif setType == "triangle"  or  setType == "quad":
                        f.write ("ATTRIBUTES\n")
                        for e in range (setElementMax):
                            f.write ("{}      {}  {}  {}\n".format (setElementID[setName][e], setMatAttribute, -1, setElementID[setName][e]))   # element ID, end element, materialID, composite attribute, cframe = elementID
                        f.write ("MATUSAGE\n")
                        f.write ("{} {}     {}\n".format (setElementID[setName][0], setElementID[setName][-1], setMatAttribute))                # start element, end element, matlawID
                
        # Close output file
        f.close()
        

        #############
        # IDISP6 FILE        
                
        # Open output file
        print (">> >> >> Writing Structure Mesh Data... IDISP6...")
        f = open (filename+"_idisp6.fem", "w")
        
        # Write header
        f.write ("IDISP6\n")

        # Write data for any nodes that have been transformed for folding
        for n in range (len (data["node"])):
            if np.linalg.norm (data["displac"][n, :]) > 0:
                f.write ("{}   {} {} {}   {} {} {}\n".format (n+1, data["displac"][n, 0], data["displac"][n, 1], data["displac"][n, 2], data["rotation"][n, 0], data["rotation"][n, 1], data["rotation"][n, 2]))
        f.close()
                
               
       ##############
        # FRAMES FILE        
        if flagWriteFrames == True  and  flagPayloadOnly == False:
                
            # Open output file
            print (">> >> >> Writing Structure Mesh Data... EFRAMES CFRAMES...")
            f = open (filename+"_frames.fem", "w")

            # Iterate over cell sets
            for setName, setData in data["cell_set"].items():
            
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                setType = setData["type"]
                setFrame = setData["frames"]

                # Write Header
                f.write ("*\n")
                f.write ("* =============================================================================\n")
                f.write ("* frames set\n")
                f.write ("* -----\n")
                f.write ("* name:         {}\n".format (setName))
                f.write ("* elementMax:   {}\n".format (setElementMax))
                f.write ("* meshType:     {}\n".format (setType))
                f.write ("*\n")    
                
                # EFRAMES or CFRAMES
                # Beam elements: write material prop coordinate frames in global coordinates - beam elements
                if setType == "line":
                    if flagWriteEFrames == True:
                        f.write ("EFRAMES\n")
                        for e in range (setElementMax):
                            #         element#  S1x        S1y        S1z          S2x        S2y        S2z          S3x        S3y        S3z
                            f.write ("{:<16d}   {:>25.16e} {:>25.16e} {:>25.16e}   {:>25.16e} {:>25.16e} {:>25.16e}   {:>25.16e} {:>25.16e} {:>25.16e}\n".format (setElementID[setName][e], setFrame[e]["s0"][0], setFrame[e]["s0"][1], setFrame[e]["s0"][2], setFrame[e]["s1"][0], setFrame[e]["s1"][1], setFrame[e]["s1"][2], setFrame[e]["s2"][0], setFrame[e]["s2"][1], setFrame[e]["s2"][2])) 
                    else:
                        f.write ("* ... Skipping EFRAMES\n")
                # Membrane elements: Write material prop coordinate frames in global coordinates - surface/membrane elements
                # Note we are setting a composite frame for each element, which is the most general case where composite does not need to be flat
                elif setType == "triangle"  or  setType == "quad":
                    f.write ("CFRAMES\n")
                    for e in range (setElementMax):
                        #         cmp_frm#  S1x        S1y        S1z          S2x        S2y        S2z          S3x        S3y        S3z
                        f.write ("{:<16d}   {:>25.16e} {:>25.16e} {:>25.16e}   {:>25.16e} {:>25.16e} {:>25.16e}   {:>25.16e} {:>25.16e} {:>25.16e}\n".format (setElementID[setName][e], setFrame[e]["s0"][0], setFrame[e]["s0"][1], setFrame[e]["s0"][2], setFrame[e]["s1"][0], setFrame[e]["s1"][1], setFrame[e]["s1"][2], setFrame[e]["s2"][0], setFrame[e]["s2"][1], setFrame[e]["s2"][2])) 
                
            # Close output file
            f.close()
       

        
        
        ####################
        # DECOMPOSITION FILE
        if flagWriteDecomp == True  and  flagPayloadOnly == False:
            print (">> >> >> Writing Structure Mesh Data... DECOMPOSITION")
            
            # Init dicts to store subdomains
            subElementID = []
            subName = []
            subCount = 0     
                
            # Save gore subdomains
            collectionNames = ["BandGores", "DiskGores"]
            for nameTmp in collectionNames:
                
                # Init vars
                goreCount = 0
                elementID = []
                
                # Set number of gores per subdomain
                if nameTmp.find ("Band") >= 0:
                    goresPerSubdomain = self.goresBandPerSubdomain
                else:
                    goresPerSubdomain = self.goresDiskPerSubdomain
                
                # Iterate over cell sets
                for setName, setData in data["cell_set"].items():
                    
                    # If set is in current collection, save subdomain data
                    if setName.find(nameTmp) >= 0:
                        goreCount = goreCount + 1                               # iterate gore count
                        elementID = elementID + setElementID[setName]           # append list of elementIDs for current gore

                        # Save subdomain data and reset gore count as needed
                        if goreCount == goresPerSubdomain:
                            subName.append ("{}_{}".format (subCount, nameTmp)) # append name to current subdomain
                            subElementID.append (elementID)                     # append current elementIDs to subdomain
                            goreCount = 0                                       # reset gore count for next subdomain
                            elementID = []                                      # reset element ID list for next subdomain
                            subCount = subCount + 1                             # iterate subdomain numbering for next subdomain
            
            # Reinit for other subdomains
            elementID = []
            
            # Iterate over cell sets
            for setName, setData in data["cell_set"].items():
                
                # Check if current set is "other" subdomain - set flag true if so
                flagOther = True
                for nameTmp in collectionNames:
                    if setName.find (nameTmp) >= 0:
                        flagOther = False
                 
                # If other subdomain, add data to list
                if flagOther == True:
                    elementID = elementID + setElementID[setName]
            
            # Split other items into apporpriate number of subdomains
            arrayTmp = np.array_split (np.array (elementID), self.numOtherSubdomain)
            for i in range (self.numOtherSubdomain):
                subElementID.append (arrayTmp[i].tolist())
                subName.append("{}_{}".format (subCount, "Other"))
                subCount = subCount + 1
            
            # Write decomposition
            f = open (filename+".dec", "w")
            f.write ("Decomposition dec for elemset\n")                                                             # write header            
            print (">> >> >> >> Created {} subdomains:".format (len (subName)))
            f.write (" {}\n".format(len(subName)))                                                                  # write number of subdomains
            for i in range (len (subName)):                                                                         # iterate over subdomains
                print(">> >> >> >> >> Subdomain {}: {} elements".format (subName[i], len(subElementID[i])))
                f.write(" {}\n".format (len (subElementID[i])))                                                     # write num elements in currend subdomain
                for e in range (len (subElementID[i])):                                                             # iterate over elements
                    f.write ("{}\n".format(subElementID[i][e]))                                                     # write current element ID
            f.close()
                
            # Get total number of elements in subdomains
            count = 0
            for i in range (len (subName)):
                count = count + len (subElementID[i])
            print (">> >> >> >> Total number of elements = {}".format (count))
                
        
        ##################
        # SURFACETOPO FILE
        if flagWriteSurfaceTopo == True  and  flagPayloadOnly == False:
            
            # Open output file
            print (">> >> >> Writing Structure Mesh Data... SURFACETOPO")
            f = open (filename+"_surfacetopo.fem", "w")

            # Init dicts to hold data collections
            collectionNames = ["BandGores", "DiskGores", "Canopy"]
            dataCollection = {}
            for nameTmp in collectionNames:
                dataCollection[nameTmp] = {"elements": []}
                
            # Iterate over cell sets to collect band and disk data into sets
            for setName, setData in data["cell_set"].items():
            
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                setType = setData["type"]
                
                # If surface with tri elements, check to save data
                for nameTmp in collectionNames:
                    
                    # Band and Disk Gores
                    if (setType == "triangle"  or  setType == "quad")  and  setName.find(nameTmp) >= 0:
                        dataCollection[nameTmp]["elements"] = dataCollection[nameTmp]["elements"] + setElement.tolist()
                        dataCollection[nameTmp]["type"] = setData["type"]
                        
            # Save Canopy data
            dataCollection["Canopy"]["elements"] = dataCollection["BandGores"]["elements"] + dataCollection["DiskGores"]["elements"]
            dataCollection["Canopy"]["type"] = dataCollection["BandGores"]["type"]
                         
            # Add dataCollection to end of existing data set
            for setName, setData in dataCollection.items():
                data["cell_set"][setName] = setData
                            
            # Iterate over cell sets to collect band and disk data into single ones
            var = self.param
            elementCount = 1
            id = 1              # start with id=1
            for setName, setData in data["cell_set"].items():
            
                # Get set info and print to screen
                setElement = np.array (setData["elements"])
                setElementMax = len (setElement)
                setType = setData["type"]
                
                # Only process surface data cell sets
                if setType == "triangle"  or  setType == "quad":
                                               
                    # Set thickness for Band, Disk, and Canopy
                    if setName.find ("BandGores") >= 0:
                        thickness = var["Band"]["t"]
                    elif setName.find ("DiskGores") >= 0:
                        thickness = var["Disk"]["t"]
                    elif setName.find ("Canopy") >= 0:
                        thickness = var["Disk"]["t"]        # assume band and disk have same thickness
                    else:
                        thickness = None                   
                    
                    # Write Header
                    f.write ("*\n")
                    f.write ("* =============================================================================\n")
                    f.write ("* surface topology set\n")
                    f.write ("* -----\n")
                    f.write ("* name:         {}\n".format (setName))
                    f.write ("* ID:           {}\n".format (id))
                    f.write ("* elementMax:   {}\n".format (setElementMax))
                    f.write ("* meshType:     {}\n".format (setType))
                    f.write ("*\n")

                    # SURFACETOPO
                    # Note only applies to surfaces, and assumes all surfaces have triangle elements (type "3")
                    if thickness != None:
                        f.write ("SURFACETOPO {} SURFACE_THICKNESS {}\n".format (id, thickness))
                    else:
                        f.write ("SURFACETOPO {}\n".format (id))
                    if setType == "triangle":
                        for e in range (setElementMax):
                            #         element#  facetype   node1   node2   node3
                            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d}\n".format (elementCount+e, 3, setElement[e, 0]+1, setElement[e, 1]+1, setElement[e, 2]+1)) 
                    elif setType == "quad":
                        for e in range (setElementMax):
                            #         element#  facetype   node1   node2   node3   node4
                            f.write ("{:<12d}   {:<6d}     {:>12d} {:>12d} {:>12d} {:>12d}\n".format (elementCount+e, 3, setElement[e, 0]+1, setElement[e, 1]+1, setElement[e, 2]+1, setElement[e, 3]+1)) 
                     
                    # Iterate elementCount and id
                    elementCount = elementCount + setElementMax
                    id = id + 1
                
            # Close output file
            f.close()

        
        

    def writeMeshTecplot (self, filename, meshType, flagWriteFrames=False):
        """ Write ascii file in Tecplot format with fluid or structure mesh data
            - Assumes mesh is meshGen data object
            - Assumes parachute mesh tree - i.e. multiple surfaces stored in a dictionary
            - Writes line and triangle data in unstructured format, as well as frame data (structure) if specified
        """
        import numpy as np
        
        # Save mesh data to local dict
        data = self.data[meshType]
        
        # Open output file
        f = open (filename+"_tecplot.dat", "w")

        # Write header
        f.write ("TITLE = MeshGen {} Mesh\n".format (meshType))
        if flagWriteFrames == True:
            f.write ('VARIABLES = "X", "Y", "Z", "S0_X", "S0_Y", "S0_Z", "S1_X", "S1_Y", "S1_Z", "S2_X", "S2_Y", "S2_Z"\n')
        else:
            f.write ('VARIABLES = "X", "Y", "Z"\n')
            
        # Copy node data to point array
        points = np.array (data["node"])
      
        
        # BLOCK SET DATA
                
        # Write all surfaces
        print (">> >> Writing Mesh Data in Tecplot format... {}".format (filename+"_tecplot.dat"))
        
        # Iterate over cell sets
        for setName, setData in data["cell_set"].items():
        
            # Get set info and print to screen
            setElement = np.array (setData["elements"])
            setElementMax = len (setElement)
            setType = setData["type"]

            # Init vars
            nodeList = []
            elementList = []     
            if flagWriteFrames == True: 
                try:
                    setFrame = setData["frames"]
                    frameList = []
                except:
                    import pdb
                    pdb.set_trace()
            # Create tmp node and element set for current block
            n = 0
            if setType == "line":
                typeLabel = "FELINESEG"
                for e in range (setElementMax):
                    nodeList.append (points[setElement[e, 0]])
                    nodeList.append (points[setElement[e, 1]])
                    if flagWriteFrames == True:
                        frameList.append (np.concatenate ((setFrame[e]["s0"], setFrame[e]["s1"], setFrame[e]["s2"])))       # note need to repeat frame data for each node of element                 
                        frameList.append (np.concatenate ((setFrame[e]["s0"], setFrame[e]["s1"], setFrame[e]["s2"])))                        
                    elementList.append ([n, n+1])
                    n = n + 2                     
            elif setType == "triangle":
                typeLabel = "FETRIANGLE"
                for e in range (setElementMax):
                    nodeList.append (points[setElement[e, 0]])
                    nodeList.append (points[setElement[e, 1]])
                    nodeList.append (points[setElement[e, 2]])
                    if flagWriteFrames == True:
                        frameList.append (np.concatenate ((setFrame[e]["s0"], setFrame[e]["s1"], setFrame[e]["s2"])))       # note need to repeat frame data for each node of element                      
                        frameList.append (np.concatenate ((setFrame[e]["s0"], setFrame[e]["s1"], setFrame[e]["s2"])))   
                        frameList.append (np.concatenate ((setFrame[e]["s0"], setFrame[e]["s1"], setFrame[e]["s2"])))
                    elementList.append ([n, n+1, n+2])
                    n = n + 3
            nodeMax = len (nodeList)
            elementMax = len (elementList)
            
            # Convert lists of np arrays to full np arrays to simplify indexing
            nodeList = np.array (nodeList)
            elementList = np.array (elementList)
            if flagWriteFrames == True:
                frameList = np.array (frameList)

            # Write zone header
            f.write ("\n")
            f.write ("ZONE T=\"{}\" N={}, E={}, DATAPACKING=BLOCK, ZONETYPE={}\n".format (setName, nodeMax, elementMax, typeLabel))
            
            # Write node var data
            for xyz in range (3):
                count = 0
                for n in range (nodeMax):
                    f.write ("{} ".format (nodeList[n, xyz]))
                    if count == 9:                             # Limit number of nodes per line to avoid exceeding 4000 char limit in tecplot
                        count = 0
                        f.write ("\n")   
                    else:
                        count = count + 1
                f.write ("\n")
                
            # Write frame var data
            if flagWriteFrames == True:
                for fxyz in range (9):
                    count = 0
                    for n in range (nodeMax):
                        f.write ("{} ".format (frameList[n, fxyz]))
                        if count == 9:                             # Limit number of nodes per line to avoid exceeding 4000 char limit in tecplot
                            count = 0
                            f.write ("\n")   
                        else:
                            count = count + 1
                    f.write ("\n")
            
            # Write element index data
            for element in elementList:
                for node in element:
                    f.write ("{} ".format (node+1))             # make node indices 1-offset for Tecplot
                f.write ("\n")
            
            
        # Close output file
        f.close()

        
    
    
    def createMeshGenData (self, meshType):
        """ Create meshGen data structure to store node and element data for output
            Inputs:
                meshType                    structure or fluid
            Saves data to class members of data[meshType]:
                "nodes"                     2D numpy array of all node coordinates of mesh (rows=nodes, cols=x,y,z)
                cell_set                    set of dictionaries (one per named mesh entity), each containing:
                    "type"                  line or triangle
                    "elements"              2D numpy array of node indicies for each element of current entity (rows=elements, cols=node ids for current element)
        """
        import numpy as np
        
        # Init data
        self.data[meshType] = {}
        cell_set = {}
        
     
        # NODES
        
        # Get node list for entire mesh
        #   note gmsh returns lists nodeTags=[1, 2, 3, ...], coord=[x1, y1, z1, x2, y2, z2, ...], and parametricCoord (not used)
        nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
        node = coord.reshape ((int(len(coord)/3), 3))                                                   # reshape list to be 2D array with columns for xyz
        
        
        # ELEMENTS
        
        # Get list of named entities for 1D (line) and 2D (surface) elements
        entityList = gmsh.model.getPhysicalGroups (dimLine) + gmsh.model.getPhysicalGroups (dimSurf)
        
        # Iterate over each entity to save as cell set data
        for entity in entityList:
        
            # Get dim, tag, name for current cell set
            dim, tagPhysical = entity[0], entity[1]
            name = gmsh.model.getPhysicalName (dim, tagPhysical)                                        # get name of set

            # Only output traingles for fluid meshes
            if meshType == "structure"  or  (meshType == "fluid"  and  dim == dimSurf):
                
                # Create cell set block for current item
                cell_set[name] = {}
                
                # Get list of entity tags for current cell set (note in gmsh tagPhysical != tagEntity)
                tagList = gmsh.model.getEntitiesForPhysicalGroup (dim, tagPhysical)                         # get line or surface tags

                # Iterate over each entity in physical group and gather element node data
                for i in range (len (tagList)):
                    elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements (dim, tagList[i])     # get elements and node id list, note gmsh stores nodeTags=[e1n1, e1n2, e1n3, e2n1, ...] for surface or only 2 nodes per eelment for lines
                    nodeTags = np.hstack (nodeTags) - 1                                                     # horz stack array to remove array inside array and make zero offset for python
                    if i == 0:
                        elementList = nodeTags                                                              # if first tag, init elementList
                    else:
                        elementList = np.append (elementList, nodeTags)                                     # append current element node ids to elementList
                
                # Convert element list to signed int dtype to avoid issues down the road with arithmetic
                elementList = elementList.astype(int)
                
                # Use elementType to determine number of nodes per element
                if elementTypes.item() == 1:
                    cell_set[name]["type"] = "line"
                    numNodes = 2
                elif elementTypes.item() == 2:
                    cell_set[name]["type"] = "triangle"
                    numNodes = 3
                elif elementTypes.item() == 3:
                    cell_set[name]["type"] = "quad"
                    numNodes = 4
                else:
                    raise Exception ("Unknown elementType! Exiting")
                
                # Copy element data to cell set
                cell_set[name]["elements"] = elementList.reshape ((int(len(elementList)/numNodes), numNodes))   # reshape list to be 2D array with columns for node ids        
            

        # Save data to class
        self.data[meshType]["node"] = node
        self.data[meshType]["cell_set"] = cell_set