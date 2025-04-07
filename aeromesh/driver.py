import gmsh
import sys, os
import yaml
import aeromesh.structs as structs
import aeromesh.terrain as tr
import aeromesh.geometry as geometry
import meshio
import numpy as np
import math
import tempfile

def toXDMF(ndim, name):
    tempdir = tempfile.TemporaryDirectory()
    gmsh.write(os.path.join(tempdir.name, "dummy.msh"))

    msh = meshio.read(os.path.join(tempdir.name, "dummy.msh"))

    def create_mesh(mesh, type, prune_dim=False):
        cells = mesh.get_cells_type(type)
        cell_data = mesh.get_cell_data("gmsh:physical", type)
        points = mesh.points if prune_dim == False else mesh.points[:, :2]
        out_mesh = meshio.Mesh(points=points, cells={type: cells}, cell_data={"facet_tags": [cell_data.astype(np.int32)]})
        return out_mesh

    if ndim == 3:
        tetras = create_mesh(msh, 'tetra')
        tris = create_mesh(msh, 'triangle')
        meshio.write(name + '.xdmf', tetras)
        meshio.write(name + '_boundary.xdmf', tris)

    else:
        tris = create_mesh(msh, 'triangle', prune_dim=True)
        lines = create_mesh(msh, 'line', prune_dim=True)
        meshio.write(name + '.xdmf', tris)
        meshio.write(name + '_boundary.xdmf', lines)

    tempdir.cleanup()

def generate2DMesh(params):

    """
    Executes a 2D meshing operation on params. Outputs results in the file out.vtk.

    :param params: A dictionary of the required parameters.
    :type params: dict()

    """
    gmsh.model.add("User Model")

    scale = params['refine']['global_scale']
    params['refine']['turbine']['length_scale'] *= scale
    params['refine']['farm']['length_scale'] *= scale
    params['refine']['background_length_scale'] *= scale
    farmType = params['domain']['type']

    domain = structs.Domain()

    wf = structs.WindFarm()
    if farmType == 'box':
        farmBorder = tr.buildTerrain2D(params, domain)
    elif farmType == 'cylinder':
        farmBorder = tr.buildTerrainCircle(params, domain)
    else:
        raise Exception("Invalid farm type specified. Farm types must be in [box, cylinder].")
    
    farm = gmsh.model.geo.addPlaneSurface([farmBorder], tag=999)
    sp1 = gmsh.model.geo.addPhysicalGroup(2, [farm], tag=0)

    fields = [999]
    fields.extend(geometry.buildFarms2D(params, wf, domain))
    fields.extend(geometry.generateCustomRefines(params))
    if params['refine']['farm']['type'] != 'none':
        fields.append(998)
        geometry.refineFarm2D(params, wf)

    mesher = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mesher, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mesher)

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.removeDuplicateElements()

    gmsh.model.mesh.generate(2)

def generate3DMesh(params):

    """
    Executes a 3D meshing operation on params. Outputs results in the file out.vtk.

    :param params: A dictionary of the required parameters.
    :type params: dict()

    """

    gmsh.model.add("User Model")

    scale = params['refine']['global_scale']
    farmType = params['domain']['type']
    params['refine']['turbine']['length_scale'] *= scale
    params['refine']['farm']['length_scale'] *= scale
    params['refine']['background_length_scale'] *= scale

    domain = structs.Domain()
    terrain = None
    if farmType == 'box':
        try:
            terrainDefined = params['domain']['terrain_path']
            terrain = tr.buildTerrainFromFile(params, domain)
        except KeyError:
            terrain = tr.buildTerrainDefault(params, domain)
    elif farmType == 'cylinder':
        terrain = tr.buildTerrainCylinder(params, domain)
    else:
        raise Exception("Invalid farm type specified. Farm types must be in [box, cylinder].")

    if farmType == 'box':
        v1 = gmsh.model.geo.addVolume([terrain], tag=1)
        vp1 = gmsh.model.geo.addPhysicalGroup(3, [v1], tag=0)
        
    wf = structs.WindFarm()
    fields = geometry.generateTurbines(params, domain, wf)
    fields.append(999) #Background field, reserved number
    farmRefine = geometry.refineFarm3D(params, wf)
    if farmRefine:
        fields.append(farmRefine)
    fields.extend(geometry.generateCustomRefines(params))

    gmsh.model.geo.synchronize()
    mesher = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mesher, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mesher)

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.removeDuplicateElements()

    gmsh.model.mesh.generate(3)

    geometry.anisotropyScale(params)

    if farmType == 'cylinder' and domain.interp:
        geometry.cylinderTerrainAdjustment(domain, params)

def main():
    gmsh.initialize()

    if len(sys.argv) < 2:
        raise Exception("Input file not specified.")
    filename = sys.argv[1]
    with open(filename) as input:
        params = yaml.safe_load(input)

    setYAMLDefaults(params)
    verifyYAML(params)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 1)
    if params['suppress_out'] > 0:
        gmsh.option.setNumber("General.Verbosity", 0)

    params['domain']['inflow_angle'] *= math.pi / 180

    if params['domain']['dimension'] == 3:
        generate3DMesh(params)
    else:
        generate2DMesh(params)
    gmsh.model.mesh.optimize()

    filename = params['filename']
    
    if params['filetype'] != 'xdmf':
        filename = filename + params['filetype']
        gmsh.write(filename)
    else:
        ndim = params['domain']['dimension']
        toXDMF(ndim, filename)

    if '-v' in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

def runAeroMesh(params):
    gmsh.initialize()

    setYAMLDefaults(params)
    verifyYAML(params)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 1)
    if params['suppress_out'] > 0:
        gmsh.option.setNumber("General.Verbosity", 0)

    params['domain']['inflow_angle'] *= math.pi / 180

    if params['domain']['dimension'] == 3:
        generate3DMesh(params)
    else:
        generate2DMesh(params)
    gmsh.model.mesh.optimize()

    filename = params['filename']
    
    if params['filetype'] != 'xdmf':
        filename = filename + params['filetype']
        gmsh.write(filename)
    else:
        ndim = params['domain']['dimension']
        toXDMF(ndim, filename)

    gmsh.finalize()

def setYAMLDefaults(params):
    refine = params['refine']
    domain = params['domain']

    params.setdefault('filetype', 'msh')
    params.setdefault('filename', 'out')
    params.setdefault('suppress_out', 1)
    params.setdefault('refine_custom', {}).setdefault('num_refines', 0)

    domain.setdefault('aspect_ratio', 1)
    domain.setdefault('upper_aspect_ratio', 1)
    domain.setdefault('aspect_distance', 1)
    domain.setdefault('inflow_angle', 0)
    domain.setdefault('type', 'box')

    refine.setdefault('global_scale', 1)
    refine.setdefault('turbine', {}).setdefault('num_turbines', 0)
    refine.setdefault('turbine', {}).setdefault('type', 'simple')

    refine.setdefault('farm', {}).setdefault('length_scale', params['refine']['background_length_scale'])
    refine.setdefault('farm', {}).setdefault('threshold_distance', 0)
    refine.setdefault('farm', {}).setdefault('type', 'none')

def verifyYAML(params):
    err = 0
    suppress = params['suppress_out']
    if not suppress:
        print("***----------------------------------------***")
        print("Validating YAML file.")
    for key in params:
        if key not in ['refine', 'domain', 'filetype', 'filename', 'refine_custom', 'suppress_out']:
            print("AeroMesh: Unknown field: " + key)
            err = 1
    domainChecks = params['domain']
    refineChecks = params['refine']
    customChecks = params['refine_custom']
    for key in domainChecks:
        valid = ['terrain_path', 'x_range', 'y_range', 'z_range', 'aspect_ratio', 'upper_aspect_ratio',
                 'aspect_distance', 'dimension', 'inflow_angle', 'type', 'center', 'radius']
        if key not in valid:
            print("AeroMesh: Unknown field: " + key)
            err = 1
    for key in refineChecks:
        valid = ['turbine', 'background_length_scale', 'farm', 'global_scale']
        if key not in valid:
            print("AeroMesh: Unknown field: " + key)
            err = 1
    turbineChecks = params['refine']['turbine']
    for key in turbineChecks:
        validNums = [i for i in range(1, turbineChecks['num_turbines'] + 1)]
        validParams = ['num_turbines', 'length_scale', 'threshold_upstream_distance', 'threshold_downstream_distance',
                       'threshold_rotor_distance', 'type']
        if key in validNums:
            validSubkeys = ['x', 'y', 'HH']
            for subkey in turbineChecks[key]:
                if subkey not in validSubkeys:
                    print("AeroMesh: Unknown field: " + str(key))
                    err = 1
        elif key not in validParams:
            print("AeroMesh: Unknown field: " + str(key))
            err = 1
    for key in customChecks:
        validNums = [i for i in range(1, customChecks['num_refines'] + 1)]
        if key == 'num_refines':
            continue
        elif key in validNums:
            validSubkeys = ['type', 'x_range', 'y_range', 'radius', 'length_scale', 'z_range', 'length']
            for subkey in customChecks[key]:
                if subkey not in validSubkeys:
                    print("AeroMesh: Unknown  refine_custom[" + str(key) + "] field: " + str(subkey))
                    err = 1
        else:
            print("AeroMesh: Unknown  refine_custom field: " + str(key))
            err = 1
    farmChecks = params['refine']['farm']
    for key in farmChecks:
        valid = ['length_scale', 'threshold_distance', 'type']
        if key not in valid:
            print("AeroMesh: Unknown field: " + key)
            err = 1
    if err == 0:
        if not suppress:
            print("YAML validated successfully.")
            print("***----------------------------------------***")
    else:
        print("***----------------------------------------***")
        raise Exception("YAML Error: Unexpected fields specified.")

if __name__ == "__main__":
    main()