import gmsh
import sys
import yaml
from src.structures import Domain, WindFarm
from src.terrain import buildTerrainFromFile, buildTerrainDefault, buildTerrain2D, buildTerrainCylinder
from src.functions2D import *
from src.functions3D import *
from src.refines import generateCustomRefines
import meshio
import numpy as np
import math

def toXDMF(ndim):
    elements = gmsh.model.mesh.getElements()
    loc = 2 if ndim == 3 else 1
    locTets = 4 if ndim == 3 else 2
    indexTris = np.where(elements[0] == loc)[0][0]
    indexTets = np.where(elements[0] == locTets)[0][0]
    points = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
    cells = elements[2][indexTets].reshape(-1, ndim + 1) - 1
    triangles = elements[2][indexTris].reshape(-1, ndim) - 1

    faces = gmsh.model.getEntities(dim=ndim - 1)
    tags = []
    for face in faces:
        target = face[1]
        tag = 0
        elements = gmsh.model.mesh.getElements(dim=ndim - 1, tag=target)
        match target:
            case 990:
                tag = 6
            case 989:
                tag = 5
            case 992:
                tag = 3
            case 993:
                tag = 2
            case 994:
                tag = 1
            case 995:
                tag = 4
            case _:
                tag = 0
        for _ in elements[1][0]:
            tags.append(tag)

    tags = np.array(tags).astype(int)

    output_mesh_name = "out.xdmf"
    output_boundary_mesh_name = output_mesh_name.split(".")[0] + "_boundary.xdmf"
    if ndim == 3:
        tetra_mesh = meshio.Mesh(points=points, cells={"tetra": cells})
        boundary_mesh = meshio.Mesh(points=points, cells=[("triangle", triangles)], cell_data={"facet_tags":[tags]})
    else:
        tetra_mesh = meshio.Mesh(points=points, cells={"triangle": cells})
        boundary_mesh = meshio.Mesh(points=points, cells=[("line", triangles)], cell_data={"facet_tags":[tags]})

    meshio.write(output_mesh_name, tetra_mesh)
    meshio.write(output_boundary_mesh_name, boundary_mesh)

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

    domain = Domain()

    wf = WindFarm()

    farmBorder = buildTerrain2D(params, domain)
    farm = gmsh.model.geo.addPlaneSurface([farmBorder], tag=999)

    fields = [999]
    fields.extend(buildFarms2D(params, wf, domain))
    fields.extend(generateCustomRefines(params))
    if params['refine']['farm']['type'] != 'none':
        fields.append(998)
        refineFarm2D(params, wf)

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

    domain = Domain()
    terrain = None
    if farmType == 'box':
        try:
            terrainDefined = params['domain']['terrain_path']
            terrain = buildTerrainFromFile(params, domain)
        except KeyError:
            terrain = buildTerrainDefault(params, domain)
    elif farmType == 'cylinder':
        terrain = buildTerrainCylinder(params, domain)
    else:
        raise Exception("Invalid farm type specified. Farm types must be in [box, cylinder].")

    if farmType == 'box':
        v1 = gmsh.model.geo.addVolume([terrain], tag=1)
        
    wf = WindFarm()
    fields = generateTurbines(params, domain, wf)
    fields.append(999) #Background field, reserved number
    farmRefine = refineFarm3D(params, wf)
    if farmRefine:
        fields.append(farmRefine)
    fields.extend(generateCustomRefines(params))

    gmsh.model.geo.synchronize()
    mesher = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mesher, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mesher)

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.removeDuplicateElements()

    gmsh.model.mesh.generate(3)

    anisotropyScale(params)

    if farmType == 'cylinder' and domain.interp:
        cylinderTerrainAdjustment(domain, params)

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

    params['domain']['inflow_angle'] *= math.pi / 180

    if params['domain']['dimension'] == 3:
        generate3DMesh(params)
    else:
        generate2DMesh(params)
    
    if params['filetype'] != 'xdmf':
        filename = 'out.' + params['filetype']
        gmsh.write(filename)
    else:
        ndim = params['domain']['dimension']
        toXDMF(ndim)

    if '-v' in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

def setYAMLDefaults(params):
    refine = params['refine']
    domain = params['domain']

    params.setdefault('filetype', 'msh')
    params.setdefault('refine_custom', {}).setdefault('num_refines', 0)

    domain.setdefault('aspect_ratio', 1)
    domain.setdefault('upper_aspect_ratio', 1)
    domain.setdefault('aspect_distance', 0)
    domain.setdefault('inflow_angle', 0)
    domain.setdefault('type', 'box')

    refine.setdefault('global_scale', 1)
    refine.setdefault('turbine', {}).setdefault('num_turbines', 0)
    refine.setdefault('turbine', {}).setdefault('type', 'rectangle')

    refine.setdefault('farm', {}).setdefault('length_scale', params['refine']['background_length_scale'])
    refine.setdefault('farm', {}).setdefault('threshold_distance', 0)
    refine.setdefault('farm', {}).setdefault('type', 'none')

def verifyYAML(params):
    err = 0
    print("***----------------------------------------***")
    print("Validating YAML file.")
    for key in params:
        if key not in ['refine', 'domain', 'filetype', 'refine_custom']:
            print("Unknown field: " + key)
            err = 1
    domainChecks = params['domain']
    refineChecks = params['refine']
    customChecks = params['refine_custom']
    for key in domainChecks:
        valid = ['terrain_path', 'x_range', 'y_range', 'height', 'aspect_ratio', 'upper_aspect_ratio',
                 'aspect_distance', 'dimension', 'inflow_angle', 'type', 'center', 'radius']
        if key not in valid:
            print("Unknown field: " + key)
            err = 1
    for key in refineChecks:
        valid = ['turbine', 'background_length_scale', 'farm', 'global_scale']
        if key not in valid:
            print("Unknown field: " + key)
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
                    print("Unknown field: " + str(key))
                    err = 1
        elif key not in validParams:
            print("Unknown field: " + str(key))
            err = 1
    for key in customChecks:
        validNums = [i for i in range(1, customChecks['num_refines'] + 1)]
        if key == 'num_refines':
            continue
        elif key in validNums:
            validSubkeys = ['type', 'x_range', 'y_range', 'radius', 'length_scale', 'z_range']
            for subkey in customChecks[key]:
                if subkey not in validSubkeys:
                    print("Unknown refine_custom[" + str(key) + "] field: " + str(subkey))
                    err = 1
        else:
            print("Unknown refine_custom field: " + str(key))
            err = 1
    farmChecks = params['refine']['farm']
    for key in farmChecks:
        valid = ['length_scale', 'threshold_distance', 'type']
        if key not in valid:
            print("Unknown field: " + key)
            err = 1
    if err == 0:
        print("YAML validated successfully.")
        print("***----------------------------------------***")
    else:
        print("***----------------------------------------***")
        raise Exception("YAML Error: Unexpected fields specified.")

if __name__ == "__main__":
    main()