import gmsh
import sys
import yaml
from src.structures import Domain, WindFarm
from src.terrain import buildTerrainFromFile, buildTerrainDefault, buildTerrain2D
from src.functions2D import *
from src.functions3D import *


def generate2DMesh(params):

    """
    Executes a 2D meshing operation on params. Outputs results in the file out.vtk.

    :param params: A dictionary of the required parameters.
    :type params: dict()

    """
        
    gmsh.initialize()

    gmsh.model.add("User Model")

    scale = params['refine']['global_scale']
    params['refine']['turbine']['length_scale'] *= scale
    params['refine']['farm']['length_scale'] *= scale
    params['refine']['background_length_scale'] *= scale

    farm = []

    domain = Domain()

    wf = WindFarm()

    farm.append(buildTerrain2D(params, domain))
    farm.extend(buildFarms2D(params, wf, domain))

    gmsh.model.geo.addPlaneSurface(farm)
    refineFarm2D(params, wf)
    gmsh.option.setNumber("Mesh.Algorithm", 8)

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(2)

def generate3DMesh(params):

    """
    Executes a 3D meshing operation on params. Outputs results in the file out.vtk.

    :param params: A dictionary of the required parameters.
    :type params: dict()

    """

    gmsh.initialize()

    gmsh.model.add("User Model")

    scale = params['refine']['global_scale']
    params['refine']['turbine']['length_scale'] *= scale
    params['refine']['farm']['length_scale'] *= scale
    params['refine']['background_length_scale'] *= scale

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    domain = Domain()
    try:
        terrainDefined = params['domain']['terrain_path']
        terrain = buildTerrainFromFile(params, domain)
    except KeyError:
        terrain = buildTerrainDefault(params, domain)

    v1 = gmsh.model.geo.addVolume([terrain], tag=999)
    wf = WindFarm()
    fields = generateTurbines(params, domain, wf)
    fields.append(999) #Background field, reserved number
    fields.append(refineFarm3D(params, wf))

    gmsh.model.geo.synchronize()
    mesher = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mesher, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mesher)

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.removeDuplicateElements()

    gmsh.model.mesh.generate(3)

    anisotropyScale(params)

def main():
    if len(sys.argv) < 2:
        raise Exception("Input file not specified.")
    filename = sys.argv[1]
    with open(filename) as input:
        params = yaml.safe_load(input)

    setYAMLDefaults(params)
    verifyYAML(params)

    if params['domain']['dimension'] == 3:
        generate3DMesh(params)
    else:
        generate2DMesh(params)
    
    filename = 'out.msh'
    gmsh.write(filename)

    gmsh.finalize()

def setYAMLDefaults(params):
    refine = params['refine']
    domain = params['domain']


    refine.setdefault('global_scale', 1)
    domain.setdefault('aspect_ratio', 1)

    domain.setdefault('aspect_distance', 0)
    refine.setdefault('turbine', {}).setdefault('num_turbines', 0)
    refine.setdefault('turbine', {}).setdefault('shudder', params['refine']['turbine']['threshold_rotor_distance'])

    refine.setdefault('farm', {}).setdefault('length_scale', params['refine']['background_length_scale'])
    refine.setdefault('farm', {}).setdefault('threshold_distance', 0)

def verifyYAML(params):
    err = 0
    print("***----------------------------------------***")
    print("Validating YAML file.")
    for key in params:
        if key not in ['refine', 'domain']:
            print("Unknown field: " + key)
            err = 1
    domainChecks = params['domain']
    refineChecks = params['refine']
    for key in domainChecks:
        valid = ['terrain_path', 'x_range', 'y_range', 'height', 'aspect_ratio', 'aspect_distance', 'dimension']
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
                       'threshold_rotor_distance', 'shudder']
        if key in validNums:
            validSubkeys = ['x', 'y', 'wake']
            for subkey in turbineChecks[key]:
                if subkey not in validSubkeys:
                    print("Unknown field: " + str(key))
                    err = 1
        elif key not in validParams:
            print("Unknown field: " + str(key))
            err = 1
    farmChecks = params['refine']['farm']
    for key in farmChecks:
        valid = ['length_scale', 'threshold_distance']
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