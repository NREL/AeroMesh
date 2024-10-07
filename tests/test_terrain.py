import gmsh
from src.terrain import buildTerrain2D, buildTerrainDefault, buildTerrainFromFile
from src.structures import Domain
from filecmp import cmp

def test_terrain3D():
    gmsh.model.add("3D Terrain Test")
    params = dict()

    params['domain'] = {
        'terrain_path': './tests/infiles/skew_terrain.txt',
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
        'height': 1000,
        'aspect_ratio': 1,
        'aspect_distance': 0
    }

    params['refine'] = {
        'background_length_scale': 100
    }

    buildTerrainFromFile(params, Domain())
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    filename = './tests/outfiles/out_terrain3D.vtk'
    gmsh.write(filename)
    gmsh.model.remove()

    assert cmp('./tests/outfiles/out_terrain3D.vtk', './tests/testcases/case_terrain3D.vtk') is True

def test_terrain3D_default():
    gmsh.model.add("3D Default Terrain Test")
    params = dict()

    params['domain'] = {
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
        'height': 1000,
        'aspect_ratio': 1,
        'aspect_distance': 0
    }

    params['refine'] = {
        'background_length_scale': 100
    }

    buildTerrainDefault(params, Domain())
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    filename = './tests/outfiles/out_terrain3D_default.vtk'
    gmsh.write(filename)
    gmsh.model.remove()

    assert cmp('./tests/outfiles/out_terrain3D_default.vtk', './tests/testcases/case_terrain3D_default.vtk') is True

def test_terrain2D():
    gmsh.model.add("2D Terrain Test")
    params = dict()

    d = Domain()
    d.setDomain([-1200, 1200], [-1200, 1200], 0)

    params['domain'] = {
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
    }

    params['refine'] = {
        'background_length_scale': 100
    }

    struct = buildTerrain2D(params, d)
    gmsh.model.geo.addPlaneSurface([struct])
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    filename = './tests/outfiles/out_terrain2D.vtk'
    gmsh.write(filename)
    gmsh.model.remove()

    assert cmp('./tests/outfiles/out_terrain2D.vtk', './tests/testcases/case_terrain2D.vtk') is True