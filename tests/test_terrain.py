import gmsh
from aeromesh.terrain.terrain import buildTerrain2D, buildTerrainDefault, buildTerrainFromFile
from aeromesh.structs.structures import Domain

def test_terrain3D():
    gmsh.model.add("3D Terrain Test")
    params = dict()

    params['domain'] = {
        'terrain_path': './tests/infiles/skew_terrain.txt',
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
        'z_range': [0, 1000],
        'aspect_ratio': 1,
        'upper_aspect_ratio': 1,
        'aspect_distance': 0
    }

    params['refine'] = {
        'background_length_scale': 100
    }

    buildTerrainFromFile(params, Domain())
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    entities2D = gmsh.model.getEntities(dim=2)
    entities3D = gmsh.model.getEntities(dim=3)

    gmsh.model.remove()
    
    assert len(entities2D) == 6
    assert len(entities3D) == 0

def test_terrain3D_default():
    gmsh.model.add("3D Default Terrain Test")
    params = dict()

    params['domain'] = {
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
        'z_range': [0, 1000],
        'aspect_ratio': 1,
        'upper_aspect_ratio': 1,
        'aspect_distance': 0
    }

    params['refine'] = {
        'background_length_scale': 100
    }

    buildTerrainDefault(params, Domain())
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    entities2D = gmsh.model.getEntities(dim=2)
    entities3D = gmsh.model.getEntities(dim=3)
    gmsh.model.remove()

    assert len(entities2D) == 6
    assert len(entities3D) == 0

def test_terrain2D():
    gmsh.model.add("2D Terrain Test")
    params = dict()

    d = Domain()
    d.setDomain([-1200, 1200], [-1200, 1200])

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

    entities1D = gmsh.model.getEntities(dim=1)
    entities2D = gmsh.model.getEntities(dim=2)
    gmsh.model.remove()

    assert len(entities2D) == 1
    assert len(entities1D) == 4