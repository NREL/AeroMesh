import gmsh
from src.functions3D import *
from src.structures import WindFarm, Domain
import pytest
from src.terrain import buildTerrainDefault

###
# Fixtures
###
@pytest.fixture
def domain():
    d = Domain()
    d.setDomain([-1200, 1200], [-1200, 1200], 1000)
    return d

@pytest.fixture
def params():
    params = dict()
    params['domain'] = {
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200],
        'height': 1000,
        'aspect_ratio': 1,
        'upper_aspect_ratio': 1,
        'aspect_distance': 0,
        'inflow_angle': 0
    }
    params['refine'] = {
        'background_length_scale': 100,
        'turbine': {},
        'farm': {}
    }
    return params


###
# Unit Tests
###
def test_turbinegen_isolated(domain, params):
    gmsh.model.add("3D Turbine Generation Test")
    
    sl = buildTerrainDefault(params, domain)
    gmsh.model.geo.addVolume([sl], tag=1)

    wf = WindFarm()

    isotropic_test = placeTurbine(500, 500, 100, 240, 300, 
                                  100, 30, 100, 50, 0, 1, wf, domain)
    
    anisotropic_test = placeTurbine(-700, -700, 100, 240, 300, 
                                  100, 30, 100, 50, 0, 3, wf, domain)
    
    gmsh.model.remove()

    assert len(isotropic_test) == 1
    
    assert len(anisotropic_test) == 4

def test_turbinegen_full(domain, params):
    gmsh.model.add("3D Turbine Generation Test (Full)")
    wf = WindFarm()

    sl = buildTerrainDefault(params, domain)
    gmsh.model.geo.addVolume([sl], tag=1)

    params['refine']['turbine']['num_turbines'] = 2
    params['refine']['turbine']['length_scale'] = 30
    params['refine']['background_length_scale'] = 100
    params['refine']['farm']['length_scale'] = 70
    params['refine']['turbine']['threshold_upstream_distance'] = 240
    params['refine']['turbine']['threshold_downstream_distance'] = 300
    params['refine']['turbine']['threshold_rotor_distance'] = 100
    params['domain']['aspect_ratio'] = 2

    params['refine']['turbine'][1] = {
        'x': 500,
        'y': 800,
        'wake': 1
    }
    params['refine']['turbine'][2] = {
        'x': -1000,
        'y': 300,
        'wake': 0
    }

    fields_test = generateTurbines(params, domain, wf)
    gmsh.model.remove()

    assert len(fields_test) == 6
    assert wf.x_range == [-1200, 800]
    assert wf.y_range == [300, 800]
    assert wf.zMax == 400

def test_anisotropy(domain, params):
    gmsh.model.add("3D Anisotropy Test")

    params['domain']['aspect_ratio'] = 2
    params['domain']['aspect_distance'] = 400

    sl = buildTerrainDefault(params, domain)
    gmsh.model.geo.addVolume([sl], tag=1)

    anisotropyScale(params)

    data = gmsh.model.mesh.getNodes()
    tags = data[0]
    coords = data[1].reshape(-1, 3)

    gmsh.model.remove()

    for node in zip(tags, coords):
        coord = node[1]
        assert coord[2] <= params['domain']['height']
