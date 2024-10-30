import pytest
from src.structures import WindFarm, Domain
from src.functions2D import *
from src.terrain import buildTerrain2D

@pytest.fixture
def domain():
    d = Domain()
    d.setDomain([-1200, 1200], [-1200, 1200], 0)
    return d

def test_turbinegen_full(domain):
    gmsh.model.add("2D Turbine Generation Test")

    params = dict()
    params['domain'] = {
        'x_range': [-1200, 1200],
        'y_range': [-1200, 1200]
    }
    params['refine'] = {
        'background_length_scale': 100,
        'turbine': {},
        'farm': {}
    }

    params['refine']['turbine']['num_turbines'] = 2
    params['refine']['turbine']['length_scale'] = 30
    params['refine']['background_length_scale'] = 100
    params['refine']['farm']['length_scale'] = 70
    params['refine']['turbine']['threshold_upstream_distance'] = 240
    params['refine']['turbine']['threshold_downstream_distance'] = 300
    params['refine']['turbine']['threshold_rotor_distance'] = 100
    params['refine']['turbine']['shudder'] = 50

    params['refine']['turbine'][1] = {
        'x': 200,
        'y': 400,
        'wake': 1
    }
    params['refine']['turbine'][2] = {
        'x': -800,
        'y': 300,
        'wake': 0
    }

    struct = buildTerrain2D(params, domain)
    gmsh.model.geo.addPlaneSurface([struct])
    wf = WindFarm()

    turbines = buildFarms2D(params, wf, domain)

    gmsh.model.remove()

    assert len(turbines) == 2
    assert wf.x_range == [-1070, 250]
    assert wf.y_range == [130, 670]

    
