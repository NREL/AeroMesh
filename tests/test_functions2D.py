import pytest
from aeromesh.structs.structures import WindFarm, Domain
from aeromesh.geometry.functions2D import *
from aeromesh.terrain.terrain import buildTerrain2D

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
        'y_range': [-1200, 1200],
        'inflow_angle': 0
    }
    params['refine'] = {
        'background_length_scale': 100,
        'turbine': {
            'type': 'circle'
        },
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
    }
    params['refine']['turbine'][2] = {
        'x': -800,
        'y': 300,
    }

    struct = buildTerrain2D(params, domain)
    gmsh.model.geo.addPlaneSurface([struct])
    wf = WindFarm()

    turbines = buildFarms2D(params, wf, domain)

    gmsh.model.remove()

    assert len(turbines) == 2
    assert wf.x_range == [-900, 300]
    assert wf.y_range == [200, 500]

    
