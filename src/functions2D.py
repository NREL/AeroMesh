import gmsh

####
## 2D Meshing Functions
####

def generateTurbine2D(x, y, lcTurbine, dRotor, dWake, shudder, wake, wf):
    curve = []
    
    if wake == 0:
        m1 = gmsh.model.geo.addPoint(x + dWake / 2, y + shudder / 2, 0, lcTurbine)
        m1p = gmsh.model.geo.addPoint(x + dWake / 2, y - shudder / 2, 0, lcTurbine)
        m2 = gmsh.model.geo.addPoint(x, y + dRotor / 2, 0, lcTurbine)
        m3 = gmsh.model.geo.addPoint(x - dWake / 2, y + shudder / 2, 0, lcTurbine)
        m3p = gmsh.model.geo.addPoint(x - dWake / 2, y - shudder / 2, 0, lcTurbine)
        m4 = gmsh.model.geo.addPoint(x, y - dRotor / 2, 0, lcTurbine)

        wf.updateXMax(x + dWake / 2)
        wf.updateXMin(x - dWake / 2)
        wf.updateYMax(y + dRotor / 2)
        wf.updateYMin(y - dRotor / 2)

    else:
        m1 = gmsh.model.geo.addPoint(x + shudder / 2, y + dWake / 2, 0, lcTurbine)
        m1p = gmsh.model.geo.addPoint(x - shudder / 2, y + dWake / 2, 0, lcTurbine)
        m2 = gmsh.model.geo.addPoint(x + dRotor / 2, y, 0, lcTurbine)
        m3 = gmsh.model.geo.addPoint(x + shudder / 2, y - dWake / 2, 0, lcTurbine)
        m3p = gmsh.model.geo.addPoint(x - shudder / 2, y - dWake / 2, 0, lcTurbine)
        m4 = gmsh.model.geo.addPoint(x - dRotor / 2, y, 0, lcTurbine)

        wf.updateXMax(x + dRotor / 2)
        wf.updateXMin(x - dRotor / 2)
        wf.updateYMax(y + dWake / 2)
        wf.updateYMin(y - dWake / 2)

    turbine = gmsh.model.geo.addPoint(x, y, 0, lcTurbine)

    curve.append(gmsh.model.geo.addLine(m1p, m1))
    curve.append(gmsh.model.geo.addLine(m1, m2))
    curve.append(gmsh.model.geo.addLine(m2, m3))
    curve.append(gmsh.model.geo.addLine(m3, m3p))
    curve.append(gmsh.model.geo.addLine(m3p, m4))
    curve.append(gmsh.model.geo.addLine(m4, m1p))

    return gmsh.model.geo.addCurveLoop(curve, reorient=True)

def buildFarms2D(params, wf, domain):
    turbines = []

    nFarms = params['refine']['turbine']['num_turbines']
    lcTurbine = params['refine']['turbine']['length_scale'] 
    upstream = params['refine']['turbine']['threshold_upstream_distance']
    downstream = params['refine']['turbine']['threshold_downstream_distance']
    rotor = params['refine']['turbine']['threshold_rotor_distance']
    shudder = params['refine']['turbine']['shudder']

    for i in range(nFarms):
        turbineData = params['refine']['turbine'][i + 1]
        x = turbineData['x'] 
        y = turbineData['y']
        wake = turbineData['wake'] 

        if not domain.withinDomain(x, y):
            raise Exception("Invalid turbine location: (" + str(x) + ", " + str(y)) + ")."
        
        turbine = generateTurbine2D(x, y, lcTurbine, rotor, upstream + downstream, shudder, wake, wf)
        turbines.append(turbine)
        gmsh.model.geo.addPlaneSurface([turbine])

    return turbines

def refineFarm2D(params, wf):
    dist = params['refine']['farm']['threshold_distance'] 
    wf.adjustDistance(dist)
    farmLC = params['refine']['farm']['length_scale'] 
    blc = params['refine']['background_length_scale']

    b = gmsh.model.mesh.field.add("Box", tag=998)
    gmsh.model.mesh.field.setNumber(b, "XMin", wf.x_range[0])
    gmsh.model.mesh.field.setNumber(b, "XMax", wf.x_range[1])
    gmsh.model.mesh.field.setNumber(b, "YMin", wf.y_range[0])
    gmsh.model.mesh.field.setNumber(b, "YMax", wf.y_range[1])
    gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
    gmsh.model.mesh.field.setNumber(b, "ZMax", 1)
    gmsh.model.mesh.field.setNumber(b, "VIn", farmLC)
    gmsh.model.mesh.field.setNumber(b, "VOut", blc)

    mesher = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mesher, "FieldsList", [998, 999])
    gmsh.model.mesh.field.setAsBackgroundMesh(mesher)