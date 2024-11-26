import gmsh
import numpy as np

####
## 2D Meshing Functions
####

def generateTurbine2D(x, y, lcTurbine, lcBackground, dRotor, wf):

    """
    Builds a single turbine in 2D space at the target (x, y) pair. Additionally,
    updates the minimum bounding region representing the farm if necessary.

    :param x: The x coordinate of the turbine center.
    :type x: double
    :param y: The y coordinate of the turbine center.
    :type y: double
    :param lcTurbine: The meshing constraint at the turbine.
    :type lcTurbine: double
    :param lcBackground: The meshing constraint at the background.
    :type lcTurbine: double
    :param dRotor: The rotor distance.
    :type dRotor: double
    :param wf: The structure representing the wind farm.
    :type wf: WindFarm
    :return: The GMESH field representation of the turbine.
    :rtype: int

    """

    points = [(x, y + dRotor), (x, y - dRotor), (x + dRotor, y), (x - dRotor, y)]

    for cords in points:
        wf.updateXMax(cords[0])
        wf.updateXMin(cords[0])
        wf.updateYMax(cords[1])
        wf.updateYMin(cords[1])

    c = gmsh.model.mesh.field.add("Cylinder")

    gmsh.model.mesh.field.setNumber(c, "Radius", dRotor)
    gmsh.model.mesh.field.setNumber(c, "VIn", lcTurbine)
    gmsh.model.mesh.field.setNumber(c, "VOut", lcBackground)
    gmsh.model.mesh.field.setNumber(c, "ZAxis", 1)
    gmsh.model.mesh.field.setNumber(c, "XCenter", x)
    gmsh.model.mesh.field.setNumber(c, "YCenter", y)

    return c

def buildFarms2D(params, wf, domain):

    """
    Builds every turbine in the range [1, num_turbines].
    It generates 2D meshes across the 1D curve loops that represent each turbine.

    :param params: The parameter dictionary.
    :type params: dict()
    :param wf: The structure representing the wind farm.
    :type wf: WindFarm
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: A list of the 2D turbine meshes.
    :rtype: list[int]

    """

    turbines = []

    nFarms = params['refine']['turbine']['num_turbines']
    lcTurbine = params['refine']['turbine']['length_scale'] 
    lcBackground = params['refine']['background_length_scale']
    rotor = params['refine']['turbine']['threshold_rotor_distance']

    for i in range(nFarms):
        turbineData = params['refine']['turbine'][i + 1]
        x = turbineData['x'] 
        y = turbineData['y']

        if not domain.withinDomain(x, y):
            raise Exception("Invalid turbine location.")
        
        turbine = generateTurbine2D(x, y, lcTurbine, lcBackground, rotor, wf)
        turbines.append(turbine)

    return turbines

def refineFarm2D(params, wf):

    """
    Initializes a 'Box' field that sets points within the minimum bounding regoin
    surrounding the farm to the farm's meshing constraint. This field is hard-coded
    to have tag 998.

    :param params: The parameter dictionary.
    :type params: dict()
    :param wf: The structure representing the wind farm.
    :type wf: WindFarm

    """

    dist = params['refine']['farm']['threshold_distance'] 
    farmLC = params['refine']['farm']['length_scale'] 
    blc = params['refine']['background_length_scale']
    ftype = params['refine']['farm']['type']

    wf.adjustDistance(dist)

    if ftype == 'box':
        b = gmsh.model.mesh.field.add("Box", tag=998)
        gmsh.model.mesh.field.setNumber(b, "XMin", wf.x_range[0])
        gmsh.model.mesh.field.setNumber(b, "XMax", wf.x_range[1])
        gmsh.model.mesh.field.setNumber(b, "YMin", wf.y_range[0])
        gmsh.model.mesh.field.setNumber(b, "YMax", wf.y_range[1])
        gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
        gmsh.model.mesh.field.setNumber(b, "ZMax", 1)
        gmsh.model.mesh.field.setNumber(b, "VIn", farmLC)
        gmsh.model.mesh.field.setNumber(b, "VOut", blc)
    else:
        centerX = (wf.x_range[0] + wf.x_range[1]) / 2
        centerY = (wf.y_range[0] + wf.y_range[1]) / 2

        center = np.array([centerX, centerY])
        corner = np.array([wf.x_range[1], wf.y_range[1]])
        radius = np.linalg.norm(center - corner)

        c = gmsh.model.mesh.field.add("Cylinder", tag=998)
        gmsh.model.mesh.field.setNumber(c, "Radius", radius)
        gmsh.model.mesh.field.setNumber(c, "VIn", farmLC)
        gmsh.model.mesh.field.setNumber(c, "VOut", blc)
        gmsh.model.mesh.field.setNumber(c, "ZAxis", 1)
        gmsh.model.mesh.field.setNumber(c, "XCenter", centerX)
        gmsh.model.mesh.field.setNumber(c, "YCenter", centerY)