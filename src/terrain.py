import gmsh
import sys
import yaml
import numpy as np
from scipy.interpolate import LinearNDInterpolator

def buildTerrainFromFile(params, domain):

    """
    A 3D terrain generator used when non-default terrain is specified. This function will
    build an interpolation function using the inputted terrain file and builds a cube where the
    bottom face is determined by the interpolating function. It also sets a 'Box' field with tag 999 
    around the domain that sets the domain-wide meshing constraints.

    :param params: The parameter dictionary.
    :type params: dict()
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: The GMESH surface loop representing the domain.
    :rtype: int

    """

    filename = params['domain']['terrain_path']
    x_range = params['domain']['x_range']
    y_range = params['domain']['y_range']
    height = params['domain']['height']
    aspect = params['domain']['aspect_ratio']
    upper_aspect = params['domain']['upper_aspect_ratio']
    aniso_dist = params['domain']['aspect_distance']
    lc = params['refine']['background_length_scale']

    N = 300

    terrain = np.loadtxt(filename)

    xTerrain = terrain[1:, 0]
    yTerrain = terrain[1:, 1]
    heightTerrain = terrain[1:, 2]

    xMin = x_range[0]
    xMax = x_range[1]
    yMin = y_range[0]
    yMax = y_range[1]

    lowerHeight = aniso_dist * aspect
    upperHeight = (height - aniso_dist) * upper_aspect
    totalHeight = lowerHeight + upperHeight

    b = gmsh.model.mesh.field.add("Box", tag=999)
    gmsh.model.mesh.field.setNumber(b, "XMin", xMin)
    gmsh.model.mesh.field.setNumber(b, "XMax", xMax)
    gmsh.model.mesh.field.setNumber(b, "YMin", yMin)
    gmsh.model.mesh.field.setNumber(b, "YMax", yMax)
    gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
    gmsh.model.mesh.field.setNumber(b, "ZMax", totalHeight)
    gmsh.model.mesh.field.setNumber(b, "VIn", lc)
    gmsh.model.mesh.field.setNumber(b, "VOut", lc * 2)

    domain.setDomain(x_range=[xMin, xMax], y_range=[yMin, yMax], height=totalHeight)

    interp = LinearNDInterpolator(list(zip(xTerrain, yTerrain)), heightTerrain)
    xPoints = np.linspace(xMin, xMax, num=N + 1)
    yPoints = np.linspace(yMin, yMax, num=N + 1)
    domain.setInterp(interp)

    def tag(i, j) -> int:
        return (N + 1) * i + j + 1

    coords = []

    nodes = []

    tris = []

    lin = [[], [], [], []]

    pnt = [tag(0, 0), tag(N, 0), tag(N, N), tag(0, N)]

    for i in range(N + 1):
        for j in range(N + 1):
            nodes.append(tag(i, j))
            coords.extend([
                xPoints[i],
                yPoints[j],
                interp(xPoints[i], yPoints[j]) * aspect
            ])
            if i > 0 and j > 0:
                tris.extend([tag(i - 1, j - 1), tag(i, j - 1), tag(i - 1, j)])
                tris.extend([tag(i, j - 1), tag(i, j), tag(i - 1, j)])
            if (i == 0 or i == N) and j > 0:
                lin[3 if i == 0 else 1].extend([tag(i, j - 1), tag(i, j)])
            if (j == 0 or j == N) and i > 0:
                lin[0 if j == 0 else 2].extend([tag(i - 1, j), tag(i, j)])

    for i in range(4):
        gmsh.model.addDiscreteEntity(0, i + 1)

    gmsh.model.setCoordinates(1, xMin, yMin, coords[3 * tag(0, 0) - 1])
    gmsh.model.setCoordinates(2, xMax, yMin, coords[3 * tag(N, 0) - 1])
    gmsh.model.setCoordinates(3, xMax, yMax, coords[3 * tag(N, N) - 1])
    gmsh.model.setCoordinates(4, xMin, yMax, coords[3 * tag(0, N) - 1])

    for i in range(4):
        gmsh.model.addDiscreteEntity(1, i + 1, [i + 1, i + 2 if i < 3 else 1])

    gmsh.model.addDiscreteEntity(2, 989, [1, 2, 3, 4])

    gmsh.model.mesh.addNodes(2, 989, nodes, coords)

    for i in range(4):
        gmsh.model.mesh.addElementsByType(i + 1, 15, [], [pnt[i]])
        gmsh.model.mesh.addElementsByType(i + 1, 1, [], lin[i])
    gmsh.model.mesh.addElementsByType(989, 2, [], tris)

    gmsh.model.mesh.reclassifyNodes()

    gmsh.model.mesh.createGeometry()

    p1 = gmsh.model.geo.addPoint(xMin, yMin, totalHeight)
    p2 = gmsh.model.geo.addPoint(xMax, yMin, totalHeight)
    p3 = gmsh.model.geo.addPoint(xMax, yMax, totalHeight)
    p4 = gmsh.model.geo.addPoint(xMin, yMax, totalHeight)
    c1 = gmsh.model.geo.addLine(p1, p2)
    c2 = gmsh.model.geo.addLine(p2, p3)
    c3 = gmsh.model.geo.addLine(p3, p4)
    c4 = gmsh.model.geo.addLine(p4, p1)
    c10 = gmsh.model.geo.addLine(p1, 1)
    c11 = gmsh.model.geo.addLine(p2, 2)
    c12 = gmsh.model.geo.addLine(p3, 3)
    c13 = gmsh.model.geo.addLine(p4, 4)
    ll1 = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
    s1 = gmsh.model.geo.addPlaneSurface([ll1], tag=990) #Top Face
    ll3 = gmsh.model.geo.addCurveLoop([c1, c11, -1, -c10]) 
    s3 = gmsh.model.geo.addPlaneSurface([ll3], tag=995) #yMin face
    ll4 = gmsh.model.geo.addCurveLoop([c2, c12, -2, -c11]) 
    s4 = gmsh.model.geo.addPlaneSurface([ll4], tag=994) #xMax face
    ll5 = gmsh.model.geo.addCurveLoop([c3, c13, 3, -c12]) 
    s5 = gmsh.model.geo.addPlaneSurface([ll5], tag=993) #yMax face
    ll6 = gmsh.model.geo.addCurveLoop([c4, c10, 4, -c13]) 
    s6 = gmsh.model.geo.addPlaneSurface([ll6], tag=992) #xMin face
    sl1 = gmsh.model.geo.addSurfaceLoop([s1, s3, s4, s5, s6, 989]) #Bottom face (989)

    p1 = gmsh.model.geo.addPhysicalGroup(2, [994], tag=1)
    p2 = gmsh.model.geo.addPhysicalGroup(2, [993], tag=2)
    p3 = gmsh.model.geo.addPhysicalGroup(2, [992], tag=3)
    p4 = gmsh.model.geo.addPhysicalGroup(2, [995], tag=4)
    p5 = gmsh.model.geo.addPhysicalGroup(2, [989], tag=5)
    p6 = gmsh.model.geo.addPhysicalGroup(2, [990], tag=6)

    gmsh.model.geo.synchronize()

    return sl1

def buildTerrainDefault(params, domain):

    """
    The default 3D terrain generator. Used when no terrain file is specified.
    It will build a cube with the specified dimensions and activate a 'Box' field with tag 999
    setting the domain meshing constraint.

    :param params: The parameter dictionary.
    :type params: dict()
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: The GMESH surface loop representing the domain.
    :rtype: int

    """

    x_range = params['domain']['x_range']
    y_range = params['domain']['y_range']
    height = params['domain']['height']
    aspect = params['domain']['aspect_ratio']
    upper_aspect = params['domain']['upper_aspect_ratio']
    aniso_dist = params['domain']['aspect_distance']
    lc = params['refine']['background_length_scale']

    xMin = x_range[0]
    xMax = x_range[1]
    yMin = y_range[0]
    yMax = y_range[1]
    
    lowerHeight = aniso_dist * aspect
    upperHeight = (height - aniso_dist) * upper_aspect
    totalHeight = lowerHeight + upperHeight
    
    b = gmsh.model.mesh.field.add("Box", tag=999)
    gmsh.model.mesh.field.setNumber(b, "XMin", xMin)
    gmsh.model.mesh.field.setNumber(b, "XMax", xMax)
    gmsh.model.mesh.field.setNumber(b, "YMin", yMin)
    gmsh.model.mesh.field.setNumber(b, "YMax", yMax)
    gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
    gmsh.model.mesh.field.setNumber(b, "ZMax", totalHeight)
    gmsh.model.mesh.field.setNumber(b, "VIn", lc)
    gmsh.model.mesh.field.setNumber(b, "VOut", lc * 2)

    domain.setDomain(x_range=[xMin, xMax], y_range=[yMin, yMax], height=totalHeight)

    b1 = gmsh.model.geo.addPoint(xMax, yMin, 0)
    b2 = gmsh.model.geo.addPoint(xMin, yMin, 0)
    b3 = gmsh.model.geo.addPoint(xMax, yMax, 0)
    b4 = gmsh.model.geo.addPoint(xMin, yMax, 0)

    lb1 = gmsh.model.geo.addLine(b1, b3)
    lb2 = gmsh.model.geo.addLine(b3, b4)
    lb3 = gmsh.model.geo.addLine(b4, b2)
    lb4 =gmsh.model.geo.addLine(b2, b1)

    farmBase = gmsh.model.geo.addCurveLoop([lb1, lb2, lb3, lb4])
    base = gmsh.model.geo.addPlaneSurface([farmBase], tag=989)

    b5 = gmsh.model.geo.addPoint(xMax, yMin, totalHeight)
    b6 = gmsh.model.geo.addPoint(xMin, yMin, totalHeight)
    b7 = gmsh.model.geo.addPoint(xMax, yMax, totalHeight)
    b8 = gmsh.model.geo.addPoint(xMin, yMax, totalHeight)

    lt1 = gmsh.model.geo.addLine(b5, b7)
    lt2 = gmsh.model.geo.addLine(b7, b8)
    lt3 = gmsh.model.geo.addLine(b8, b6)
    lt4 =gmsh.model.geo.addLine(b6, b5)

    farmTop = gmsh.model.geo.addCurveLoop([lt1, lt2, lt3, lt4])
    top = gmsh.model.geo.addPlaneSurface([farmTop], tag=990)

    lc1 = gmsh.model.geo.addLine(b1, b5)
    lc2 = gmsh.model.geo.addLine(b2, b6)
    lc3 = gmsh.model.geo.addLine(b3, b7)
    lc4 = gmsh.model.geo.addLine(b4, b8)

    face1 = gmsh.model.geo.addCurveLoop([lc1, lt4, lc2, lb4], reorient=True) #YMin Face
    face2 = gmsh.model.geo.addCurveLoop([lc2, lt3, lc4, lb3], reorient=True) #XMin Face
    face3 = gmsh.model.geo.addCurveLoop([lc4, lt2, lc3, lb2], reorient=True) #YMax Face
    face4 = gmsh.model.geo.addCurveLoop([lc3, lt1, lc1, lb1], reorient=True) #XMax Face

    f1 = gmsh.model.geo.addPlaneSurface([face1], tag=995)
    f2 = gmsh.model.geo.addPlaneSurface([face2], tag=992)
    f3 = gmsh.model.geo.addPlaneSurface([face3], tag=993)
    f4 = gmsh.model.geo.addPlaneSurface([face4], tag=994)

    s1 = gmsh.model.geo.addPhysicalGroup(2, [994], tag=1)
    s2 = gmsh.model.geo.addPhysicalGroup(2, [993], tag=2)
    s3 = gmsh.model.geo.addPhysicalGroup(2, [992], tag=3)
    s4 = gmsh.model.geo.addPhysicalGroup(2, [995], tag=4)
    s5 = gmsh.model.geo.addPhysicalGroup(2, [989], tag=5)
    s6 = gmsh.model.geo.addPhysicalGroup(2, [990], tag=6)

    return gmsh.model.geo.addSurfaceLoop([base, f1, f2, f3, f4, top])

def buildTerrainCylinder(params, domain):
    centerLocation = params['domain']['center']
    radius = params['domain']['radius']
    height = params['domain']['height']
    aspect = params['domain']['aspect_ratio']
    upper_aspect = params['domain']['upper_aspect_ratio']
    aspect_distance = params['domain']['aspect_distance']
    lc = params['refine']['background_length_scale']
    filename = params['domain'].get('terrain_path')

    lowerHeight = aspect_distance * aspect
    upperHeight = (height - aspect_distance) * upper_aspect
    totalHeight = lowerHeight + upperHeight

    domain.setDomain(radius=radius, center=centerLocation, height=totalHeight)

    if filename:
        terrain = np.loadtxt(filename)

        xTerrain = terrain[1:, 0]
        yTerrain = terrain[1:, 1]
        heightTerrain = terrain[1:, 2]

        interp = LinearNDInterpolator(list(zip(xTerrain, yTerrain)), heightTerrain)
        domain.setInterp(interp)

    c = gmsh.model.mesh.field.add("Cylinder", tag=999)
    gmsh.model.mesh.field.setNumber(c, "Radius", radius)
    gmsh.model.mesh.field.setNumber(c, "XCenter", centerLocation[0])
    gmsh.model.mesh.field.setNumber(c, "YCenter", centerLocation[1])
    gmsh.model.mesh.field.setNumber(c, "ZAxis", height)
    gmsh.model.mesh.field.setNumber(c, "VIn", lc)
    gmsh.model.mesh.field.setNumber(c, "VOut", lc * 2)

    center = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1], 0)
    posX = gmsh.model.geo.addPoint(centerLocation[0] + radius, centerLocation[1], 0)
    negX = gmsh.model.geo.addPoint(centerLocation[0] - radius, centerLocation[1], 0)
    posY = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1] + radius, 0)
    negY = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1] - radius, 0)

    a1 = gmsh.model.geo.addCircleArc(posX, center, negY)
    a2 = gmsh.model.geo.addCircleArc(negY, center, negX)
    a3 = gmsh.model.geo.addCircleArc(negX, center, posY)
    a4 = gmsh.model.geo.addCircleArc(posY, center, posX)

    loop = gmsh.model.geo.addCurveLoop([a1, a2, a3, a4])
    surface = gmsh.model.geo.addPlaneSurface([loop], tag=999)

    height = params['domain']['height']
    v1 = gmsh.model.geo.extrude([(2, surface)], 0, 0, totalHeight)

    base = gmsh.model.geo.addPhysicalGroup(2, [999], tag=7)
    top = gmsh.model.geo.addPhysicalGroup(2, [1021], tag=8)
    outflow = gmsh.model.geo.addPhysicalGroup(2, [1020, 1008], tag=5)
    inflow = gmsh.model.geo.addPhysicalGroup(2, [1016, 1012], tag=6)

    vol = gmsh.model.geo.addPhysicalGroup(3, [1], tag=0)

def buildTerrainCircle(params, domain):
    centerLocation = params['domain']['center']
    radius = params['domain']['radius']
    lc = params['refine']['background_length_scale']

    domain.setDomain(radius=radius, center=centerLocation)

    c = gmsh.model.mesh.field.add("Cylinder", tag=999)
    gmsh.model.mesh.field.setNumber(c, "Radius", radius)
    gmsh.model.mesh.field.setNumber(c, "XCenter", centerLocation[0])
    gmsh.model.mesh.field.setNumber(c, "YCenter", centerLocation[1])
    gmsh.model.mesh.field.setNumber(c, "ZAxis", 1)
    gmsh.model.mesh.field.setNumber(c, "VIn", lc)
    gmsh.model.mesh.field.setNumber(c, "VOut", lc * 2)

    center = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1], 0)
    posX = gmsh.model.geo.addPoint(centerLocation[0] + radius, centerLocation[1], 0)
    negX = gmsh.model.geo.addPoint(centerLocation[0] - radius, centerLocation[1], 0)
    posY = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1] + radius, 0)
    negY = gmsh.model.geo.addPoint(centerLocation[0], centerLocation[1] - radius, 0)

    a1 = gmsh.model.geo.addCircleArc(posX, center, negY)
    a2 = gmsh.model.geo.addCircleArc(negY, center, negX)
    a3 = gmsh.model.geo.addCircleArc(negX, center, posY)
    a4 = gmsh.model.geo.addCircleArc(posY, center, posX)

    outflow = gmsh.model.geo.addPhysicalGroup(1, [a1, a4], tag=7)
    inflow = gmsh.model.geo.addPhysicalGroup(1, [a2, a3], tag=8)

    loop = gmsh.model.geo.addCurveLoop([a1, a2, a3, a4])
    
    return loop

def buildTerrain2D(params, domain):

    """
    The default 2D terrain generator. It will build a rectangular region
    as specified by the parameters. It will also initialize a 'Box' field
    with tag 999 to enforce the domain-wide meshing constraints.

    :param params: The parameter dictionary.
    :type params: dict()
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: The GMESH curve loop representing the domain.
    :rtype: int

    """

    x_range = params['domain']['x_range']
    y_range = params['domain']['y_range']
    lc = params['refine']['background_length_scale']
    
    domain.setDomain(x_range=x_range, y_range=y_range, height=0)

    b1 = gmsh.model.geo.addPoint(x_range[1], y_range[0], 0)
    b2 = gmsh.model.geo.addPoint(x_range[0], y_range[0], 0)
    b3 = gmsh.model.geo.addPoint(x_range[1], y_range[1], 0)
    b4 = gmsh.model.geo.addPoint(x_range[0], y_range[1], 0)

    lb1 = gmsh.model.geo.addLine(b1, b3, tag=994) #XMax
    lb2 = gmsh.model.geo.addLine(b3, b4, tag=993) #YMax
    lb3 = gmsh.model.geo.addLine(b4, b2, tag=992) #XMin
    lb4 =gmsh.model.geo.addLine(b2, b1, tag=995) #YMin

    farmBorder = gmsh.model.geo.addCurveLoop([lb1, lb2, lb3, lb4])

    s1 = gmsh.model.geo.addPhysicalGroup(1, [994], tag=1)
    s2 = gmsh.model.geo.addPhysicalGroup(1, [993], tag=2)
    s3 = gmsh.model.geo.addPhysicalGroup(1, [992], tag=3)
    s4 = gmsh.model.geo.addPhysicalGroup(1, [995], tag=4)

    b = gmsh.model.mesh.field.add("Box", tag=999)
    gmsh.model.mesh.field.setNumber(b, "XMin", x_range[0])
    gmsh.model.mesh.field.setNumber(b, "XMax", x_range[1])
    gmsh.model.mesh.field.setNumber(b, "YMin", y_range[0])
    gmsh.model.mesh.field.setNumber(b, "YMax", y_range[1])
    gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
    gmsh.model.mesh.field.setNumber(b, "ZMax", 1)
    gmsh.model.mesh.field.setNumber(b, "VIn", lc)
    gmsh.model.mesh.field.setNumber(b, "VOut", lc * 2)

    return farmBorder