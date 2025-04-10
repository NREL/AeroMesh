import gmsh
import math
import numpy as np

def generateTurbines(params, domain, wf):

    """
    Builds every turbine in the range [1, num_turbines].
    It also adjusts the minimum bounding regions of the wind farm to surround
    the turbines.

    :param params: The parameter dictionary.
    :type params: dict()
    :param domain: The structure representing the domain.
    :type domain: Domain
    :param wf: The structure representing the wind farm.
    :type wf: WindFarm
    :return: A list of all the GMESH fields that represent the turbines.
    :rtype: list[int]

    """

    fields = []
    nFarms = params['refine']['turbine']['num_turbines']
    lc = params['refine']['turbine']['length_scale'] 
    lcb = params['refine']['background_length_scale']
    lcf = params['refine']['farm']['length_scale']
    rotor = params['refine']['turbine']['threshold_rotor_distance']
    aspect = params['domain']['aspect_ratio']
    aspect_upper = params['domain']['upper_aspect_ratio']
    aspect_distance = params['domain']['aspect_distance']
    inflow = params['domain']['inflow_angle']
    turbineType = params['refine']['turbine']['type']

    interp = domain.interp
    for i in range(nFarms):
        turbineData = params['refine']['turbine'][i + 1]
        x = turbineData['x'] 
        y = turbineData['y']
        hh = 100
        if turbineData.get('HH') is not None:
            hh = turbineData.get('HH')
        if interp is not None:
            z = (interp(x, y) + hh) * aspect
        else:
            z = hh * aspect

        if not domain.withinDomain(x, y, z):
            raise Exception("Invalid turbine location.")
        
        if turbineType == 'wake':
            upstream = params['refine']['turbine']['threshold_upstream_distance']
            downstream = params['refine']['turbine']['threshold_downstream_distance']
            fields.extend(placeTurbineWake(x, y, z, upstream, downstream, rotor, lc, lcb, lcf, inflow, aspect, wf, domain))
        elif turbineType == 'simple':
            ground = domain.calculateGround(x, y)
            fields.extend(placeTurbineSimple(x, y, rotor, hh, lc, lcb, ground, aspect, aspect_upper, aspect_distance, wf))
        elif turbineType == 'sphere':
            fields.extend(placeTurbineSphere(x, y, z, rotor, lc, lcb, wf))
        else:
            raise Exception("Invalid turbine type. Must be simple, wake, or sphere.")

    return fields

def placeTurbineSimple(x, y, rotor, hh, lc, lcb, ground, aspect, upperAspect, aspectDist, wf):
    radius = rotor / 2
    top = ground + hh + (rotor / 2)
    bottom = ground
    bottomDist = aspectDist - bottom
    topDist = top - aspectDist

    wf.updateXMax(x)
    wf.updateXMin(x)
    wf.updateYMax(y)
    wf.updateYMin(y)
    wf.updateZMax(top + (bottomDist * (aspect - 1)) - (topDist * (1 - upperAspect)))
    
    c = gmsh.model.mesh.field.add("Cylinder")
    gmsh.model.mesh.field.setNumber(c, "Radius", radius)
    gmsh.model.mesh.field.setNumber(c, "XCenter", x)
    gmsh.model.mesh.field.setNumber(c, "YCenter", y)
    gmsh.model.mesh.field.setNumber(c, "ZAxis", top + (bottomDist * (aspect - 1)) - (topDist * (1 - upperAspect)))
    gmsh.model.mesh.field.setNumber(c, "VIn", lc)
    gmsh.model.mesh.field.setNumber(c, "VOut", lcb)

    fields = [c]
    return fields

def placeTurbineSphere(x, y, hh, rotor, lc, lcb, wf):
    wf.updateXMax(x)
    wf.updateXMin(x)
    wf.updateYMax(y)
    wf.updateYMin(y)
    wf.updateZMax(hh)

    s = gmsh.model.mesh.field.add("Ball")
    radius = rotor / 2

    gmsh.model.mesh.field.setNumber(s, "Radius", radius)
    gmsh.model.mesh.field.setNumber(s, "XCenter", x)
    gmsh.model.mesh.field.setNumber(s, "YCenter", y)
    gmsh.model.mesh.field.setNumber(s, "ZCenter", hh)
    gmsh.model.mesh.field.setNumber(s, "VIn", lc)
    gmsh.model.mesh.field.setNumber(s, "VOut", lcb)

    return s
        
def placeTurbineWake(x, y, z, upstream, downstream, rotor, lc, lcb, lcf, inflow, aspect, wf, domain):

    """
    Builds a single turbine in 3D space at the target (x, y, z) triple.
    Updates the minimum bounding region representing the farm if necessary.
    Adjusts the meshing region in the z-direction to account for anisotropic effects if needed.

    :param x: The x coordinate of the turbine center.
    :type x: double
    :param y: The y coordinate of the turbine center.
    :type y: double
    :param z: The z coordinate of the turbine center.
    :type z: double
    :param upstream: The extension of the wake in the negative direction.
    :type upstream: double
    :param downstream: The extension of the wake in the positive direction.
    :type downstream: double
    :param rotor: The radius of the rotor.
    :type rotor: double
    :param lc: The turbine level meshing restriction.
    :type lc: double
    :param lcb: The domain level meshing restriction.
    :type lcb: double
    :param lcf: The farm level meshing restriction.
    :type lcf: double
    :param inflow: The incoming wind direction.
    :type inflow: double
    :param aspect: The aspect ratio.
    :type aspect: int
    :param wf: The structure representing the wind farm.
    :type wf: WindFfarm
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: A list of GMESH fields that define the rotor's meshing region plus any anisotropic meshing in the z-direction.
    :rtype: list[int]

    """

    ###
    # Initialize proper data structures.
    ###

    increment = rotor / 2
    downPoints = math.ceil(downstream / increment)
    upPoints = math.ceil(upstream / increment)
    turbine = [gmsh.model.geo.addPoint(x, y, z)]
    anisoPoints = []

    ###
    # These loops build the body of the turbine.
    ###
    for i in range(1, downPoints + 1):
        turbine.append(gmsh.model.geo.addPoint(x + increment * i, y, z))
   
    for i in range(1, upPoints + 1):
        turbine.append(gmsh.model.geo.addPoint(x - increment * i, y, z))
    
    
    if aspect == 1:
        aspect = 0 # Locally turns off anisotropy loops.

    ###
    # These loops add the points needed to maintain the proper aspect ratio.
    ###
    for i in range(1, math.floor(aspect + 1)):
        level = []
        if z + (rotor * i) < domain.z_range[1]: 
            gmsh.model.geo.addPoint(x, y, z + (rotor * i))
        if z - (rotor * i) > domain.calculateGround(x, y):
            gmsh.model.geo.addPoint(x, y, z - (rotor * i))
        for j in range(1, downPoints + 1):
            if z + (rotor * i) < domain.z_range[1]: 
                level.append(gmsh.model.geo.addPoint(x + increment * j, y, z + (rotor * i)))
            if z - (rotor * i) > domain.calculateGround(x + increment * j, y):
                level.append(gmsh.model.geo.addPoint(x + increment * j, y, z - (rotor * i)))

        for j in range(1, upPoints + 1):
            if z + (rotor * i) < domain.z_range[1]: 
                level.append(gmsh.model.geo.addPoint(x - increment * j, y, z + (rotor * i)))
            if z - (rotor * i) > domain.calculateGround(x - increment * j, y):
                level.append(gmsh.model.geo.addPoint(x - increment * j, y, z - (rotor * i)))
        
        anisoPoints.append(level)

    ###
    # Embed all the points into the proper surfaces and volumes.
    ###

    turbineTags = list(tag for tag in zip([0] * len(turbine), turbine))
    gmsh.model.geo.rotate(turbineTags, x, y, z, 0, 0, 1, inflow)
    gmsh.model.geo.synchronize()

    for point in turbine:
        coords = gmsh.model.getValue(0, point, [])
        if domain.withinDomain(coords[0], coords[1], coords[2]):
            wf.updateXMax(coords[0])
            wf.updateXMin(coords[0])
            wf.updateYMax(coords[1])
            wf.updateYMin(coords[1])
            wf.updateZMax(coords[2])
        else:
            gmsh.model.geo.remove([(0, point)])
            turbine.remove(point)
    gmsh.model.geo.synchronize()

    for level in anisoPoints:
        levelTags = list(tag for tag in zip([0] * len(level), level))
        gmsh.model.geo.rotate(levelTags, x, y, z, 0, 0, 1, inflow)

        gmsh.model.geo.synchronize()

        for point in level:
            coords = gmsh.model.getValue(0, point, [])
            if domain.withinDomain(coords[0], coords[1], coords[2]):
                wf.updateXMax(coords[0])
                wf.updateXMin(coords[0])
                wf.updateYMax(coords[1])
                wf.updateYMin(coords[1])
                wf.updateZMax(coords[2])
            else:
                gmsh.model.geo.remove([(0, point)])
                level.remove(point)

        gmsh.model.geo.synchronize()

    gmsh.model.geo.synchronize()

    ###
    # Initialize the distance fields.
    ###
    f = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f, "PointsList", turbine)

    t = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(t, "InField", f)
    gmsh.model.mesh.field.setNumber(t, "SizeMin", lc)
    gmsh.model.mesh.field.setNumber(t, "SizeMax", lcb)
    gmsh.model.mesh.field.setNumber(t, "DistMin", rotor)
    gmsh.model.mesh.field.setNumber(t, "DistMax", rotor + 0.5 * (lc + lcf) * 4)

    fields = [t]

    anisoScale = (aspect - 1) / 2
    a = rotor
    b = rotor + (rotor * 2 * anisoScale)

    for i, level in enumerate(anisoPoints):
        fa = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fa, "PointsList", level)

        zScaled = rotor * (i + 1) #Shift so that z is centered at 0.
        
        anisoRadius = calcEllipse(a, b, zScaled)
    
        ta = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(ta, "InField", fa)
        gmsh.model.mesh.field.setNumber(ta, "SizeMin", lc)
        gmsh.model.mesh.field.setNumber(ta, "SizeMax", lcb)
        gmsh.model.mesh.field.setNumber(ta, "DistMin", anisoRadius)
        gmsh.model.mesh.field.setNumber(ta, "DistMax", anisoRadius + 0.5 * (lc + lcf) * 4)

        fields.append(ta)

    return fields

def anisotropyScale(params):

    """
    Compresses the mesh nodes to generate the anisotropic effects.

    :param params: The parameter dictionary.
    :type params: dict()
    """

    aspect = params['domain']['aspect_ratio']
    upper_aspect = params['domain']['upper_aspect_ratio']
    z_range = params['domain']['z_range']
    base, height = z_range[0], z_range[1]
    dist = params['domain']['aspect_distance']

    if aspect == 1:
        return

    data = gmsh.model.mesh.getNodes()
    tags = data[0]
    coords = data[1].reshape(-1, 3)

    for node in zip(tags, coords):
        tag = node[0]
        coord = node[1]
        if coord[2] <= dist * aspect:
            coord[2] /= aspect
        else:
            coord[2] -= (dist * (aspect))
            coord[2] /= ((height - dist) * upper_aspect)
            coord[2] *= (height - dist)
            coord[2] += dist
        
        gmsh.model.mesh.setNode(tag, coord, [])

def calcEllipse(a, b, z):

    """
    Calculates the distance to the border of the radius of the
    ellipse created by any anisotropic effects.

    :param a: The ellipse minor axis.
    :type a: double
    :param b: The ellipse major axis.
    :type b: double
    :param z: The elevation of the ellipse formed by the turbine, scaled so the center is at 0.
    :type z: double
    :return: The distance to the border of the anisotropic ellipse.
    :rtype: double
    """

    term = (z ** 2) / (b ** 2)
    delta = 1 - term
    delta *= (a ** 2)
    return math.sqrt(delta)

def refineFarm3D(params, wf):

    """
    Initializes a 'Box' or 'Cylinder' field that sets points within the minimum bounding regoin
    surrounding the farm to the farm's meshing constraint.

    :param params: The parameter dictionary.
    :type params: dict()
    :param wf: The structure representing the wind farm.
    :type wf: WindFarm
    :return: The GMESH field that represents the box.
    :rtype: int
    """

    lc = params['refine']['farm']['length_scale']
    jitter = params['refine']['farm']['threshold_distance']
    lcb = params['refine']['background_length_scale']
    farmType = params['refine']['farm']['type']
    
    if farmType == 'box':
        wf.adjustDistance(jitter)
        b = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(b, "XMin", wf.x_range[0])
        gmsh.model.mesh.field.setNumber(b, "XMax", wf.x_range[1])
        gmsh.model.mesh.field.setNumber(b, "YMin", wf.y_range[0])
        gmsh.model.mesh.field.setNumber(b, "YMax", wf.y_range[1])
        gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
        gmsh.model.mesh.field.setNumber(b, "ZMax", wf.zMax)
        gmsh.model.mesh.field.setNumber(b, "VIn", lc)
        gmsh.model.mesh.field.setNumber(b, "VOut", lcb)

        return b
    elif farmType == 'cylinder':
        wf.adjustDistance(jitter)
        centerX = (wf.x_range[0] + wf.x_range[1]) / 2
        centerY = (wf.y_range[0] + wf.y_range[1]) / 2

        center = np.array([centerX, centerY])
        corner = np.array([wf.x_range[1], wf.y_range[1]])
        radius = np.linalg.norm(center - corner)

        c = gmsh.model.mesh.field.add("Cylinder")
        gmsh.model.mesh.field.setNumber(c, "Radius", radius)
        gmsh.model.mesh.field.setNumber(c, "VIn", lc)
        gmsh.model.mesh.field.setNumber(c, "VOut", lcb)
        gmsh.model.mesh.field.setNumber(c, "ZAxis", wf.zMax)
        gmsh.model.mesh.field.setNumber(c, "XCenter", centerX)
        gmsh.model.mesh.field.setNumber(c, "YCenter", centerY)

        return c
    else:
        return None
    
def cylinderTerrainAdjustment(domain, params):

    z_range = params['domain']['z_range']
    gradient = lambda z, h : (h - z) / h
    base, height = z_range[0], z_range[1]

    nodes = gmsh.model.mesh.getNodes()
    tags = nodes[0]
    coords = nodes[1].reshape(-1, 3)
    for tag, coord in zip(tags, coords):
        x, y, z = coord[0], coord[1], coord[2]
        interpolation = domain.interp(x, y)
        coord[2] = (interpolation * gradient(z, height)) + z

        gmsh.model.mesh.setNode(tag, coord, [])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.removeDuplicateElements()
    gmsh.model.mesh.removeDuplicateNodes()