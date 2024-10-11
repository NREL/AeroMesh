import gmsh
import math

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
    upstream = params['refine']['turbine']['threshold_upstream_distance']
    downstream = params['refine']['turbine']['threshold_downstream_distance']
    rotor = params['refine']['turbine']['threshold_rotor_distance']
    aspect = params['domain']['aspect_ratio']

    interp = domain.interp
    for i in range(nFarms):
        turbineData = params['refine']['turbine'][i + 1]
        x = turbineData['x'] 
        y = turbineData['y'] 
        if interp is not None:
            z = (interp(x, y) + 100) * aspect
        else:
            z = 100

        if not domain.withinDomain(x, y, z):
            raise Exception("Invalid turbine location: (" + str(x) + ", " + str(y) + ", " + str(z)) + ")."
        
        # 0: points placed in x-direction 1: points placed in y-direction.
        wake = turbineData['wake']
        eUpstream, eDownstream = checkBounds(domain, upstream, downstream, (x, y), wake)
        fields.extend(placeTurbine(x, y, z, eUpstream, eDownstream, rotor, lc, lcb, lcf, wake, aspect, domain))
        if wake == 0:
            wf.updateXMax(x + eDownstream)
            wf.updateXMin(x - eUpstream)
            wf.updateYMax(y)
            wf.updateYMin(y)
        else:
            wf.updateYMin(y - eUpstream)
            wf.updateYMax(y + eDownstream)
            wf.updateXMax(x)
            wf.updateXMin(x)
        wf.updateZMax(z + rotor)
    return fields
        
def placeTurbine(x, y, z, upstream, downstream, rotor, lc, lcb, lcf, wake, aspect, domain):

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
    :param wake: The direction of the wake.
    :type wake: int
    :param aspect: The aspect ratio.
    :type aspect: int
    :param domain: The structure representing the domain.
    :type domain: Domain
    :return: A list of GMESH fields that define the rotor's meshing region plus any anisotropic meshing in the z-direction.
    :rtype: list[int]

    """

    ###
    # Initialize proper data structures.
    ###
    goY = 0
    goX = 0
    if wake == 0:
        goX = 1
    else:
        goY = 1

    increment = rotor / 2
    downPoints = math.ceil(downstream / increment)
    upPoints = math.ceil(upstream / increment)
    turbine = [gmsh.model.geo.addPoint(x, y, z)]
    anisoPoints = []
    xMinB = []
    xMaxB = []
    yMinB = []
    yMaxB = []


    ###
    # Check boundary conditions.
    ###
    xMin, xMax, yMin, yMax = 0, 0, 0, 0

    if x + goX * increment * downPoints > domain.x_range[1]:
        downPoints -= 1
        xMax = 1
    if y + goY * increment * downPoints > domain.y_range[1]:
        downPoints -= 1
        yMax = 1
    if x - goX * increment * upPoints < domain.x_range[0]:
        upPoints -= 1
        xMin = 1
    if y - goX * increment * upPoints < domain.y_range[0]:
        upPoints -= 1
        yMin = 1


    ###
    # These loops build the body of the turbine.
    ###
    for i in range(1, downPoints + 1):
        turbine.append(gmsh.model.geo.addPoint(x + goX * increment * i, y + goY * increment * i, z))
    
    if xMax == 1:
        point = gmsh.model.geo.addPoint(domain.x_range[1], y, z)
        turbine.append(point)
        xMaxB.append(point)
    if yMax == 1:
        point = gmsh.model.geo.addPoint(x, domain.y_range[1], z)
        turbine.append(point)
        yMaxB.append(point)
   
    for i in range(1, upPoints + 1):
        turbine.append(gmsh.model.geo.addPoint(x - goX * increment * i, y - goY * increment * i, z))
    
    if xMin == 1:
        point = gmsh.model.geo.addPoint(domain.x_range[0], y, z)
        turbine.append(point)
        xMinB.append(point)
    if yMin == 1:
        point = gmsh.model.geo.addPoint(x, domain.y_range[0], z)
        turbine.append(point)
        yMinB.append(point)

    if aspect == 1:
        aspect = 0 # Locally turns off anisotropy loops.

    ###
    # These loops add the points needed to maintain the proper aspect ratio.
    ###
    for i in range(1, aspect + 1):
        level = [gmsh.model.geo.addPoint(x, y, z + (rotor * i)), gmsh.model.geo.addPoint(x, y, z - (rotor * i))]
        for j in range(1, downPoints + 1):
            level.append(gmsh.model.geo.addPoint(x + goX * increment * j, y + goY * increment * j, z + (rotor * i)))
            level.append(gmsh.model.geo.addPoint(x + goX * increment * j, y + goY * increment * j, z - (rotor * i)))
        
        if xMax == 1:
            upper = gmsh.model.geo.addPoint(domain.x_range[1], y, z + (rotor * i))
            lower = gmsh.model.geo.addPoint(domain.x_range[1], y, z - (rotor * i))
            level.append(upper)
            level.append(lower)
            xMaxB.append(upper)
            xMaxB.append(lower)
        if yMax == 1:
            upper = gmsh.model.geo.addPoint(x, domain.y_range[1], z + (rotor * i))
            lower = gmsh.model.geo.addPoint(x, domain.y_range[1], z - (rotor * i))
            level.append(upper)
            level.append(lower)
            yMaxB.append(upper)
            yMaxB.append(lower)

        for j in range(1, upPoints + 1):
            level.append(gmsh.model.geo.addPoint(x - goX * increment * j, y - goY * increment * j, z + (rotor * i)))
            level.append(gmsh.model.geo.addPoint(x - goX * increment * j, y - goY * increment * j, z - (rotor * i)))

        if xMin == 1:
            upper = gmsh.model.geo.addPoint(domain.x_range[0], y, z + (rotor * i))
            lower = gmsh.model.geo.addPoint(domain.x_range[0], y, z - (rotor * i))
            level.append(upper)
            level.append(lower)
            xMinB.append(upper)
            xMinB.append(lower)
        if yMin == 1:
            upper = gmsh.model.geo.addPoint(x, domain.y_range[0], z + (rotor * i))
            lower = gmsh.model.geo.addPoint(x, domain.y_range[0], z - (rotor * i))
            level.append(upper)
            level.append(lower)
            yMinB.append(upper)
            yMinB.append(lower)
        
        anisoPoints.append(level)

    ###
    # Embed all the points into the proper surfaces and volumes.
    ###
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, turbine, 3, 999)

    for level in anisoPoints:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, level, 3, 999)
    gmsh.model.geo.synchronize()

    if len(xMinB) > 0:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, xMinB, 2, 992)
        gmsh.model.geo.synchronize()

    if len(xMaxB) > 0:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, xMaxB, 2, 994)
        gmsh.model.geo.synchronize()
    
    if len(yMinB) > 0:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, yMinB, 2, 995)
        gmsh.model.geo.synchronize()

    if len(yMaxB) > 0:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, yMaxB, 2, 993)
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

def checkBounds(domain, upstream, downstream, origin, wake):

    """
    Checks if the extension of a turbine's wake interests the borders of the domain.
    If this case occurs, the wake is stopped at the domain boundary.

    :param domain: The structure representing the domain.
    :type domain: Domain
    :param upstream: The extension of the wake in the negative direction.
    :type upstream: double
    :param downstream: The extension of the wake in the positive direction.
    :type upstream: double
    :param origin: The turbine center.
    :type origin: (double, double)
    :param wake: The direction of the wake.
    :type wake: int
    :return: An adjusted upstream and downstream distance that prevent domain violations.
    :rtype: (double, double)

    """

    eUpstream = upstream
    eDownstream = downstream
    if wake == 0:
        if origin[0] - upstream < domain.x_range[0]:
            eUpstream = origin[0] - domain.x_range[0]
        if origin[0] + downstream > domain.x_range[1]:
            eDownstream = domain.x_range[1] - origin[0]
    else:
        if origin[1] - upstream < domain.y_range[0]:
            eUpstream = origin[1] - domain.y_range[0]
        if origin[1] + downstream > domain.y_range[1]:
            eDownstream = domain.y_range[1] - origin[1]
    return eUpstream, eDownstream

def anisotropyScale(params):

    """
    Compresses the mesh nodes to generate the anisotropic effects.

    :param params: The parameter dictionary.
    :type params: dict()
    """

    aspect = params['domain']['aspect_ratio']
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
            coord[2] -= dist
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
    Initializes a 'Box' field that sets points within the minimum bounding regoin
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