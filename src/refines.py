import gmsh

def generateCustomRefines(params):
    """
    Generates all the custom refinements specified and returns them as a list of field
    IDs.

    :param params: The parameter dictionary.
    :type params: dict()
    :return: A list of the fields representing the custom refinements.
    :rtype: list[int]

    """
    n_refines = params['refine_custom']['num_refines']
    dim = params['domain']['dimension']
    blc = params['refine']['background_length_scale']
    lower_aspect = params['domain']['aspect_ratio']
    threshold = params['domain']['aspect_distance']
    upper_aspect = params['domain']['upper_aspect_ratio']

    fields = []
    for i in range(1, n_refines + 1):
        shape = params['refine_custom'][i]['type']
        x = params['refine_custom'][i]['x_range']
        y = params['refine_custom'][i]['y_range']
        z = [0, 1] if dim == 2 else _getAdjustedHeight(lower_aspect, upper_aspect, threshold, params['refine_custom'][i]['z_range'])
        lc = params['refine_custom'][i]['length_scale']
        if shape == 'box':
            fields.append(_customBox(x, y, z, lc, blc))
        elif shape == 'cylinder':
            radius = params['refine_custom'][i]['radius']
            fields.append(_customCylinder(x, y, radius, z, lc, blc))
            
    return fields

def _getAdjustedHeight(lower_aspect, upper_aspect, threshold, z_range):
    """
    Adjusts the lower and upper z-values of a custom refinement in order to account
    for changes caused by the upper or lower aspect ratios.

    :param lower_aspect: The aspect ratio below the threshold distance.
    :type params: double
    :param upper_aspect: The aspect ratio above the threshold distance.
    :type upper_aspect: double
    :param threshold: The threshold that defines the upper and lower aspect split.
    :type threshold: double
    :param z_range: The original, unadjusted z-range of the refinement.
    :type z_range: list[double, double]
    :return: The adjusted heights, in a list.
    :rtype: list[double, double]

    """
    bottom, top = z_range[0], z_range[1]
    if bottom > threshold:
        return [bottom * upper_aspect, top * upper_aspect]
    
    if top <= threshold:
        return [bottom * lower_aspect, top * lower_aspect]
    
    return  [bottom * lower_aspect, top * upper_aspect]

def _customBox(x, y, z, lc, blc):
    """
    Generates a field to represent a single box refinement.

    :param x: The x-range.
    :type x: list[double, double]
    :param y: The y-range.
    :type y: list[double, double]
    :param z: The z-range.
    :type z: list[double, double]
    :param lc: The length scale within the box refinement.
    :type lc: double
    :param blc: The background length scale.
    :type blc: double
    :return: A list of the fields representing the custom refinements.
    :rtype: list[int]

    """
    b = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(b, "XMin", x[0])
    gmsh.model.mesh.field.setNumber(b, "XMax", x[1])
    gmsh.model.mesh.field.setNumber(b, "YMin", y[0])
    gmsh.model.mesh.field.setNumber(b, "YMax", y[1])
    gmsh.model.mesh.field.setNumber(b, "ZMin", z[0])
    gmsh.model.mesh.field.setNumber(b, "ZMax", z[1])
    gmsh.model.mesh.field.setNumber(b, "VIn", lc)
    gmsh.model.mesh.field.setNumber(b, "VOut", blc)

    return b

def _customCylinder(x, y, radius, height, lc, blc):
    """
    Generates a field to represent a single box refinement.

    :param x: The x-coordinate of the center.
    :type x: double
    :param y: The y-coordinate of the center.
    :type y: double
    :param radius: The radius of the cylinder.
    :type radius: double
    :param height: The upper and lower z-coordinates of the cylinder.
    :type height: list[double, double]
    :param lc: The length scale within the box refinement.
    :type lc: double
    :param blc: The background length scale.
    :type blc: double
    :return: A list of the fields representing the custom refinements.
    :rtype: list[int]

    """
    c = gmsh.model.mesh.field.add("Cylinder")

    gmsh.model.mesh.field.setNumber(c, "Radius", radius)
    gmsh.model.mesh.field.setNumber(c, "VIn", lc)
    gmsh.model.mesh.field.setNumber(c, "VOut", blc)
    gmsh.model.mesh.field.setNumber(c, "ZAxis", height[1])
    gmsh.model.mesh.field.setNumber(c, "XCenter", x)
    gmsh.model.mesh.field.setNumber(c, "YCenter", y)
    gmsh.model.mesh.field.setNumber(c, "ZCenter", height[0])

    return c