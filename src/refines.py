import gmsh

def generateCustomRefines(params):
    n_refines = params['refine_custom']['num_refines']
    dim = params['domain']['dimension']
    blc = params['refine']['background_length_scale']

    fields = []
    for i in range(1, n_refines + 1):
        shape = params['refine_custom'][i]['type']
        if shape == 'box':
            x = params['refine_custom'][i]['x_range']
            y = params['refine_custom'][i]['y_range']
            z = 0 if dim == 2 else params['refine_custom'][i]['height']
            lc = params['refine_custom'][i]['length_scale']
            fields.append(_customBox(x, y, z, lc, blc))
        elif shape == 'cylinder':
            x = params['refine_custom'][i]['x_range']
            y = params['refine_custom'][i]['y_range']
            height = 0 if dim == 2 else params['refine_custom'][i]['height']
            lc = params['refine_custom'][i]['length_scale']
            radius = params['refine_custom'][i]['radius']
            fields.append(_customCylinder(x, y, radius, height, lc, blc))
            
    return fields


def _customBox(x, y, z, lc, blc):
    b = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(b, "XMin", x[0])
    gmsh.model.mesh.field.setNumber(b, "XMax", x[1])
    gmsh.model.mesh.field.setNumber(b, "YMin", y[0])
    gmsh.model.mesh.field.setNumber(b, "YMax", y[1])
    gmsh.model.mesh.field.setNumber(b, "ZMin", 0)
    gmsh.model.mesh.field.setNumber(b, "ZMax", z)
    gmsh.model.mesh.field.setNumber(b, "VIn", lc)
    gmsh.model.mesh.field.setNumber(b, "VOut", blc)

    return b

def _customCylinder(x, y, radius, height, lc, blc):
    c = gmsh.model.mesh.field.add("Cylinder")

    gmsh.model.mesh.field.setNumber(c, "Radius", radius)
    gmsh.model.mesh.field.setNumber(c, "VIn", lc)
    gmsh.model.mesh.field.setNumber(c, "VOut", blc)
    gmsh.model.mesh.field.setNumber(c, "ZAxis", height)
    gmsh.model.mesh.field.setNumber(c, "XCenter", x)
    gmsh.model.mesh.field.setNumber(c, "YCenter", y)

    return c