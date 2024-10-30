import gmsh

def generateCustomRefines(params):
    n_refines = params['refine_custom']['num_refines']
    dim = params['domain']['dimension']
    blc = params['refine']['background_length_scale']

    fields = []
    for i in range(n_refines):
        shape = params['refine_custom'][i]['shape']
        if shape == 'box':
            x = params['refine_custom'][i]['x_range']
            y = params['refine_custom'][i]['y_range']
            z = 0 if dim == 2 else params['refine_custom'][i]['height']
            lc = params['refine_custom'][i]['length_scale']
            fields.append(_customBox(x, y, z, lc, blc))
            
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