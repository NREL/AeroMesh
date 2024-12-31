Parameters
====================================

.. _yaml_params:

YAML and Directory Structure
-----------------------------

The YAML structure below details all possible parameters, their expected data types, and input formats.

.. code-block:: yaml
    :linenos:

    filetype: string
    domain:
      terrain_path: string
      x_range: [double, double]
      y_range: [double, double]
      height: double
      aspect_ratio: int
      upper_aspect_ratio: double
      aspect_distance: double
      dimension: int
      inflow_angle: double
      type: string
      center: [double, double]
      radius: double
    refine:
      global_scale: double
      background_length_scale: double
      farm:
        type: string
        length_scale: double
        threshold_distance: double
      turbine:
        type: string
        length_scale: double
        threshold_upstream_distance: double
        threshold_downstream_distance: double
        threshold_rotor_distance: double
        num_turbines: int
        1:
          x: double
          y: double
          HH: double
        2:
          x: double
          y: double
          HH: double

          # ...

        n:
          x: double
          y: double
          HH: double
    refine_custom:
      num_refines: int
      1:
        type: string
        x_range: [double, double] or double
        y_range: [double, double] or double
        z_range: [double, double]
        radius: double
        length_scale: double

        # ...

      n:
        type: string
        x_range: [double, double] or double
        y_range: [double, double] or double
        z_range: [double, double]
        radius: double
        length_scale: double

Parameter Descriptions
-----------------------------

The schema above contains a combination of 3D, 2D, and general parameters.
Additionally, some of these parameters are optional. The tables below will provide further
details on the parameters and their functions.

2D and 3D Parameters
~~~~~~~~~~~~~~~~~~~~~~
The following list details the parameters used across all meshing procedures. If the parameters are not optional, they can never be safely omitted.

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - filetype
      - The file type of the GMESH output. By default, the output will be a .msh file.
      - Yes.
    * - type
      - Box by default. Can be set to cylinder for a cylindrical outer domain (or a circle in 2D).
      - Yes.
    * - x_range
      - The minimum and maximum x-coordinates that define the meshing domain.
      - No.
    * - y_range
      - The minimum and maximum y-coordinates that define the meshing domain.
      - No.
    * - dimension
      - The dimension of the mesh. 2 for 2D meshes, 3 for 3D meshes.
      - No.
    * - global_scale
      - A scaling factor applied uniformly to all length scales. 1 by default, to create no effect.
      - Yes.
    * - threshold_upstream_distance
      - The extent to which the turbine meshes should extend in the negative wake direction.
      - No.
    * - threshold_downstream_distance
      - The extent to which the turbine meshes should extend in the positive wake direction.
      - No.
    * - threshold_rotor_distance
      - The radius formed by a rotation of the turbine's rotor.
      - No.
    * - background_length_scale
      - The length scale of the meshing across the entire domain.
      - No.
    * - length_scale (turbine)
      - The length scale of the meshing near the turbines.
      - No.
    * - num_turbines
      - The expected number of turbines.
      - No.
    * - inflow_angle
      - If there is one incoming wind vector, this flag will orient the turbines to face the wind. Input the angle in degrees.
      - Yes.

Creating Turbines
~~~~~~~~~~~~~~~~~~~~~~

As discussed above, the number of turbines must be specified using the appropriate field.
For each turbine, create a numerical field and populate it with the three relevant turbine parameters.
Examples are shown :ref:`here <yaml_params>`, under the turbine header. The table below describes the three
turbine parameters.

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - x
      - The x-coordinate of the center of the turbine.
      - No.
    * - y
      - The y-coordinate of the center of the turbine.
      - No.
    * - HH
      - The turbine's hub height. This parameter determines how high off the terrain the turbine center is placed. Set to 100 by default.
      - Yes.

Creating Farm Refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AeroMesh defines a "farm" as, at a minimum, a bounding geometry that surrounds all the turbines and an associated length scale for this region.
Note that farms are optional and do not need to be included in the YAML if not desired. Examples are shown :ref:`here <yaml_params>`, under the farm header. The table below describes the farm parameters.


.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - type
      - None (disabled) by default. Can be set to cylinder or box if desired.
      - Yes.
    * - threshold_distance
      - By default, a farm region is defined by a minimum bounding rectangle surrounding all the turbines. This parameter extends the bounding region by its value.
      - Yes.
    * - length_scale
      - The length scale of the points contained within the farms.
      - Yes.

Creating Custom Refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating custom refinements is almost identical to creating turbines. The number of custom refinements must be specified. Depending on the type of
refinement (box or cylinder) desired, different parameters must be supplied. In 2D meshes, the 3D refinement types will be replaced with their
2D analogs without the need for any flags from the user. Examples are shown :ref:`here <yaml_params>`, under the refine_custom header. The table below describes the relevant
sub-parameters. 

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - type
      - Either "box" for a box refine or "cylinder" for a cylindrical one.
      - No.
    * - x_range
      - In box refines, the x range of the refinement. In cylinders, the x-coordinate of the center.
      - No.
    * - y_range
      - In box refines, the y range of the refinement. In cylinders, the y-coordinate of the center.
      - No.
    * - z_range
      - The z range of the refinement.
      - No.
    * - radius
      - The radius of a cylindrical refinement. Unused by box meshes and does not need to be included.
      - No.
    * - length_scale
      - The length scale across the custom refinement.
      - No.


3D Parameters Only
~~~~~~~~~~~~~~~~~~~~~~
The following list details the parameters used only in 3D simulations. If the simulation is in 2D,
these parameters can be safely omitted.

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - terrain_path
      - The path to a file containing valid terrain data, where the terrain is a function f(x, y) = z. If omitted, the domain will have a smooth bottom face.
      - Yes.
    * - height
      - The extension of the wind farm in the z-direction.
      - No.
    * - aspect_ratio
      - The ratio of nodes in the z-direction to nodes in the x-y plane. Used to create anisotropic effects, if desired.
      - Yes.
    * - upper_aspect_ratio
      - Similar to aspect_ratio, but applied above the threshold distance instead to create expanded nodes in the z-direction. Note that the resultant size after the application of both aspect ratios must be at least the original height of the domain to prevent undefined behavior.
      - Yes.
    * - aspect_distance
      - The z-distance up to which the anisotropic effects generated by aspect_ratio will extend.
      - Yes.


2D Parameters Only
~~~~~~~~~~~~~~~~~~~~~~

In 2D meshes, the type of refinement that defines the turbine may be specified. This customizability is done using the
refine[turbine][type] flag. The table below describes the flag.

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
      - Optional
    * - type
      - Determines whether turbines are meshed using a rectangular wake or as a large circle. The options are rectangle and circle respectively. Rectangles are used by default.
      - Yes.

Output Format
~~~~~~~~~~~~~~~~~~~~~~

AeroMesh produces an output file named out.extension that contains the mesh data. The extension type can be controlled by the
"filetype" field. By default, AeroMesh produces msh files. However, it supports any other data formats 
`handled by GMSH <https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-command_002dline-interface>`_ and the additional xdmf output type if desired.