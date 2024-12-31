Developer API
====================================

Meshing in **AeroMesh** is achieved through calls to several under-the-hood functions. For practical use, these functions should never be
referenced by another process. However, their functionality is detailed here so that developers who want to fork the repository
have a resource to understand **AeroMesh**'s backend.

3D Meshing functions
-----------------------------------

These functions, found in ``src/functions3D.py``, serve as the backend for the 3D meshing protocols.

.. autofunction:: src.functions3D.generateTurbines
    
.. autofunction:: src.functions3D.placeTurbine

.. autofunction:: src.functions3D.anisotropyScale

.. autofunction:: src.functions3D.calcEllipse

.. autofunction:: src.functions3D.refineFarm3D

2D Meshing functions
-----------------------------------

These functions, found in ``src/functions2D.py``, serve as the backend for the 2D meshing protocols.

.. autofunction:: src.functions2D.buildFarms2D

.. autofunction:: src.functions2D.refineFarm2D

Terrain Meshing functions
-----------------------------------

These functions, found in ``src/terrain.py``, generate the terrain for the various meshing programs.

.. autofunction:: src.terrain.buildTerrainDefault

.. autofunction:: src.terrain.buildTerrainFromFile

.. autofunction:: src.terrain.buildTerrain2D