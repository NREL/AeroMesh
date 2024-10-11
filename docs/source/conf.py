import pathlib
import sys
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

project = 'AeroMesh'
copyright = '2024, National Renewable Energy Laboratory'
author = 'National Renewable Energy Laboratory'
release = 'v0.1'

extensions = ['sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = []


html_theme = 'furo'
html_static_path = ['_static']
