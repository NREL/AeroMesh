import pathlib
import sys
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

project = 'AeroMesh'
copyright = '2025, Alliance for Sustainable Energy, LLC'
author = 'Alliance for Sustainable Energy, LLC'
release = 'v0.1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme'
    ]

templates_path = ['_templates']
exclude_patterns = []


html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
