Getting Started
====================================

Installation
------------

To use AeroMesh, clone `the GitHub repository <https://github.nrel.gov/agushin/AeroMesh>`_ onto your local machine.
Change directory into AeroMesh and run the command:

.. code-block:: console

    $ conda env create -f requirements/environment.yml

Activate the environment with:

.. code-block:: console

    $ conda activate AeroMesh

Usage
------------

To start, create a YAML file that details the specifications of your wind farm. A sample structure for the YAML
file can be viewed :ref:`here <yaml_params>`. Once this file is created,
the program can be executed via:

.. code-block:: console

    $ python aeromesh.py file.yaml