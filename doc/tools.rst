.. _sec:tools:

Tools
=====

Mesh manipulation
-----------------

meshConvert
...........
.. program-output:: convertMesh -h
    :nostderr:

meshMerge
.........
.. program-output:: meshMerge -h
    :nostderr:

createParameterMesh
...................
.. program-output:: createParameterMesh -h
    :nostderr:

pyCreateSurface
...............
.. program-output:: pyCreateSurface -h
    :nostderr:

recountPara
...........
.. program-output:: recountPara -h
    :nostderr:

Data manipulation
-----------------

datamerge
.........

.. program-output:: dataMerge -h
    :nostderr:

Visualisation
-------------

pytripatch
..........

Generic two-dimensional mesh/model viewer:

.. program-output:: pytripatch -h
    :nostderr:

.. admonition:: Examples

    * Show a 2d mesh only

    ::

        pytripatch mesh.bms

    * Color the mesh based on data. Data can be a file with a list of values i.e., one data for each cell.


    ::

        pytripatch -d data.txt mesh.bms

    * Color the mesh based on data. Data can be a string of a data vector shipped with the :gimliapi:`GIMLI::Mesh` see :gimliapi:`GIMLI::Mesh::addExportData` and binary format v2. Or internal attributes like ``attribute`` or ``marker``


    ::

        pytripatch --data=marker mesh.bms

.. Move to pybert?
.. showSens
.. ........
.. .. program-output:: showSens -h
..     :nostderr:
