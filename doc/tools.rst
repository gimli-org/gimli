.. _sec:tools:

Tools
=====

Mesh manipulation
-----------------

meshconvert
...........
.. program-output:: meshconvert -h
    :nostderr: 

dctriangle
..........

Polytools
.........

meshmerge
.........
.. program-output:: meshmerge -h
    :nostderr: 

Data manipulation
-----------------

datamerge
.........

.. program-output:: datamerge -h
    :nostderr: 

dataedit
........


Visualisation
-------------

pytripatch
..........

Generic 2 dimensional mesh/model viewer:

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
