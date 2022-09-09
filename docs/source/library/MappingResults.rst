Results
=======

This class captures additional results from the :code:`compile` method.

    .. currentmodule:: mqt.qmap
    .. autoclass:: MappingResults

In addition to the mapped circuit

    .. autoattribute:: MappingResults.mapped_circuit

it also includes runtime information of the mapping procedure.

    .. autoattribute:: MappingResults.time
    .. autoattribute:: MappingResults.timeout

Information about the input and output circuits are also stored.

    .. autoattribute:: MappingResults.input
    .. autoattribute:: MappingResults.output

If specified in the call to the :code:`compile` method, the `weighted max-SAT formula <http://www.maxhs.org/docs/wdimacs.html>`_ is also tracked.

    .. autoattribute:: MappingResults.wcnf

In addition, the class provides methods to export to other formats.

    .. automethod:: MappingResults.json
    .. automethod:: MappingResults.csv
