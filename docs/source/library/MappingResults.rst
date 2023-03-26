Results
=======

This class captures additional results from the :code:`compile` method.

    .. currentmodule:: mqt.qmap
    .. autoclass:: MappingResults

In addition to the mapped circuit

    .. autoattribute:: MappingResults.mapped_circuit

it also includes runtime information of the mapping procedure

    .. autoattribute:: MappingResults.time
    .. autoattribute:: MappingResults.timeout

and information about the input and output circuits.

    .. autoattribute:: MappingResults.input
    .. autoattribute:: MappingResults.output

If specified in the call to the :code:`compile` method, the `weighted MaxSAT formula <http://www.maxhs.org/docs/wdimacs.html>`_ is also tracked.

    .. autoattribute:: MappingResults.wcnf

If heuristic mapping is used and :code:`debug` is set to :code:`true`, some benchmark information is provided for the whole circuit

    .. autoattribute:: MappingResults.heuristic_benchmark

and for each layer separately.

    .. autoattribute:: MappingResults.layer_heuristic_benchmark

In addition, the class provides methods to export to other formats.

    .. automethod:: MappingResults.json
    .. automethod:: MappingResults.csv
