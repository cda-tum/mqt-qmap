Visualization
=============

.. currentmodule:: mqt.qmap

All visualization functionality is bundled in the sub-module :mod:`mqt.qmap.visualization`.

Search graph visualization
##########################

.. currentmodule:: mqt.qmap.visualization

The recommended way to visualize search graphs of a mapping process is via a :class:`SearchVisualizer` object, which can be passed to the :func:`mqt.qmap.compile` function to enable data logging during mapping, after which this data can be visualized by calling the method :meth:`SearchVisualizer.visualize_search_graph`.

.. note::
    :meth:`SearchVisualizer.visualize_search_graph` returns an IPython display object and therefore requires to be executed in a jupyter notebook or jupyter lab environment.

.. note::
    Automatic layouting of architecture or search nodes requires `Graphviz <https://graphviz.org/>`_ to be installed (except for the layouting method :code:`walker`). If Graphviz is called without it being installed, it will ensue in an error such as:

    :code:`FileNotFoundError: [Errno 2] "sfdp" not found in path.`

    Consequently, the only way to use :meth:`SearchVisualizer.visualize_search_graph` without Graphviz is by passing explicit architecture node positions or hiding the architecture graph by passing :code:`show_layout=None`.

.. autoclass:: SearchVisualizer
    :special-members: __init__
    :members:
