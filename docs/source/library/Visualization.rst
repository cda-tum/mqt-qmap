Visualization
=============

.. currentmodule:: mqt.qmap

All visualization functionality is bundled in the sub-module :mod:`mqt.qmap.visualization`.

Search graph visualization
##########################

.. currentmodule:: mqt.qmap.visualization

The recommended way to visualize search graphs of a mapping process is via a :class:`SearchVisualizer` object, which can be passed to the :func:`compile` function to enable data logging during mapping, after which this data can be visualized by calling the method :meth:`SearchVisualizer.visualize_search_graph`.

    .. note::
        :meth:`visualize_search_graph` returns an IPython display object and therefore requires to be executed in a jupyter notebook or jupyter lab environment.

    .. autoclass:: SearchVisualizer
        :special-members: __init__
        :members:
