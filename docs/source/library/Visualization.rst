Visualization
===========================

All visualization functionality is bundled in the sub-module :code:`mqt.qmap.visualization`.

Search graph visualization
############################

The recommended way to visualize search graphs of a mapping process is via a :code:`SearchVisualizer` object, which can be passed to the :code:`compile` function when mapping a circuit to enable data logging during mapping, after which this data can be visualized by calling the method :code:`visualize_search_graph` of the :code:`SearchVisualizer` object.
Closing the :code:`SearchVisualizer` object will cleanup the data logged during mapping, which will however also be done automatically when the python process terminates (if not custom data logging path is used).

    .. note::
        :code:`visualize_search_graph` returns an IPython display object and therefore requires to be executed in a jupyter notebook or jupyter lab environment.

    .. currentmodule:: mqt.qmap.visualization
    .. autoclass:: SearchVisualizer
