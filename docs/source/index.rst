Welcome to QMAP's documentation!
================================

QMAP is a tool for quantum circuit compilation developed as part of the `Munich Quantum Toolkit <https://mqt.readthedocs.io>`_ (*MQT*) by the `Chair for Design Automation <https://www.cda.cit.tum.de/>`_ at the `Technical University of Munich <https://www.tum.de>`_. It builds upon `MQT Core <https://github.com/cda-tum/mqt-core>`_, which forms the backbone of the MQT.

We recommend you to start with the :doc:`installation instructions <Installation>` or by reading our overview paper :cite:labelpar:`wille2023qmap`.
Then proceed to the :doc:`mapping page <Mapping>` or the :doc:`synthesis/optimization page <Synthesis>` and read the :doc:`reference documentation <library/Library>`.
If you are interested in the theory behind QMAP, have a look at the publications in the :doc:`publication list <Publications>`.

We appreciate any feedback and contributions to the project. If you want to contribute, you can find more information in the :doc:`Contribution <Contributing>` guide. If you are having trouble with the installation or the usage of QMAP, please let us know at our :doc:`Support <Support>` page or by reaching out to us at `quantum.cda@xcit.tum.de <mailto:quantum.cda@xcit.tum.de>`_.

----

.. toctree::
    :hidden:

    self

.. toctree::
    :maxdepth: 2
    :caption: User Guide
    :glob:

    Installation
    Mapping
    Synthesis
    Publications

.. toctree::
    :maxdepth: 2
    :caption: Developers
    :glob:

    Contributing
    DevelopmentGuide
    Support

.. toctree::
    :maxdepth: 6
    :caption: API Reference
    :glob:

    library/Library
