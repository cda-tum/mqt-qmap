Publications
============

*QMAP* is academic software. Thus, many of its built-in algorithms have been published as scientific papers.
See :cite:labelpar:`wille2023qmap` for a general overview of *QMAP* and its features.
If you want to cite this article, please use the following BibTeX entry:

    .. code-block:: bibtex

        @inproceedings{qmap,
            title = {{QMAP: A Quantum Circuit Mapping Tool}},
            booktitle = {International Symp. on Physical Design},
            author = {Wille, Robert and Burgholzer, Lukas},
            year = {2023}
        }

*QMAP* is part of the Munich Quantum Toolkit, which is described in :cite:labelpar:`willeMQTHandbookSummary2024`.
If you want cite the Munich Quantum Toolkit, please use the following BibTeX entry:

    .. code-block:: bibtex

        @inproceedings{mqt,
            title = {The {{MQT}} Handbook: {{A}} Summary of Design Automation Tools and Software for Quantum Computing},
            shorttitle = {{The MQT Handbook}},
            booktitle = {IEEE International Conference on Quantum Software (QSW)},
            author = {Wille, Robert and Berent, Lucas and Forster, Tobias and Kunasaikaran, Jagatheesan and Mato, Kevin and Peham, Tom and Quetschlich, Nils and Rovara, Damian and Sander, Aaron and Schmid, Ludwig and Schoenberger, Daniel and Stade, Yannick and Burgholzer, Lukas},
            date = {2024},
            doi = {10.1109/QSW62656.2024.00013},
            eprint  = {2405.17543},
            eprinttype = {arxiv},
            addendum = {A live version of this document is available at \url{https://mqt.readthedocs.io}}
        }

If you use *QMAP* in your work, we would appreciate if you cited

- :cite:labelpar:`zulehnerEfficientMethodologyMapping2019` when using the heuristic mapper,
- :cite:labelpar:`willeMappingQuantumCircuits2019` when using the exact mapper,
- :cite:labelpar:`peham2023DepthOptimalSynthesis` when using the Clifford circuit synthesis approach,
- :cite:labelpar:`schmid2024HybridCircuitMapping` when using the hybrid mapper for neutral atom quantum computers,
- :cite:labelpar:`stadeAbstractModelEfficient2024` when using the neutral atom logical array compiler (NALAC), and
- :cite:labelpar:`stadeOptimalStatePreparation2024` when using the optimal state preparation for neutral atoms (NASP).

Furthermore, if you use any of the particular algorithms such as

- the heuristic mapping scheme using teleportation :cite:labelpar:`hillmichExlpoitingQuantumTeleportation2021`
- the search space limitation techniques of the exact mapper (some of which are enabled per default) :cite:labelpar:`burgholzer2022limitingSearchSpace`
- the method for finding (near-)optimal subarchitectures :cite:labelpar:`peham2023OptimalSubarchitectures`

please consider citing their respective papers as well. A full list of related papers is given below.

.. bibliography::
