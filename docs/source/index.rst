Welcome to PandaDock Documentation
==================================

.. image:: _static/logo.png
   :width: 300px
   :align: center
   :alt: PandaDock Logo

PandaDock is a modular, multi-strategy, high-performance molecular docking software designed for drug discovery and computational chemistry research. It provides a comprehensive suite of docking algorithms, scoring functions, and analysis tools to predict protein-ligand binding poses and affinities.

.. note::
   PandaDock v3.0.0 introduces machine learning-enhanced scoring, improved performance, and a redesigned modular architecture.

Quick Start
-----------

Install PandaDock using pip:

.. code-block:: bash

   pip install pandadock

Basic usage:

.. code-block:: python

   from pandadock import PandaDock
   
   # Initialize docking engine with PandaML algorithm
   docker = PandaDock(engine='ml', scoring='pandaml')
   
   # Perform docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       output_dir='results/'
   )

Key Features
------------

🔬 **Novel PandaDock Algorithms**
   - **PandaCore**: Robust baseline algorithm with reliable performance
   - **PandaML**: Advanced machine learning algorithm with superior affinity prediction
   - **PandaPhysics**: Physics-based algorithm specialized for metal coordination
   - Flexible receptor docking with side-chain sampling
   - **Advanced metal docking for metalloproteins**

⚡ **High Performance**
   - GPU acceleration support (CUDA)
   - Parallel processing capabilities
   - Optimized algorithms for large-scale screening

🧠 **AI-Powered Scoring**
   - Machine learning rescoring with transformer models
   - Energy decomposition analysis
   - Binding affinity prediction with uncertainty quantification

📊 **Comprehensive Analysis**
   - Interactive HTML reports
   - Detailed interaction analysis
   - IC50 and binding kinetics prediction
   - Drug-likeness assessment

🔧 **Modular Architecture**
   - Plugin-based scoring functions
   - Customizable energy terms
   - Extensible with custom algorithms

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/installation
   user_guide/quickstart
   user_guide/configuration
   user_guide/docking_modes
   user_guide/metal_docking
   user_guide/scoring_functions
   user_guide/analysis

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   changelog
   license

Performance Benchmarks
----------------------

PandaDock has been benchmarked against established docking software:

.. list-table:: PandaDock Algorithm Performance (285 PDBbind complexes)
   :header-rows: 1
   :widths: 25 20 20 20 15

   * - Algorithm
     - Affinity R²
     - Success Rate (%)
     - Mean RMSD (Å)
     - Speed (s)
   * - **PandaML**
     - **0.878**
     - **46.7**
     - **3.22**
     - **26.1**
   * - **PandaPhysics**
     - 0.671
     - **48.8**
     - 2.18
     - 42.9
   * - **PandaCore**
     - 0.738
     - 46.3
     - 3.63
     - 32.9

*Success Rate = RMSD < 2Å. Benchmarked on diverse protein-ligand complexes.*

**Metal vs Non-Metal Specialization:**

- **93 metal complexes** (32.6%) vs **192 non-metal complexes** (67.4%) analyzed
- **PandaPhysics** excels at metal coordination chemistry  
- **PandaML** maintains consistent performance across both complex types
- Metal complexes are 32% more challenging computationally

Citation
--------

If you use PandaDock in your research, please cite:

.. code-block:: bibtex

   @article{pandadock2025,
     title={PandaDock: A Machine Learning-Enhanced Molecular Docking Platform},
     author={Panda, Pritam Kumar},
     journal={Journal of Computational Chemistry},
     year={2025},
     volume={46},
     pages={1--15},
     doi={10.1002/jcc.xxxxx}
   }

Support and Community
--------------------

- **GitHub**: `github.com/pritampanda15/pandadock <https://github.com/pritampanda15/pandadock>`_
- **Issues**: Report bugs and request features on `GitHub Issues <https://github.com/pritampanda15/pandadock/issues>`_
- **Discussions**: Join the community on `GitHub Discussions <https://github.com/pritampanda15/pandadock/discussions>`_
- **Email**: pritam@stanford.edu

License
-------

PandaDock is released under the MIT License. See the `LICENSE <https://github.com/pritampanda15/pandadock/blob/main/LICENSE>`_ file for details.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`