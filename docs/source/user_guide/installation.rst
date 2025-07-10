Installation Guide
=================

PandaDock can be installed through multiple methods depending on your needs and system configuration.

Prerequisites
-------------

**System Requirements:**
- Python 3.8 or higher
- 4 GB RAM minimum (16 GB recommended)
- 1 GB free disk space
- CUDA-compatible GPU (optional, for acceleration)

**Operating Systems:**
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS (10.14+)
- Windows 10/11 (via WSL recommended)

Basic Installation
------------------

Install via pip (recommended):

.. code-block:: bash

   pip install pandadock

Install with all dependencies:

.. code-block:: bash

   pip install pandadock[all]

Development Installation
------------------------

For development or to access the latest features:

.. code-block:: bash

   git clone https://github.com/pritampanda15/pandadock.git
   cd pandadock
   pip install -e .

Optional Dependencies
---------------------

PandaDock supports several optional dependency groups:

**Machine Learning Features:**

.. code-block:: bash

   pip install pandadock[ml]

Includes:
- PyTorch for neural network models
- Transformers for advanced ML scoring
- XGBoost for ensemble methods

**Visualization:**

.. code-block:: bash

   pip install pandadock[viz]

Includes:
- PyMOL for 3D structure visualization
- RDKit for chemical structure handling
- py3Dmol for interactive visualization

**Chemical Processing:**

.. code-block:: bash

   pip install pandadock[chem]

Includes:
- RDKit for cheminformatics
- BioPython for structure processing
- OpenMM for molecular mechanics

**Web Interface:**

.. code-block:: bash

   pip install pandadock[web]

Includes:
- Flask for web interface
- FastAPI for REST API
- Uvicorn for async server

**GPU Acceleration:**

.. code-block:: bash

   pip install pandadock[gpu]

Includes:
- CuPy for GPU-accelerated computations
- CUDA toolkit integration

**Development Tools:**

.. code-block:: bash

   pip install pandadock[dev]

Includes:
- Testing frameworks (pytest)
- Code formatting (black, flake8)
- Documentation tools (Sphinx)

GPU Setup (Optional)
--------------------

For GPU acceleration, install CUDA toolkit:

**Linux/WSL:**

.. code-block:: bash

   # Install CUDA 11.8 (recommended)
   wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
   sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
   wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda-repo-ubuntu2004-11-8-local_11.8.0-520.61.05-1_amd64.deb
   sudo dpkg -i cuda-repo-ubuntu2004-11-8-local_11.8.0-520.61.05-1_amd64.deb
   sudo cp /var/cuda-repo-ubuntu2004-11-8-local/cuda-*-keyring.gpg /usr/share/keyrings/
   sudo apt-get update
   sudo apt-get -y install cuda

**macOS:**
GPU acceleration is not available on macOS due to CUDA limitations. CPU-only installation is recommended.

Conda Installation
------------------

PandaDock can also be installed via conda:

.. code-block:: bash

   conda install -c conda-forge pandadock

Create a dedicated environment:

.. code-block:: bash

   conda create -n pandadock python=3.9
   conda activate pandadock
   conda install -c conda-forge pandadock

Docker Installation
-------------------

Run PandaDock in a container:

.. code-block:: bash

   docker pull pritampanda15/pandadock:latest
   docker run -it --gpus all -v $(pwd):/workspace pritampanda15/pandadock:latest

For CPU-only:

.. code-block:: bash

   docker run -it -v $(pwd):/workspace pritampanda15/pandadock:latest

Verification
------------

Test your installation:

.. code-block:: python

   import pandadock
   print(f"PandaDock version: {pandadock.__version__}")
   
   # Test basic functionality
   from pandadock import PandaDock
   docker = PandaDock()
   print("Installation successful!")

Run the test suite:

.. code-block:: bash

   python -m pytest tests/

Troubleshooting
---------------

**Common Issues:**

1. **Import Error:**
   
   .. code-block:: text
   
      ImportError: No module named 'pandadock'
   
   **Solution:** Ensure pip installation completed successfully and you're using the correct Python environment.

2. **CUDA Not Found:**
   
   .. code-block:: text
   
      RuntimeError: CUDA not available
   
   **Solution:** Install CUDA toolkit or use CPU-only mode by setting ``device='cpu'``.

3. **Memory Error:**
   
   .. code-block:: text
   
      MemoryError: Unable to allocate array
   
   **Solution:** Reduce batch size or use a machine with more RAM.

4. **Permission Error:**
   
   .. code-block:: text
   
      PermissionError: [Errno 13] Permission denied
   
   **Solution:** Use ``pip install --user pandadock`` or run with appropriate permissions.

**Getting Help:**

- Check the `FAQ <https://pandadock.readthedocs.io/en/latest/faq.html>`_
- Search `GitHub Issues <https://github.com/pritampanda15/pandadock/issues>`_
- Join `GitHub Discussions <https://github.com/pritampanda15/pandadock/discussions>`_

Environment Variables
---------------------

Configure PandaDock behavior with environment variables:

.. code-block:: bash

   # Set default GPU device
   export PANDADOCK_DEVICE=cuda:0
   
   # Set temporary directory
   export PANDADOCK_TMPDIR=/tmp/pandadock
   
   # Enable debug logging
   export PANDADOCK_DEBUG=1
   
   # Set number of CPU threads
   export PANDADOCK_THREADS=8

Next Steps
----------

After installation, check out the :doc:`quickstart` guide to begin using PandaDock.