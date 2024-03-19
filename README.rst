hidive
""""""""""""

|GitHub release| |PyPI version hidive|

.. |GitHub release| image:: https://img.shields.io/github/release/broadinstitute/hidive.svg
   :target: https://github.com/broadinstitute/hidive/releases/

.. |PyPI version hidive| image:: https://img.shields.io/pypi/v/hidive.svg
   :target: https://pypi.python.org/pypi/hidive/

Hidive is a targeted genome co-assembler for biobank-scale long-read and short-read data.

Documentation for the API can be found on the `documentation page <https://broadinstitute.github.io/hidive/>`_.


Prerequisites
-------------

Hidive is designed to access local files or data in Google Cloud Storage (GCS). Within certain cloud-computing environments (i.e. Terra, All of Us Researcher Workbench), access to GCS is already configured. For accessing files in GCS on your local machine, you will also need to install the `Google Cloud CLI <https://cloud.google.com/sdk/docs/install-sdk>`_. Then, configure your `Application Default Credentials (ADC) <https://cloud.google.com/docs/authentication/provide-credentials-adc#local-dev>`_.


Installation
------------

At the moment, Hidive is under active development and is not yet available through standard bioinformatic software channels (pip, conda, cargo, etc.).


Building from source
--------------------

To build the primary Rust source, follow the procedure below.

.. code-block:: bash
   # Install Rust.
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

   # Clone repository.
   git clone https://github.com/broadinstitute/hidive.git
   cd hidive

   # Build Hidive's Rust codebase.
   cargo build --release

To optionally build the Python bindings to the Rust codebase, follow the procedure below after performing the above steps.

.. code-block:: bash
   # Create a Python virtual environment and install Maturin, the tool that
   # will compile the Rust and Python code into a complete library.
   # For more information on Maturin, visit https://github.com/PyO3/maturin .
   python -mvenv venv
   . venv/bin/activate
   pip install maturin

   # Build the library (with release optimizations) and install it in
   # the currently active virtual environment.
   cd src/pydive/
   maturin develop --release

Supported platforms
-------------------

Hidive is compiled for Linux and MacOSX. Windows is not currently supported.

Getting help
------------

If you encounter bugs or have questions/comments/concerns, please file an issue on our `Github page <https://github.com/broadinstitute/hidive/issues>`_.

Developers' guide
-----------------

For information on contributing to Hidive development, visit our `developer documentation <DEVELOP.rst>`_.
