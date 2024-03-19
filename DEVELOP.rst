Hidive development guide
""""""""""""""""""""""""""""""

We welcome contributions to the codebase. This document provides guidance
on setting up your development environment to test new features locally, then
contributing your changes to the main repository.


About the codebase
------------------
Hidive is primarily written in Rust for fast operations: multi-threaded data
retrieval from Google Cloud Storage, genome assembly, columnar storage, etc.).
Optionally, a Python library binding to the Rust codebase is available for
interactive use of the Rust code, including visualization and notebook-based
exploration.

Developers need not know both languages to meaningfully contribute to Hidive
development. We welcome improvements and additional features to both the
Rust and Python parts of the codebase.

Setting up your development environment
---------------------------------------


We recommend using Visual Studio Code (VSCode) as your development environment.
Follow the steps below to set up your development environment:

1. Install Rust
   Open a terminal and run the following command to install Rust:

   .. code-block:: bash

       curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh


   This will download a script and start the installation. Once the installation
   is complete, close the terminal and open a new one to ensure the changes take effect.
   You can verify that you have installed Rust correctly by running the following
   command in your terminal:

   .. code-block:: bash

       rustc --version

   If Rust is installed correctly, this command should display the version of Rust currently installed.

2. Install Visual Studio Code
   You can download it from `here <https://code.visualstudio.com/download>`_.
   Follow the instructions on the website to install it.

3. Install Rust and Python extensions for VSCode
   Open VSCode, go to the Extensions view by clicking on the Extensions icon in the Activity Bar on the side of the window. In the Extensions view, search for and install the following extensions:
   
   - 'rust-analyzer': This extension provides advanced language support for Rust.
   - 'CodeLLDB': This extension is an optional but recommended tool for developing and debugging Rust code.
   - 'Jupyter': This extension provides Jupyter notebook support for VSCode.

4. (Optional: development in the Python codebase.) Create a virtual environment for development
   Python's built-in `venv` module allows you to create virtual environments. These environments have their own installation directories and don't share libraries with other virtual environments. Here's how to create one:

   - Open a terminal and navigate to the directory where you want to create the virtual environment.
   - Run the following command to create a new virtual environment:

   .. code-block:: bash

       python3 -m venv hidive-env

   This will create a new directory called `hidive-env` in your current directory. This directory will contain the Python executable files and a copy of the pip library which you can use to install other packages.

   - To activate the virtual environment, run the following command:

   .. code-block:: bash

       source hidive-env/bin/activate

   Once the virtual environment is activated, the name of your virtual environment will appear on left side of the terminal prompt. This will let you know that the virtual environment is currently active. 

   - In the virtual environment, you can use the command `pip` to install packages that will be isolated from the global Python installation. Install the required packages for Hidive development by running the following command:

   .. code-block:: bash

       pip install -r dev-requirements.txt

   - When you are done with your work, you can deactivate the virtual environment by running the following command:

   .. code-block:: bash

       deactivate

   This command will deactivate the virtual environment and you will return to your normal shell.

5. Compile and install Hidive into the virtual environment
   To compile and install Hidive, we will use the `maturin develop --release` command. `maturin` is a build system for Rust-based Python extensions, and the `develop` command compiles and installs the package into the current Python interpreter. The `--release` flag is used to compile the package in release mode, which includes optimizations.

   First, change into the `src/pydive` directory.

   .. code-block:: bash

       cd src/pydive

   Run the following command in your terminal:

   .. code-block:: bash

       maturin develop --release

   This command will compile the Hidive Rust code and install the resulting Python package into your active virtual environment. This means you can now import and use the Hidive library in your Python scripts and Jupyter notebooks.

6. Open an interactive Python session
   At the python terminal, import the Hidive library by running the following code:

   .. code-block:: python

       (hidive-env) user@host pydive % python
       Python 3.12.2 (main, Feb  6 2024, 20:19:44) [Clang 15.0.0 (clang-1500.1.0.2.5)] on darwin
       Type "help", "copyright", "credits" or "license" for more information.
       >>> import pydive

   If the library imports successfully, you are ready to start using Hidive in your notebook.

Each time you make changes to the codebase, recompile the library by rerunning
step 5, then trying out the changes in step 6.

Now, you are ready to start developing with Pydive!


Contributing to Hidive
----------------------------

1. Fork the Hidive repository
   Go to the `Hidive repository <https://github.com/broadinstitute/hidive>`_ and click on the "Fork" button. This will create a copy of the repository in your own GitHub account.

2. Clone the forked repository
   On your GitHub account, navigate to the forked repository and click on the "Clone or download" button. Copy the URL.
   Open a terminal and run the following git command:
   
   .. code-block:: bash

       git clone "url you just copied"

3. Create a new branch
   Change to the repository directory on your computer (if you are not already there):

   .. code-block:: bash

       cd hidive

   Now create a branch using the `git checkout` command:

   .. code-block:: bash

       git checkout -b your-new-branch-name

4. Make necessary changes and commit those changes
   Now you can go ahead and make changes to the files. Once you've made changes or added files, you can see them listed with `git status`. Add those changes with `git add` and then commit those changes:

   .. code-block:: bash

       git add .
       git commit -m "commit message"

5. Push changes to GitHub
   Push your changes using the command `git push`:

   .. code-block:: bash

       git push origin your-new-branch-name

6. Submit your changes for review
   If you go to your repository on GitHub, you'll see a "Compare & pull request" button. Click on that button and describe the changes you made. Once you submit the pull request, a Hidive reviewer will review your changes.


Code Style Guidelines
---------------------

We follow the official style guides for our code. For Rust, we adhere to the `Rust Style Guide <https://rust-lang.github.io/rfcs/1607-style-guide.html>`_. For Python, we follow the `PEP 8 Style Guide <https://pep8.org/>`_. Please ensure your contributions adhere to these standards.

Testing
-------

We use pytest for our Python tests and cargo test for our Rust tests. Please add tests for new features and ensure all tests pass before submitting a pull request.

Documentation
-------------

Please update the documentation to reflect any changes you make to the codebase. This includes comments in the code, docstrings, and our user guides.

Issue Tracking
--------------

We use GitHub issues to track work on Hidive. If you're adding a new feature or fixing a bug, please create an issue describing the work.

Communication
-------------

If you have any questions or want to discuss your work, please join our community chat on Slack or by email. Our team is always happy to help.

Code Review Process
-------------------

After you submit your pull request, it will be reviewed by at least one core contributor. We'll provide feedback and may request changes. Once your pull request is approved, it will be merged into the main codebase and automatically released as an incremental version update on Pypi.