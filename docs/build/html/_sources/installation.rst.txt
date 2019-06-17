.. _installation:

Installation
============

This pipeline was build upon many dependencies and `Python <http://www.python.org/>`_ code. Therefore there was the need to use a package/environment manager. For this task we highly recommend `CONDA <https://docs.conda.io/en/latest/>`_ . The following installation instructions assume the user wants to install and use `CONDA <https://docs.conda.io/en/latest/>`_. If you already have any recent conda distribution on your system do not use the `Quick install Guide`_, as it is it will install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. If you do not which to use `CONDA <https://docs.conda.io/en/latest/>`_ please look at all the dependencies needed to use SNAPPy in the `'environment.yaml' file <https://github.com/PMMAraujo/snappy/blob/master/environment.yaml>`_.

 .. _quick:

Quick install Guide
^^^^^^^^^^^^^^^^^^^

1) Download SNAPPy from `here <https://github.com/PMMAraujo/snappy>`_ and put it in a desired location

2) Open a terminal and move to the downloaded folder::

    cd snappy

3) Give permission to the snappy installation script to run, and run it::

    chmod +x snappy_installer_x86_64.sh 
    ./snappy_installer_x86_64.sh

This task is time consuming and may take several minutes.

4) After the SNAPPy's dependencies installation a series of tests will be executed using `pytest <https://docs.pytest.org/en/latest/>`_. Please ensure no errors occur.
    a) If there was any errors clean the downloaded folder and start from step 1. If the problem persist please try to follow the `Install point-by-point`_. If you are still unable to install SNAPPy after this contact us to EMAIL.

5) SNAPPy is now ready to use! When you want to use SNAPPy open a new terminal inside SNAPPy's folder and activate the conda environment by typing::

    conda activate snappy
 

For a quick start on SNAPPy's usage look at the :ref:`tutorials` section.

 .. _point-by-point:

Install point-by-point
^^^^^^^^^^^^^^^^^^^^^^

1) Open a terminal, download Miniconda3 version 4.6.14 and install it::

    wget https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ./miniconda.sh
    bash ./miniconda.sh -b -p $HOME/miniconda

This task may take several minutes.

2) Download SNAPPy from `here <https://github.com/PMMAraujo/snappy>`_ and put it in a desired location

3) Move to the downloaded folder::

    cd snappy

4) Activate conda and add it tho the bashrc::

    export PATH="$HOME/miniconda/bin:$PATH"
    echo -e '\n# added by SNAPPy to run miniconda' >> ~/.bashrc
    MISSING

Note: It is not actually needed to add conda path to the bashrc but it makes it easier to use.

5) Create SNAPPy's conda environment::

    conda-env create -f environment.yaml

6) Run the tests to ensure the installation was successful. Please ensure no errors occur::

    py.test

7) SNAPPy is now ready to use! When you want to use SNAPPy open a new terminal inside SNAPPy's folder and activate the conda by typing::

    conda activate snappy

For a quick start on SNAPPy's usage look at the :ref:`tutorials` section.
