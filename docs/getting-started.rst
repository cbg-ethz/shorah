Getting started
===============

There are two options to set up the environment required for ``shorah``:

* Docker

    Navigate to the project's root directory and run the following command to
    start a console inside the container.

    .. code-block:: bash
        
        # the following will build and run a Docker container with shorah
        docker run --rm -w="/usr/app" -it $(docker build -q .) bash
        # once inside the container:
        poetry run shorah -h # to get the help
        # the for example:
        poetry run shorah shotgun <args...> # replace <args...> 

* On your own machine 

    Make sure you have the following software installed locally:

        * Poetry (a Python package manager)
        * HTSlib
        * Boost C++ library (only ``random`` is necessary)

    .. code-block:: bash

        # in the project's root directory
        poetry install
        poetry run shorah -h # to get the help
        # the for example:
        poetry run shorah shotgun <args...> # replace <args...> 