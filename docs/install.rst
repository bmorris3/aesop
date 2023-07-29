.. _install:

.. include:: references.txt

****************
Installing aesop
****************

Clone the ``aesop`` repository from GitHub, and go into the top-level ``aesop``
directory::

    git clone https://github.com/bmorris3/aesop
    cd aesop

Install ``aesop`` with::

    pip install .

You can ensure that the package was installed successfully by running the tests
from the command line with::

    tox -e test
