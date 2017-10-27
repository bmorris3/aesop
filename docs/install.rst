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

    python setup.py install

You can ensure that the package was installed successfully by running the tests
from a Python interpreter with the following commands:

.. code-block:: python

    >>> import aesop
    >>> aesop.test()
