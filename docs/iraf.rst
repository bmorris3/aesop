.. include:: references.txt

.. _iraf:

****
IRAF
****

Contents
========

* :ref:`iraf-down`
* :ref:`iraf-prep`
* :ref:`iraf-reduce`

.. _iraf-down:

Downloading and installing AstroConda with IRAF
-----------------------------------------------
Go to the `AstroConda installation page <http://astroconda.readthedocs.io/en/latest/installation.html?highlight=iraf#legacy-software-stack-with-iraf>`_,
and follow the instructions to download AstroConda with IRAF. Unfortunately this only works with Python 2.7.

.. _iraf-prep:

Preparing your data
-------------------

``ReduceARCES`` will expect your flats in the red to start with "redflat",
your flats in the blue to start with "blueflat", your biases to start with "bias", and your ThAr lamp exposures to start
with "ThAr".

.. _iraf-reduce:

Reducing your data
------------------
First, make sure your iraf27 (or whatever you named it when installing AstroConda with IRAF) environment is activated by
doing ::

    source activate iraf27

from a bash terminal.

Copy all of the files in the ``aesop/iraf`` directory into the directory where your raw data and cals are. From that
directory, do ::

    mkiraf

which will create ``login.cl`` in that directory. If prompted to specify a terminal, say ``xgterm``.
Edit ``login.cl`` to include the following imports ::

    imred
    ccdred
    echelle
    crutil
    astutil
    images
    imgeom
    twodspec
    apextract
    onedspec

Next, open an xgterm by doing ``xgterm``, and do ::

    cl
    cl < ReduceARCES.cl

This will extract the spectrum from your data. When prompted to edit the parameters of the tasks that ``ReduceARCES.cl``
uses, simply do ``:q``. You'll have to do this twice at the beginning of the reduction and a few more times in the middle.

Now in order to use ``aesop``, follow the tutorial in :ref:`getting_started`