cvolume
=======

``cvolume`` is a `SageMath <https://www.sagemath.org>`_ module to compute completed
and Masur-Veech volumes of strata of quadratic differentials with odd zeros.
You can install it on top of a SageMath installation on your computer (see the instructions
below). 

``cvolume`` is based on `SageMath <https://www.sagemath.org>`_ module 
`admcycles <https://gitlab.com/jo314schmitt/admcycles>`_
to compute intersection numbers on moduli space of complex curves. This
README file is partially borrowed from ``admcycles`` page.

Installation
------------

Prerequisites
^^^^^^^^^^^^^

Make sure that `SageMath 9 <https://www.sagemath.org>`_ is installed on your
computer. Detailed installation instructions for different operating systems
are available `here
<http://doc.sagemath.org/html/en/installation/binary.html>`_ and on the
SageMath website.

WARNING: Module will not work with SageMath 8 as that version doesn't support some
functionality of Python 3. Make sure to use the last version of SageMath!

All the command below must be run inside a console (the last character of the
prompt is usally a ``$``). On Windows this is called ``SageMath Shell`` while
on Linux and MacOS this is often referred to as a ``Terminal``.

Inside the console, we assume that the command ``sage`` launches a Sage
session (whose prompt is usually ``sage:``). To quit the Sage session
type ``quit`` and press ``Enter``.

Manual installation with pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the procedure below:

- Download the admcycles package as a ``.zip`` file from `github
  <https://github.com/eduryev/cvolume/archive/v1.0.zip>`__.

- Inside a shell console run::

      $ sage -pip install /where/is/the/package.zip --user

  Here, the ``--user`` is an optional argument to ``pip`` which, when
  provided, will install admcycles not inside the Sage folder but in your home
  directory.
 
Cloning with git
^^^^^^^^^^^^^^^^

You can also clone this repository to your computer, if you have
versioning software ``git`` installed. Go to the folder where you
want ``cvolume`` to be stored and run the following commands::

    $ git init
    $ git clone https://github.com/eduryev/cvolume.git
    $ cd cvolume
    $ git fetch --tags
    $ git checkout tags/v1.0 -b v1.0
    
Use without installation
^^^^^^^^^^^^^^^^^^^^^^^^

To use the package without installing, download the package as a ``.zip`` file either
from `github
<https://github.com/eduryev/cvolume/archive/v1.0.zip>`__.
Unpack the ``.zip`` file. This creates a folder which should
contain a file ``setup.py``. In order to use the
module, you need to run Sage from this folder. For example, if the full path of
this folder is ``/Amazing/Organized/Folder/System/cvolume``, you could do::

    $ cd /Amazing/Organized/Folder/System/cvolume
    $ sage

Or directly inside a Sage session::

    sage: cd /Amazing/Organized/Folder/System/admcycles


Examples
-------

To start using cvolume, start a Sage session: either in the command line::

    $ sage
  
or a Jupyter notebook, or inside one of the online services). Then type::

    sage: from cvolume import *

To try a first computation, you can compute the completed volume of Q(1,1,1,1), 
which coincides with its Masur-Veech volume, since the stratum is principle::

    sage: completed_volume([1,1,1,1])
    1/15*pi^6

You can also compute the completed volume of Q(3,1), which has Masur-Veech volume 0,
since it is empty::

    sage: completed_volume([3,1])
    23/90*pi^4
  
Apart from computing volumes you can generate labeled stable graph. Here is a an example,
which generated all of labeled stable graphs for stratum [3,1,1,-1]::

    sage: stable_lab_graphs([3,1,1,-1])
    {Labeled Stable Graph with edges = ((0, 1, 1),), loops = (0, 0), kappa = ((1, 1), (3, -1)),
     Labeled Stable Graph with edges = ((0, 1, 1),), loops = (0, 1), kappa = ((1, 1), (3, -1)),
     Labeled Stable Graph with edges = ((0, 1, 1),), loops = (1, 0), kappa = ((1, 1), (3, -1)),
     Labeled Stable Graph with edges = ((0, 1, 1),), loops = (1, 1), kappa = ((1, 1), (3, -1)),
     Labeled Stable Graph with edges = ((0, 1, 2),), loops = (0, 0), kappa = ((1, -1), (3, 1)),
     Labeled Stable Graph with edges = ((0, 1, 2),), loops = (0, 1), kappa = ((1, -1), (3, 1)),
     Labeled Stable Graph with edges = ((0, 1, 3),), loops = (0, 0), kappa = ((1, 1), (3, -1)),
     Labeled Stable Graph with edges = (), loops = (1,), kappa = ((3, 1, 1, -1),),
     Labeled Stable Graph with edges = (), loops = (2,), kappa = ((3, 1, 1, -1),)}




