cvolume
=======

``cvolume`` is a `SageMath <https://www.sagemath.org>`_ module to compute completed
and Masur-Veech volumes of strata of quadratic differentials with odd zeros.
You can install it on top of a SageMath installation on your computer (see the instructions
below). 

``cvolume`` is based on `SageMath <https://www.sagemath.org>`_ module 
_`admcycles <https://gitlab.com/jo314schmitt/admcycles>`_
to compute intersection numbers on moduli spae of complex curves. This
README file is partially borrowed from ``admcycles`` page.

Installation
------------

Prerequisites
^^^^^^^^^^^^^

Make sure that `SageMath <https://www.sagemath.org>`_ is installed on your
computer. Detailed installation instructions for different operating systems
are available `here
<http://doc.sagemath.org/html/en/installation/binary.html>`_ and on the
SageMath website.

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
  <https://github.com/eduryev/cvolume/archive/master.zip>`__.

- Inside a shell console run::

      $ sage -pip install /where/is/the/package.zip --user

  Here, the ``--user`` is an optional argument to ``pip`` which, when
  provided, will install admcycles not inside the Sage folder but in your home
  directory.
