.. _install:

Installation
============

Prior to using pyRVT, Python and a number of packages need to be installed. In
addition to Python, the following packages need to be installed:

*Required:*

- numpy -- fast vector operations

- scipy -- indefinite integration

- nose -- required for the unit tests

*Optional:*

-  matplotlib -- used for plotting

-  xlrd/xlwt -- reading and writing xls files

-  openpyxl -- reading xlsx files

Install Python dependencies is best accomplished with a package manager. On
Windows or OS-X, I recommend using Miniconda.

Miniconda has installers for `Windows 32-bit`_, `Windows 64-bit`_, and `OS-X`_.
See the full list of installers `here` -- make sure you select Miniconda3 for
Python3.

.. _Windows 32-bit: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-Windows-x86.exe
.. _Windows 64-bit: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-Windows-x86_64.exe
.. _OS-X: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-MacOSX-x86_64.sh

After the installer is finished, install the required dependencies by opening a
terminal. On Windows, this is best accomplished with ``Windows Key + r``, enter
"cmd". Next enter the following command:

::
 
  conda install --yes setuptools numpy scipy matplotlib nose openpyxl xlrd pip

On Windows, the text can copied and pasted if "Quick Edit" mode is enabled. To
enable this feature, right click on the icon in the upper left portion of the
window, and select "Properties", and then check the "Quick Edit Mode" check box
within the "Edit Options" group. Copy the text, and then paste it by click the
right mouse button.

Next, install or upgrade pyRVT:

::

  pip install --upgrade https://github.com/arkottke/pyrvt/archive/master.zip
