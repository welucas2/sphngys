# README #

sphngys is a Python class with which you can quickly read and store data for use from an sphNG binary file.

### Requirements ###

* Python 2.6+ (no idea about Python 3)

* numpy

* fortranfile

* colorama (optionally - helps with interpreting prints to STDOUT, comment out if you don't care)

### How to use ###

sphngys.py should be placed inside the working directory, or else placed in your PYTHONPATH environment variable. Make sure if changing it to similarly force a change in the compiled sphngys.pyc file by deleting it or otherwise!

Import to your Python program with:
from sphngys import SphngBin

Then, create a reference to a SphngBin object with, eg.:
binfile = SphngBin(fname = "TEST001")

fname should always be provided to SphngBin. Optional arguments are:

* contiguous - logical, default True - see below

* nsinkmax - integer, default 2000 - size of sink arrays

* verbosity - 0, 1, 2 for increasing levels of verbosity - 0 for silent running

* igradh - logical, default True - whether to expect grad-h arrays

* imhd - logical, default False - whether to expect MHD arrays

* iexf - integer, default 0 - which external force was used, 0 meaning none

* imigrate - logical, default False - whether planetesimal stuff was going on

Care should be taken with the latter three arguments, as they aren't implemented, although the code still has to allow for the possibility. Generally, always leave them with their default values.

Once the variable referencing a SphngBin instance is there, you may read the header with e.g.:
binfile.read_header()
and the data block(s) with:
binfile.read_data()

The logical variable 'contiguous' should be considered when working with files using MPI blocks. When True, the class will read the MPI blocks sequentially, storing the arrays from each block in single contiguous arrays. As such, the entire dataset can be accessed without further read operations. If the particle number is large enough that your machine is not able to read the entire file, set contiguous to False. In this case, control will pass back to you after each block is read. In this case you will be reading one block, performing your calculations, then reading the next block, and so on - that is to say that in the example above, one binfile.read_data() call would be made for each MPI block.

### For the future ###

Create a new method to write sphNG binaries. This should allow use of Python as a potentially costly but simple replacement for extract.f or even setpart.F.

### Who do I talk to? ###

Email wel2@st-andrews.ac.uk for help.