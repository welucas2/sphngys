#!/usr/bin/env python
"""
William Lucas, University of St Andrews
Read an sphNG binary dump and make the data accessible to Python.
This is based on my Fortran module rbin.f90 rather than directly on
rdump.F from within the sphNG code itself.

Created: 25 February 2015
Last modified: 17 June 2016

Revision list: Added verbosity option.
Bugs:
Future: - allow for different default real size? specify as an input parameter.
"""

import sys
import colorama
import numpy as np
import fortranfile as ff

class SphngBin:
    def __init__(self, fname=None, contiguous=None, verbosity=None, nsinkmax=None, igradh=None, imhd=None, iexf=None, imigrate=None):
        """
        Called at the time of the object's creation. Sets important parameters
        such as how to treat MPI files and whether gradh SPH was used.
        """

        # Firstly, set the default arguments.
        # fname is the name of the binary file we're going to open
        if fname is None:
            sys.exit("Error in SphngBin, filename must be provided as an argument!")
        # contiguous determines how MPI files are going to be treated:
        #   if true, the whole file will be read and stored at once
        #   if false, the file will be read block-by-block
        if contiguous is None:
            contiguous = True
        # verbosity determines how much information to print to STDOUT while reading.
        #   if verbosity = 0, silent operation: nothing at all is printed
        #   if verbosity = 1, print filenames and blocks
        #   if verbosity = 2, print filenames, blocks and information from arrays
        if verbosity is None:
            verbosity = 1
        elif verbosity < 0 or verbosity > 2:
            sys.exit("verbosity must be 0, 1 or 2")
        # nsinkmax is the size of sink arrays, hence the maximum number which may be stored.
        if nsinkmax is None:
            nsinkmax = 2000
        # igradh indicates if the file was produced with gradh code
        if igradh is None:
            igradh = True
        # imhd indicates if MHD was used in the simulation
        if imhd is None:
            imhd = False
        # iexf indicates which external force was applied; 0 means none
        if iexf is None:
            iexf = 0
        # imigrate indicates if planetesimal... bits were used
        if imigrate is None:
            imigrate = False

        # Now set the class variables to the arguments.
        self.fname = fname
        self.contiguous = contiguous
        self.verbosity = verbosity
        self.nsinkmax = nsinkmax
        self.igradh = igradh
        self.imhd = imhd
        self.iexf = iexf
        self.imigrate = imigrate

        # Sink arrays are created here.
        self.listpm = np.empty(self.nsinkmax, dtype="i4")
        self.spinx = np.empty(self.nsinkmax, dtype="f8")
        self.spiny = np.empty(self.nsinkmax, dtype="f8")
        self.spinz = np.empty(self.nsinkmax, dtype="f8")
        self.angaddx = np.empty(self.nsinkmax, dtype="f8")
        self.angaddy = np.empty(self.nsinkmax, dtype="f8")
        self.angaddz = np.empty(self.nsinkmax, dtype="f8")
        self.angaddz = np.empty(self.nsinkmax, dtype="f8")
        self.spinadx = np.empty(self.nsinkmax, dtype="f8")
        self.spinady = np.empty(self.nsinkmax, dtype="f8")
        self.spinadz = np.empty(self.nsinkmax, dtype="f8")

        # Declare ALL the variables here in the constructor, just because it's good practice (apparently).
        self.binfile = None
        self.fileident = None
        self.tagged = None
        self.nparttot = None
        self.n1 = None
        self.n2 = None
        self.nreassign = None
        self.naccrete = None
        self.nkill = None
        self.nblocks = None
        self.npartblocks = None
        self.iuniquemax = None
        self.gt = None
        self.dtmax = None
        self.gamma = None
        self.rhozero = None
        self.RK2 = None
        self.escap = None
        self.tkin = None
        self.tgrav = None
        self.tterm = None
        self.anglostx = None
        self.anglosty = None
        self.anglostz = None
        self.specang = None
        self.ptmassin = None
        self.tmag = None
        self.Bextx = None
        self.Bexty = None
        self.Bextz = None
        self.hzero = None
        self.uzero_n2 = None
        self.hmass = None
        self.gapfac = None
        self.sdprof = None
        self.rorbit_orig = None
        self.pmrate = None
        self.rorbitmax = None
        self.min_rplan = None
        self.max_rplan = None
        self.planetesimalmass = None
        self.coremass_orig = None
        self.coremass = None
        self.udist = None
        self.umass = None
        self.utime = None
        self.umagfd = None
        self.numberarray = None
        self.icount = None
        self.icountsink = None
        self.iblock = None

        # Variable notes whether potens are available
        self.potenavail = None

        # Open the file to get things going.
        self.__open_file()

    def __verbprint(self, targetverb, *args):
        """
        Prints a list of arguments to screen only if self.verbosity >= targetverb.
        A call with no *args will simply move to the next line.
        """
        if self.verbosity >= targetverb:
            # Print all given arguments on the same line.
            for arg in args:
                print arg,
            # Then move onto a new line.
            print

    def __open_file(self):
        """
        Opens the file as a FortranFile object.
        """
        colorama.init()
        self.__verbprint(1, colorama.Fore.RED + "OPENING BINARY: " + self.fname + colorama.Style.RESET_ALL)
        colorama.deinit()
        try:
            self.binfile = ff.FortranFile(self.fname)
        except IOError:
            sys.exit("Error opening file - check it exists!")
        self.__verbprint(1, "File open")

    def read_header(self):
        """
        Reads the header section of the binary, containing npart, the number of
        MPI blocks, etc.
        """
        colorama.reinit()
        self.__verbprint(1, colorama.Fore.RED + "Starting header read sequence." + colorama.Style.RESET_ALL)
        colorama.deinit()
        self.binfile.readRecord()
        self.fileident = self.binfile.readString()
        self.__verbprint(1, "IDENTITY:", self.fileident.strip())
        if self.fileident[0] != "F":
            sys.exit("File does not contain a full dump, exiting!")
        if self.fileident[1] == "T":
            self.tagged = True
            self.__verbprint(2, "File is tagged.")
        else:
            self.tagged = False
            self.__verbprint(2, "File is NOT tagged.")

        number = self.binfile.readInts()[0]
        if self.tagged:
            self.binfile.readRecord()
        buff = self.binfile.readInts()
        self.nparttot = buff[0]
        self.n1 = buff[1]
        self.n2 = buff[2]
        self.nreassign = buff[3]
        self.naccrete = buff[4]
        self.nkill = buff[5]
        if number == 6:
            self.nblocks = 1
        else:
            self.nblocks = buff[6]
        self.__verbprint(1, "NPARTTOT", self.nparttot, "NBLOCKS", self.nblocks)
        self.npartblocks = np.zeros(self.nblocks, dtype="i4")

        self.binfile.readRecord()
        self.binfile.readRecord()
        self.binfile.readRecord()
        number = self.binfile.readInts()[0]
        if number == 1:
            if self.tagged:
                self.binfile.readRecord()
            self.iuniquemax = self.binfile.readInts()[0]
        else:
            self.iuniquemax = self.nparttot
        self.__verbprint(2, "iuniquemax", self.iuniquemax)

        number = self.binfile.readInts()[0]
        if self.tagged:
            self.__verbprint(2, "Reading", number, "tags from file.")
            tagsreal = self.binfile.readString().split()
        else:
            self.__verbprint(2, "Simulating tags for non-tagged file.")
            self.__verbprint(2, "number =", number)
            tagsreal = ["gt", "dtmax", "gamma", "rhozero", "RK2", "escap",
                        "tkin", "tgrav", "tterm", "anglostx", "anglosty",
                        "anglostz", "specang", "ptmassin", "tmag", "Bextx",
                        "Bexty", "Bextz", "hzero", "uzero_n2", "hmass",
                        "gapfac", "pmassinitial", "sdprof", "rorbit_orig",
                        "min_rplan", "max_rplan", "planetesimalmasscoremass_orig",
                        "coremass"]
        self.__verbprint(2, "TAGS:", tagsreal)

        # Read the block of reals containing header data. Note we have
        # to specify double precision here.
        rheader = self.binfile.readReals("d")
        try:
            self.gt = rheader[tagsreal.index("gt")]
        except ValueError:
            sys.exit("gt not found in header, exiting")
        self.__verbprint(2, "Time is", self.gt)
        try:
            self.dtmax = rheader[tagsreal.index("dtmax")]
        except ValueError:
            sys.exit("dtmax not found in header, exiting")
        try:
            self.gamma = rheader[tagsreal.index("gamma")]
        except ValueError:
            sys.exit("gamma not found in header, exiting")
        try:
            self.rhozero = rheader[tagsreal.index("rhozero")]
        except ValueError:
            sys.exit("rhozero not found in header, exiting")
        try:
            self.RK2 = rheader[tagsreal.index("RK2")]
        except ValueError:
            sys.exit("RK2 not found in header, exiting")
        self.escap = rheader[tagsreal.index("escap")]
        self.tkin = rheader[tagsreal.index("tkin")]
        self.tgrav = rheader[tagsreal.index("tgrav")]
        self.tterm = rheader[tagsreal.index("tterm")]
        self.anglostx = rheader[tagsreal.index("anglostx")]
        self.anglosty = rheader[tagsreal.index("anglosty")]
        self.anglostz = rheader[tagsreal.index("anglostz")]
        self.specang = rheader[tagsreal.index("specang")]
        self.ptmassin = rheader[tagsreal.index("ptmassin")]
        if self.imhd:
            self.tmag = rheader[tagsreal.index("tmag")]
            self.Bextx = rheader[tagsreal.index("Bextx")]
            self.Bexty = rheader[tagsreal.index("Bexty")]
            self.Bextz = rheader[tagsreal.index("Bextz")]
            self.__verbprint(2, "External field found, Bextx =", self.Bextx, self.Bexty, self.Bextz)
        self.hzero = rheader[tagsreal.index("hzero")] if len(rheader) > tagsreal.index("hzero") else 0.
        self.uzero_n2 = rheader[tagsreal.index("uzero_n2")] if len(rheader) > tagsreal.index("uzero_n2") else 0.
        if self.uzero_n2 > 0.:
            self.__verbprint(2, "u for surrounding medium =", self.uzero_n2)
        self.hmass = rheader[tagsreal.index("hmass")] if len(rheader) > tagsreal.index("hmass") else 0.
        self.gapfac = rheader[tagsreal.index("gapfac")] if len(rheader) > tagsreal.index("gapfac") else 0.
        try:
            self.sdprof = rheader[tagsreal.index("sdprof")] if len(rheader) > tagsreal.index("sdprof") else 0.
        except ValueError:
            self.sdprof = -0.5
            self.__verbprint(2, "Surface density pre vary (goes as r^-0.5)")
        self.rorbit_orig = rheader[tagsreal.index("rorbit_orig")] if len(rheader) > tagsreal.index("rorbit_orig") else 0.
        # CHECK THE NEXT BLOCK OF CODE IF ACTUALLY USING PLANETESIMALS
        if self.imigrate:
            self.rorbitmax = (self.rorbitmax - self.rorbit_orig) / self.pmrate
        # I'm actually going to lock this whole section behind an if -
        # delete if it's ever actually needed.
        if self.imigrate:
            try:
                self.min_rplan = rheader[tagsreal.index("min_rplan")]
                self.__verbprint(2, "r_min =", self.min_rplan)
            except IndexError:
                pass
            try:
                self.max_rplan = rheader[tagsreal.index("max_rplan")]
                self.__verbprint(2, "r_max =", self.max_rplan)
            except IndexError:
                pass
            if not self.tagged:
                # At the moment this prints no matter the verbosity. Change in future
                # if it's a problem.
                self.__verbprint(0, "DO NOT USE CODE WITH DUMPS MADE FROM APRIL 2011")
                self.__verbprint(0, "AND BEFORE 16/05/2011. For these, reals 26, 27 are")
                self.__verbprint(0, "planetesimal radius and density respectively.")
            try:
                self.planetesimalmass = rheader[tagsreal.index("planetesimalmass")]
                self.__verbprint(2, "Planetesimal mass", self.planetesimalmass)
            except IndexError:
                pass
            try:
                self.coremass_orig = rheader[tagsreal.index("coremass_orig")]
                self.__verbprint(2, "Core mass orig", self.coremass_orig)
            except IndexError:
                pass
            try:
                self.coremass = rheader[tagsreal.index("coremass")]
                self.__verbprint(2, "Core mass running record", self.coremass)
            except IndexError:
                pass

        self.binfile.readRecord()
        number = self.binfile.readInts()[0]
        if number < 3:
            sys.exit("Error in rbin, nreal8 too small in header section")
        if self.tagged:
            self.binfile.readRecord()
        buff = self.binfile.readReals("d")
        self.udist = buff[0]
        self.umass = buff[1]
        self.utime = buff[2]
        if self.imhd:
            if number > 3:
                self.umagfd = buff[3]
            else:
                self.__verbprint(0, "WARNING: no mag field units in rdump")
        self.__verbprint(2, "Distance, mass, time units are:", self.udist, self.umass, self.utime)

        number = self.binfile.readInts()[0]
        self.numberarray = number / self.nblocks
        self.__verbprint(2, "Array types", number, self.numberarray, "\n")

        # AT THIS POINT, WE ARE AT THE START OF THE DATA
        self.icount = 0
        self.icountsink = 0
        # AND AT THE FIRST MPI BLOCK
        self.iblock = 0

        return 0

    def __create_hydro_arrays(self, nentries):
        """
        Create the arrays to hold particle data. We don't create them on the fly
        as if reading multiple blocks at once it's preferable to fill them up
        slowly rather than constantly make them larger whenever we move onto
        the next block.
        """
        self.__verbprint(2, "ALLOCATING HYDRO ARRAYS WITH", nentries, "ENTRIES")
        self.isteps = np.empty(nentries, dtype="i4")
        self.iphase = np.empty(nentries, dtype="i1")
        self.iunique = np.empty(nentries, dtype="i8")
        self.xyzmh = np.empty((nentries,5), dtype="f8")
        self.vxyzu = np.empty((nentries,4), dtype="f8")
        self.rho = np.empty(nentries, dtype="f8")
        self.alphaMM = np.empty(nentries, dtype="f4")
        self.poten = np.empty(nentries, dtype="f4")
        if self.igradh:
            self.gradh = np.empty(nentries, dtype="f4")
            self.gradhs = np.empty(nentries, dtype="f4")
        else:
            self.dgrav = np.empty(nentries, dtype="f4")

    def __create_rt_arrays(self, nentries):
        """
        Create RT arrays much the same way as in __create_hydro_arrays.
        """
        self.__verbprint(2, "ALLOCATING RT ARRAYS WITH", nentries, "ENTRIES")
        self.nneigh = np.empty(nentries, dtype="i2")
        self.e = np.empty(nentries, dtype="f8")
        self.rkappa = np.empty(nentries, dtype="f8")
        self.cv = np.empty(nentries, dtype="f8")
        self.rlambda = np.empty(nentries, dtype="f8")
        self.edd = np.empty(nentries, dtype="f8")
        self.force = np.empty((nentries,3), dtype="f8")
        self.dlnTdlnP = np.empty(nentries, dtype="f4")
        self.adiabaticgradient = np.empty(nentries, dtype="f4")
        self.pressure = np.empty((nentries,3), dtype="f4")
        self.viscosity = np.empty((nentries,3), dtype="f4")
        self.gravity = np.empty((nentries,3), dtype="f4")
        self.radpres = np.empty((nentries,3), dtype="f4")

    def __read_block(self):
        """
        Read a single block of data into memory.
        """
        colorama.reinit()
        self.__verbprint(1, colorama.Fore.RED + "Reading block", self.iblock, colorama.Style.RESET_ALL)
        colorama.deinit()

        # The record is mixed between integer*8 for number8 and integer*4 for the nums,
        # so we need to extract them manually from the record.
        buff = self.binfile.readRecord()
        number8 = np.fromstring(buff[0:8], dtype=self.binfile.ENDIAN + "q")[0]
        nums = np.fromstring(buff[8:40], dtype=self.binfile.ENDIAN + "i")
        self.__verbprint(2, self.nparttot, self.nblocks, self.iblock, number8, number8 * self.nblocks)
        self.npart = number8
        self.npartblocks[self.iblock] = number8
        self.__verbprint(2, "Block contains", number8, "particles, current icount", self.icount)
        self.__verbprint(2, "nums", nums)

        # Number of sinks
        buff = self.binfile.readRecord()
        number8 = np.fromstring(buff[0:8], dtype=self.binfile.ENDIAN + "q")[0]
        numssink = np.fromstring(buff[8:40], dtype=self.binfile.ENDIAN + "i")
        self.__verbprint(2, "numssink", numssink)
        self.nptmass = number8

        # Radiative transfer bits
        if self.numberarray == 3:
            buff = self.binfile.readRecord()
            number8 = np.fromstring(buff[0:8], dtype=self.binfile.ENDIAN + "q")[0]
            numsRT = np.fromstring(buff[8:40], dtype=self.binfile.ENDIAN + "i")
            self.__verbprint(2, "numsRT", numsRT)

        ic1 = self.icount
        ic2 = self.icount + 3

        # If this is the first block, the arrays haven't yet been created.
        if not self.contiguous or self.iblock == 0:
            if self.contiguous and self.nblocks > 1:
                # In this case, the arrays will hold the data for ALL particles
                nhydro = self.nparttot
            else:
                # while here, they will only hold particles from the current block.
                nhydro = self.npart
            # Then create the arrays.
            self.__create_hydro_arrays(nhydro)

        # Finally, get reading the actually useful data!
        # Start with the timesteps
        if self.tagged:
            tagi = self.binfile.readRecord()
        self.isteps[self.icount:self.icount + self.npart] = self.binfile.readInts(prec="i")
        self.__verbprint(2, "isteps", self.isteps[ic1:ic2])

        # Skip reading listinactive
        if nums[0] >= 2:
            if self.tagged:
                self.binfile.readRecord()
            self.__verbprint(2, "Skipping listinactive")
            self.binfile.readRecord()

        # Read particle phases
        # Another bit of trickery is needed since iphase is integer*1
        if self.tagged:
            tagi = self.binfile.readRecord()
        buff = self.binfile.readRecord()
        self.iphase[self.icount:self.icount + self.npart] = np.fromstring(buff, dtype=self.binfile.ENDIAN + "b")
        self.__verbprint(2, "iphase", self.iphase[ic1:ic2])

        # Particle iuniques, only in versions post-MPI incorporation I believe
        if nums[4] >= 1:
            if self.tagged:
                self.binfile.readRecord()
            buff = self.binfile.readRecord()
            self.iunique[self.icount:self.icount + self.npart] = np.fromstring(buff, dtype=self.binfile.ENDIAN + "q")
            self.__verbprint(2, "iunique", self.iunique[ic1:ic2])

        # Position, mass, smoothing length
        # NB the arrays are now stored in row-major order (Python/C) in order
        # to ensure interoperability with other Python code, even though numpy
        # arrays CAN be switched to column-major (Fortran) mode.
        for i in range(5):
            if self.tagged:
                self.binfile.readRecord()
            self.xyzmh[self.icount:self.icount + self.npart, i] = self.binfile.readReals("d")
        self.__verbprint(2, "x:", self.xyzmh[ic1:ic2, 0])
        self.__verbprint(2, "m:", self.xyzmh[ic1:ic2, 3])
        self.__verbprint(2, "h:", self.xyzmh[ic1:ic2, 4])

        # Velocity, internal energy
        for i in range(4):
            if self.tagged:
                self.binfile.readRecord()
            self.vxyzu[self.icount:self.icount + self.npart, i] = self.binfile.readReals("d")
        self.__verbprint(2, "vx:", self.vxyzu[ic1:ic2, 0])
        self.__verbprint(2, "u:", self.vxyzu[ic1:ic2, 3])

        # Rho
        if self.tagged:
            self.binfile.readRecord()
        self.rho[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
        self.__verbprint(2, "rho:", self.rho[ic1:ic2])
        iread = 1

        # gradh and gradhs if sphNG ran in gradh mode;
        # dgrav if not.
        if self.igradh:
            if nums[6] >= 2:
                if self.tagged:
                    self.binfile.readRecord()
                self.gradh[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
                self.__verbprint(2, "gradh:", self.gradh[ic1:ic2])
                iread = iread + 1
            if nums[6] >= 3:
                if self.tagged:
                    tagi = self.binfile.readString().strip()
                self.gradhs[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
                if self.tagged:
                    # Following line is only good for Python 2.7+
                    # print "{}:".format(tagi), self.gradhs[ic1:ic2]
                    self.__verbprint(2, tagi.replace(" ", "") + ":", self.gradhs[ic1:ic2])
                else:
                    self.__verbprint(2, "gradhs:", self.gradhs[ic1:ic2])
                iread = iread + 1
        else:
            if self.tagged:
                self.binfile.readRecord()
            self.dgrav[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
            self.__verbprint(2, "dgrav:", self.dgrav[ic1:ic2])
            iread = 2

        # alphaMM
        if (nums[6] > iread):
            if self.tagged:
                self.binfile.readRecord()
            self.alphaMM[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
            self.__verbprint(2, "alphaMM:", self.alphaMM[ic1:ic2])

        # potens for Duncan, if they exist
        if nums[6] >= 5:
            if self.tagged:
                self.binfile.readRecord()
            self.poten[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
            self.__verbprint(2, "poten:", self.poten[ic1:ic2])
            self.potenavail = True
        else:
            self.potenavail = False

        # SINK PARTICLE DATA
        self.__verbprint(2, "Reading sink particle data for", self.nptmass, "sinks.")
        if self.tagged:
            self.binfile.readRecord()
        self.listpm[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readInts()
        self.listpm = self.listpm - 1  # Change to Python indexing starting at 0.
        if self.tagged:
            self.binfile.readRecord()
        self.spinx[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.spiny[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.spinz[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.angaddx[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.angaddy[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.angaddz[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.spinadx[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.spinady[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")
        if self.tagged:
            self.binfile.readRecord()
        self.spinadz[self.icountsink:self.icountsink + self.nptmass] = self.binfile.readReals("d")

        # This point onwards is a mystery; I've commented out the next few lines as they seem to
        # not be represented in a 'modern' rdump.F.
        #  print "numssink[5] =", numssink[5]
        #  for i in range(numssink[5] - 9):
        #      self.binfile.readRecord()
        for i in range(numssink[7]):
            if self.tagged:
                self.binfile.readRecord()
            self.binfile.readRecord()

        # THE NEXT SECTION OF CODE IS ONLY USED WITH RT RUNS and IS NOT updated for tags.
        # I think it's been wrong since the old original F77 code and no one's bothered with
        # it because none of us use RT.
        if self.numberarray == 3:
            self.__verbprint(2, "Reading RT data...", self.numberarray)
            if numsRT[2] == 1:
                self.nneigh[self.icount:self.icount + self.npart] = self.binfile.readInts("h")
            self.e[self.icount:self.icount + self.npart] = self.binfile.readReals("d")
            self.rkappa[self.icount:self.icount + self.npart] = self.binfile.readReals("d")
            self.cv[self.icount:self.icount + self.npart] = self.binfile.readReals("d")
            self.rlambda[self.icount:self.icount + self.npart] = self.binfile.readReals("d")
            self.edd[self.icount:self.icount + self.npart] = self.binfile.readReals("d")
            if numsRT[5] == 8:
                for i in range(3):
                    self.force[self.icount:self.icount + self.npart, i] = self.binfile.readReals("d")
            if numsRT[6] == 1:
                self.dlnTdlnP[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
            if numsRT[6] == 13 or numsRT[6] == 14:
                self.dlnTdlnP[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
                if numsRT[6] == 11:
                    self.adiabaticgradient[self.icount:self.icount + self.npart] = self.binfile.readReals("f")
                for i in range(3):
                    self.pressure[self.icount:self.icount + self.npart, i] = self.binfile.readReals("f")
                for i in range(3):
                    self.viscosity[self.icount:self.icount + self.npart, i] = self.binfile.readReals("f")
                for i in range(3):
                    self.gravity[self.icount:self.icount + self.npart, i] = self.binfile.readReals("f")
                for i in range(3):
                    self.radpres[self.icount:self.icount + self.npart, i] = self.binfile.readReals("f")

        # Finally, we correctly update the icounts and iblock.
        if self.contiguous:
            self.icount = self.icount + self.npart
            self.icountsink = self.icountsink + self.nptmass
        else:
            self.icount = 0
            self.icountsink = 0
        self.iblock += 1

        self.__verbprint(2, "Block finished.\n")

        return 0

    def read_data(self):
        """
        Get actual particle data from the binary. This function only acts as a
        wrapper for __readBlock, and will call it once or multiple times 
        depending on whether self.contiguous is set true or false and whether 
        more than one block actually exists.
        """
        if self.iblock >= self.nblocks:
            self.__verbprint(0, "Error: entire file has already been read!")
            return 1
        if not self.contiguous or self.nblocks == 1:
            # If reading only one block at a time, or if only one block exists,
            # then only one call needs to be made.
            ierr = self.__read_block()
        else:
            # If multiple blocks exist and the mode IS contiguous, read the entire
            # file, storing each blocks data in the arrays successively.
            for i in range(self.nblocks):
                ierr = self.__read_block()
                if ierr != 0:
                    sys.exit("Error reading block " + str(i))
            # Copy icountsink into nptmass for easier use
            self.nptmass = self.icountsink
        self.__verbprint(1)
        self.__verbprint(2)
        return ierr

if __name__ == "__main__":
    print "SPHNGYS DEMO CODE\n"
    ftest = "/Users/william/PycharmProjects/density_on_sky/data/RUNI180"
    # ftest = "../data/EMCD007"  # tagged file
    # ftest = "../data/TES2497"
    # ftest = "CLOD010"  # normal file
    # ftest = "../data/GCBG061"  # MPI file
    testfile = SphngBin(fname=ftest, contiguous=True)
    iout = testfile.read_header()
    print "readHeader returned", iout

    iout = testfile.read_data()
    print "readData returned", iout
