##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates parsing functions of pyReefCore XmL input file.
"""
import os
import glob
import numpy
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict
from decimal import Decimal

class xmlParser:
    """
    This class defines XmL input file variables.

    Parameters
    ----------
    string : inputfile
        The XmL input file name.
    """

    def __init__(self, inputfile = None, makeUniqueOutputDir=True):
        """
        If makeUniqueOutputDir is set, we create a uniquely-named directory for
        the output. If it's clear, we blindly accept what's in the XML file.
        """

        if inputfile==None:
            raise RuntimeError('XmL input file name must be defined to run a pyReef simulation.')
        if not os.path.isfile(inputfile):
            raise RuntimeError('The XmL input file name cannot be found in your path.')
        self.inputfile = inputfile

        self.tStart = None
        self.tEnd = None
        self.tCarb = None
        self.laytime = None

        self.depth0 = None
        self.speciesNb = None
        self.speciesName = None
        self.malthusParam = None
        self.speciesPopulation = None
        self.maxpop = None
        self.speciesProduction = None
        self.communityMatrix = None

        self.seaOn = False
        self.seaval = 0.
        self.seafile = None

        self.flowOn = False
        self.flowval = 0.
        self.flowfile = None
        self.flowfunc = None
        self.flowdecay = None
        self.flowlinb = None
        self.flowlina = None
        self.flowdepth = None

        self.sedOn = False
        self.sedval = 0.
        self.sedfile = None
        self.sedfunc = None
        self.seddecay = None
        self.sedlinb = None
        self.sedlina = None
        self.seddepth = None

        self.enviDepth = None
        self.enviSed = None
        self.enviFlow = None

        self.makeUniqueOutputDir = makeUniqueOutputDir
        self.outDir = None

        #self.h5file = 'h5/surf.time'
        #self.xmffile = 'xmf/surf.time'
        #self.xdmffile = 'surf.series.xdmf'
        self._get_XmL_Data()

        return

    def _get_XmL_Data(self):
        """
        Main function used to parse the XmL input file.
        """

        # Load XmL input file
        tree = ET.parse(self.inputfile)
        root = tree.getroot()

        # Extract time structure information
        time = None
        time = root.find('time')
        if time is not None:
            element = None
            element = time.find('start')
            if element is not None:
                self.tStart = float(element.text)
            else:
                raise ValueError('Error in the definition of the simulation time: start time declaration is required')
            element = None
            element = time.find('end')
            if element is not None:
                self.tEnd = float(element.text)
            else:
                raise ValueErself.ror('Error in the definition of the simulation time: end time declaration is required')
            if self.tStart > self.tEnd:
                raise ValueError('Error in the definition of the simulation time: start time is greater than end time!')
            element = None
            element = time.find('tcarb')
            if element is not None:
                self.tCarb = float(element.text)
            else:
                raise ValueError('Error in the definition of the simulation time: simulation time step for carbonate module is required')
            element = None
            element = time.find('laytime')
            if element is not None:
                self.laytime = float(element.text)
            else:
                self.laytime = self.tCarb
            if Decimal(self.laytime) % Decimal(self.tCarb) != 0.:
                raise ValueError('Error in the XmL file: stratal layer interval needs to be an exact multiple of the carbonate interval!')
            if Decimal(self.tEnd-self.tStart) % Decimal(self.laytime) != 0.:
                raise ValueError('Error in the XmL file: layer time interval needs to be an exact multiple of the simulation time interval!')
        else:
            raise ValueError('Error in the XmL file: time structure definition is required!')

        # Extract habitats structure information
        litho = None
        litho = root.find('habitats')
        if litho is not None:
            element = None
            element = litho.find('depth')
            if element is not None:
                self.depth0 = float(element.text)
                if self.depth0<0:
                    raise ValueError('Error in the initial depth needs to be positive!')
            else:
                raise ValueError('Error in the definition of the initial depth is required!')
            element = None
            element = litho.find('speciesNb')
            if element is not None:
                self.speciesNb = int(element.text)
            else:
                raise ValueError('Error you need to define at least 1 species!')
            element = None
            element = litho.find('maxPopulation')
            if element is not None:
                self.maxpop = int(element.text)
            else:
                self.maxpop = 20
            self.speciesName = numpy.empty(self.speciesNb, dtype="S14")
            self.malthusParam = numpy.zeros(self.speciesNb, dtype=float)
            self.speciesPopulation = numpy.zeros(self.speciesNb, dtype=float)
            self.speciesProduction = numpy.zeros(self.speciesNb, dtype=float)
            self.communityMatrix = numpy.zeros((self.speciesNb,self.speciesNb), dtype=float)
            id = 0
            for facies in litho.iter('species'):
                if id >= self.speciesNb:
                    raise ValueError('The number of species does not match the number of defined ones.')
                element = None
                element = facies.find('name')
                if element is not None:
                    self.speciesName[id] = element.text
                else:
                    raise ValueError('Species name %d is missing in the habitats structure.'%id)
                element = None
                element = facies.find('malthus')
                if element is not None:
                    self.malthusParam[id] = float(element.text)
                else:
                    raise ValueError('Malthusian parameter for species %d is missing in the habitats structure.'%id)
                element = None
                element = facies.find('population')
                if element is not None:
                    self.speciesPopulation[id] = float(element.text)
                else:
                    raise ValueError('Initial population for species %d is missing in the habitats structure.'%id)
                element = None
                element = facies.find('production')
                if element is not None:
                    self.speciesProduction[id] = float(element.text)
                else:
                    raise ValueError('Production for species %d is missing in the habitats structure.'%id)
                id += 1

            if id != self.speciesNb:
                raise ValueError('The number of species declared does not match the number of defined ones.')
            element = None
            element = litho.find('communityMatrix')
            if element is not None:
                data = defaultdict(list)
                # Group into rows of (col, val) tuples
                for val in litho.iter('value'):
                    data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                # Sort columns and format into a space separated string
                rows = []
                for row in data:
                    rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                # Build array from matrix string
                self.communityMatrix = numpy.array(numpy.mat(';'.join(rows)))
            else:
                raise ValueError('Error definition of the community matrix interaction is missing in the habitats structure!')
        else:
            raise ValueError('Error in the XmL file: habitats structure definition is required!')

        # Extract sea-level structure information
        sea = None
        sea = root.find('sea')
        if sea is not None:
            self.seaOn = True
            element = None
            element = sea.find('val')
            if element is not None:
                self.seaval = float(element.text)
            else:
                self.seaval = 0.
            element = None
            element = sea.find('curve')
            if element is not None:
                self.seafile = element.text
                if not os.path.isfile(self.seafile):
                    raise ValueError('Sea level file is missing or the given path is incorrect.')
            else:
                self.seafile = None
        else:
            self.seapos = 0.
            self.seafile = None

        # Extract ocean flow information
        flow = None
        flow = root.find('flow')
        if flow is not None:
            self.flowOn = True
            element = None
            element = flow.find('val')
            if element is not None:
                self.flowval = float(element.text)
            else:
                self.flowval = 0.
            element = None
            element = flow.find('curve')
            if element is not None:
                self.flowfile = element.text
                if not os.path.isfile(self.flowfile):
                    raise ValueError('Flow velocity file is missing or the given path is incorrect.')
            else:
                self.flowfile = None
            element = None
            fun = flow.find('function')
            if fun is not None:
                self.flowfunc = 0
                linear = fun.find('linear')
                if linear is not None:
                    element = linear.find('fmax')
                    if element is not None:
                        self.flowdepth = float(element.text)
                    element = linear.find('a')
                    if element is not None:
                        self.flowlina = float(element.text)
                    element = linear.find('b')
                    if element is not None:
                        self.flowlinb = float(element.text)
                    if self.flowlinb is not None and self.flowlina is None:
                         raise ValueError('Flow velocity linear function is declared but is missing a.')
                    if self.flowlina is not None and self.flowlinb is None:
                         raise ValueError('Flow velocity linear function is declared but is missing b.')
                edec = fun.find('expdecay')
                if edec is not None:
                    self.flowdecay = numpy.zeros((2,3), dtype=float)
                    data = defaultdict(list)
                    # Group into rows of (col, val) tuples
                    for val in edec.iter('fdvalue'):
                        data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                    # Sort columns and format into a space separated string
                    rows = []
                    for row in data:
                        rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                    # Build array from matrix string
                    self.flowdecay = numpy.array(numpy.mat(';'.join(rows)))
                if self.flowdecay is None and self.flowlina is None:
                    raise ValueError('Flow velocity function is declared but is missing some parameters.')
            else:
                self.flowfunc = None
        else:
            self.flowval = 0.
            self.flowfile = None

        # Extract sedimentation input information
        sed = None
        sed = root.find('sedinput')
        if sed is not None:
            self.sedOn = True
            element = None
            element = sed.find('val')
            if element is not None:
                self.sedval = float(element.text)
            else:
                self.sedval = 0.
            element = None
            element = sed.find('curve')
            if element is not None:
                self.sedfile = element.text
                if not os.path.isfile(self.sedfile):
                    raise ValueError('Sediment input file is missing or the given path is incorrect.')
            else:
                self.sedfile = None

            element = None
            fct = sed.find('function')
            if fct is not None:
                self.sedfunc = 0
                linear = fct.find('linear')
                if linear is not None:
                    element = linear.find('dmax')
                    if element is not None:
                        self.seddepth = float(element.text)
                    element = linear.find('a')
                    if element is not None:
                        self.sedlina = float(element.text)
                    element = linear.find('b')
                    if element is not None:
                        self.sedlinb = float(element.text)
                    if self.sedlinb is not None and self.sedlina is None:
                         raise ValueError('Flow velocity linear function is declared but is missing fa.')
                    if self.sedlina is not None and self.sedlinb is None:
                         raise ValueError('Flow velocity linear function is declared but is missing fb.')
                edec = fct.find('expdecay')
                if edec is not None:
                    self.seddecay = numpy.zeros((2,3), dtype=float)
                    data = defaultdict(list)
                    # Group into rows of (col, val) tuples
                    for val in edec.iter('sdvalue'):
                        data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                    # Sort columns and format into a space separated string
                    rows = []
                    for row in data:
                        rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                    # Build array from matrix string
                    self.seddecay = numpy.array(numpy.mat(';'.join(rows)))
                if self.seddecay is None and self.sedlina is None:
                    raise ValueError('Sediment input function is declared but is missing some parameters.')
            else:
                self.sedfunc = None
        else:
            self.sedval = 0.
            self.sedfile = None

        # Extract environmental function trapezoidal shape parameter
        envi = None
        envi = root.find('envishape')
        if envi is not None:
            element = None
            element = envi.find('depthshape')
            if element is not None:
                self.enviDepth = numpy.zeros((self.speciesNb,4), dtype=float)
                data = defaultdict(list)
                # Group into rows of (col, val) tuples
                for val in envi.iter('dvalue'):
                    data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                # Sort columns and format into a space separated string
                rows = []
                for row in data:
                    rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                # Build array from matrix string
                self.enviDepth = numpy.array(numpy.mat(';'.join(rows)))
            element = None
            element = envi.find('flowshape')
            if element is not None:
                self.enviFlow = numpy.zeros((self.speciesNb,4), dtype=float)
                data = defaultdict(list)
                # Group into rows of (col, val) tuples
                for val in envi.iter('fvalue'):
                    data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                # Sort columns and format into a space separated string
                rows = []
                for row in data:
                    rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                # Build array from matrix string
                self.enviFlow = numpy.array(numpy.mat(';'.join(rows)))
            element = None
            element = envi.find('sedshape')
            if element is not None:
                self.enviSed = numpy.zeros((self.speciesNb,4), dtype=float)
                data = defaultdict(list)
                # Group into rows of (col, val) tuples
                for val in envi.iter('svalue'):
                    data[int(val.attrib['row'])].append((int(val.attrib['col']), val.text))
                # Sort columns and format into a space separated string
                rows = []
                for row in data:
                    rows.append(' '.join([cols[1] for cols in sorted(data[row])]))
                # Build array from matrix string
                self.enviSed = numpy.array(numpy.mat(';'.join(rows)))

        # Get output directory
        out = None
        out = root.find('outfolder')
        if out is not None:
            self.outDir = out.text
        else:
            self.outDir = os.getcwd()+'/out'

        if self.makeUniqueOutputDir:
            if os.path.exists(self.outDir):
                self.outDir += '_'+str(len(glob.glob(self.outDir+str('*')))-1)

            os.makedirs(self.outDir)
            shutil.copy(self.inputfile,self.outDir)

        return
