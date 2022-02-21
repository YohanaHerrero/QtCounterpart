import numpy as np
import matplotlib.pyplot as plt

import os, sys
import copy
import pylab as p
from astropy.utils.data import get_pkg_data_filename

from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy import units as u
import scipy.ndimage as ndimage
from acstools import acszpt

import matplotlib
import math
from photutils import aperture_photometry
from photutils import CircularAperture
from decimal import Decimal

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
from PyQt5.QtWidgets import QRadioButton

# import own functions for reading and writing fits files and getting colours
import read_write_files as rwf
import colors as cl
from colormaps import viridis_simple as viridis #I could also use inferno_simple, magma_simple or plasma_simple instead

matplotlib.rcParams.update({'font.size': 20})

#GUI default: black background and white plot for pyqtgraph
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


################################################################################
######################## Change file names here ################################

#names of catalogues, info to create NB from MUSE data cube and directories
#This are default file names and directories

MUSE_DATA_CUBE = 'MUSE_cubes/DATACUBE_candels-cdfs-01_v2.0_dcsub_effnoised.fits' #MUSE cube

dict_cats = {
             'MUSE':'MW_60fields_main_table_v1.1.fits', # MUSE catalogue by default 
             #'MUSEnarrowband':'Field1/narrowband_MUSEWIDE.fits', #MUSE Narrowband, which is created in the code from the MUSE cube
             'UV':'CANDELS.GOODSS.F160W.v1_1.photom.cat', # photometric catalogue 
             'output':'test_output.fits', # output table where counterparts info is saved
             'HST1':'HST_images/acs_435w_candels-cdfs-01_cut_v2.0.fits', # HST images inside HST directory
             'HST2':'HST_images/acs_606w_candels-cdfs-01_cut_v2.0.fits',  
             'HST3':'HST_images/acs_814w_candels-cdfs-01_cut.fits',  
             'HST4':'HST_images/wfc3_160w_candels-cdfs-01_cut_v2.0.fits'   
             } 

width_window = [8.]  # width of the thumbnail cut-out windows in arcsec (for the main GUI)  
NBwindow_size = 8 #size of NB cut-out in arcsec 

######################## from here on, nothing else to change ##################
################################################################################




################################################################################
######################## Important functions for GUI ###########################


def read_hdu(infile,hdunum,nans_to_value=False,nanvalue=0,memmap=True):
    """
    image_data, header = read_hdu(infile,hdunum,nans_to_value=False,nanvalue=0)
    Read in a HDU <hdunum> from FITS-File <infile>, returns the
    data <image_data> and the header <header>.
    Additionally with nans_to_value = True (default: False) one can set
    all NaNs to a specific <nanvalue> (default: 0)
    This does not read data cubes!!
    """
    hdu = fits.open(infile,memmap=memmap)
    header = hdu[hdunum].header
    image_data = hdu[hdunum].data 
    hdu.close()
    if nans_to_value == True:
        # change NaNs to nanvalue
        select = np.isnan(image_data)
        image_data[select] = nanvalue
    return image_data,header

#Sets the cut/contrast of the small HST images
def vminvmax(image, scale=0.99, nans=0):
    """
    vmin,vmax = vminvmax(data, scale=0.99, nans=0)

    In:
    ---
    data - 2D array (intensity image)
    scale (=0.99) - scale factor for colorbar
    nans (=0) value with wich NaNs in the array are replaced

    Out:
    ---
    vmin, vmax -- values to be used in matplotlib.imshow as vmin & vmax
    """
    assert scale <= 1 and scale > 0 
    data = image.copy()    

    data[np.isnan(data)] = nans #Test element-wise for NaN and return result as a boolean array
    if len(data.shape) > 1:
        data = data.flatten() #get a copy of an given array collapsed into one dimension
    data = np.sort(data)

    length = len(data)

    vmax_index = int(np.floor(scale*length)) #np.floor([-1.7,-0.2,0.1,1.1,2])=-2,-1,0,1,2
    vmin_index = int(np.ceil((1-scale)*length)) #np.ceil([-1.7,-0.2,0.1,1.1,2])=-1,0,1,2,2

    vmin = data[vmin_index]
    vmax = data[vmax_index]

    return vmin, vmax


def pix_scale(header):
    """
    Prints out averaged pixel scale in arcseconds at reference pixel. 
    Ignores disstortions, but this OK considering our small field of view.

    In:
    ---
    header ... an astropy.io.fits header object

    Out:
    ---
    pix_scale ... average linear extent of a pixel in arc-seconds
    """

    wcs_obj = wcs.WCS(header,relax=False)
    scale_deg_xy = wcs.utils.proj_plane_pixel_scales(wcs_obj) #returns pixel scales along each axis of the image pixel at the CRPIX location once it is projected onto the “plane of intermediate world coordinates”
    scale_deg = p.sum(scale_deg_xy[:2]) / 2.

    scale_deg = Angle(scale_deg,unit=u.deg)
    scale_asec = scale_deg.arcsec

    return scale_asec


def radec_to_pix(ra_coords,dec_coords,header,origin=0):
    """
    Converts RA & DEC world coordinates to pixel coordinates.

    In:
    ---
    ra_coords ... 1D array of RA in degrees (float)
    dec_coords ... 1D array of corresponding DEC in degrees (float)
    header ... an astropy.io.fits header object
    origin ... output coordinates 0-indexed if 0, or 1-indexed if 1
               (default: 0)

    Out:
    ---
    x_coords, y_coords ... 2-tuple of 1D arrays with pixel coordinates
    """

    wcs_obj = wcs.WCS(header, relax=False)  # no HST only spec
                                            # allowed, only WCS
                                            # standard
    coords = wcs_obj.wcs_world2pix(ra_coords,dec_coords,origin)
    x_coords = coords[0]
    y_coords = coords[1]
    return x_coords, y_coords

def pix_to_radec(ra_coords,dec_coords,header,origin=0):
    """
    Converts pixel coordinates to RA & DEC world coordinates.
    """

    wcs_obj = wcs.WCS(header, relax=False)  # no HST only spec
                                            # allowed, only WCS
                                            # standard
    coords = wcs_obj.wcs_pix2world(ra_coords,dec_coords,origin)
    ra_coords = coords[0]
    dec_coords = coords[1]
    return ra_coords, dec_coords


#function to make the image of the small HST images clickable
#what we do here is to create a filter on the widget that we want that would detect
#a double click
def clickable(widget): 
    class Filter(QtCore.QObject):
        clicked = QtCore.pyqtSignal()
        def eventFilter(self, obj, event):
            if obj == widget:
                if event.type() == QtCore.QEvent.MouseButtonDblClick:
                    if obj.rect().contains(event.pos()):
                        self.clicked.emit()
                        return True
            return False
    filter = Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked  

    
################################################################################
    
### Read MUSE and photometric catalogues

try:
    print('Reading MUSE catalogue:',dict_cats['MUSE'])
    info_MUSE,info_MUSE_ids = rwf.read_fits_table2dict(dict_cats['MUSE'],1)
    print('Successfully')
except:
    print('Could not read input catalogue:',dict_cats['MUSE'],' please check and run again.')
    exit()

ids_all = np.asarray(sorted(np.asarray(list(info_MUSE_ids.keys()))),dtype=str)
#Dictionary in Python is an unordered collection of data values, used to store data values like a map
#keys() method in Python Dictionary, returns a view object that displays a list of all the keys in the dictionary    

'''
#Read MUSE narrowband (not needed because I create the NB in the code)
try:
    print('Reading MUSE narrowband:',dict_cats['MUSEnarrowband'])
    print('Successfully')
except:
    print('Could not read MUSE narrowband:',dict_cats['MUSEnarrowband'],' please check and run again.')
    exit()
'''

try:
    print('Reading photometric catalogue (ID, RA, Dec):',dict_cats['UV'])
    UV_open     = open(dict_cats['UV'],'r')
    UV_data     = np.genfromtxt(UV_open) #Load data from a text file
    UV_ID       = [row[0] for row in UV_data] # ID of Object in cataloge
    UV_RA       = [row[2] for row in UV_data]
    UV_DEC      = [row[3] for row in UV_data]
    UV_open.close()
    info_UV= {'ID':np.asarray(UV_ID),
              'RA':np.asarray(UV_RA),
              'DEC':np.asarray(UV_DEC)}
    print('Successfully')
except:
    print('Could not read input catalogue:',dict_cats['UV'],' please check and run again.')
    exit()
    

#Read HST images
try:
    print('Reading HST band 435 image:',dict_cats['HST1'])
    print('Successfully')
except:
    print('Could not read band 435:',dict_cats['HST1'],' please check and run again.')
    exit()
    
try:
    print('Reading HST band 606 image:',dict_cats['HST2'])
    print('Successfully')
except:
    print('Could not read band 606:',dict_cats['HST2'],' please check and run again.')
    exit()

try:
    print('Reading HST band 775/814 image:',dict_cats['HST3'])
    print('Successfully')
except:
    print('Could not read band 775/814:',dict_cats['HST3'],' please check and run again.')
    exit()

try:
    print('Reading HST band 160 image:',dict_cats['HST4'])
    print('Successfully')
except:
    print('Could not read band 160:',dict_cats['HST4'],' please check and run again.')
    exit()


################################################################################
######################## Main GUI to identify counterparts #####################

counterparts = {} # initialise dictionary with counterparts (empty)

class main_GUI(QtGui.QWidget):
    
    def __init__(self, *args, **kwargs):
        super(main_GUI, self).__init__(*args, **kwargs)
        self.numDegrees436 = 0 
        self.numDegrees606 = 0
        self.numDegrees775 = 0
        self.numDegrees160 = 0
        self.numDegreesMuse = 0
        self.numDegreesHST = 0
        self.initUI() #creates window, the main GUI one
        
    def initUI(self):

        self.setGeometry(0, 80, 600, 500) #setGeometry(left start point on screen, top start point, width, height) in px
        self.setWindowTitle('Identifying Counterparts') 
        
        self.field = 0 # initialise MUSE field (I pass a MUSE data file (i.e. 60 fields file) with gal in different fields)
        self.j = 0 # initialise counter for the whole GUI
        self.get_infos() #I call the funcion (defined below). It passes the info from the MUSE file, ID, SN, z, lead line, etc
                
        self.grid = QtGui.QGridLayout()
        self.grid.setSpacing(10)
        
        #editable scale cuts for the small HST images with their default values
        self.le_resolution436 = QtGui.QLineEdit('0.99',self) #default scale, 0.99
        self.le_resolution436.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.le_resolution436,4,8,1,1)
        
        self.le_resolution606 = QtGui.QLineEdit('0.99',self)
        self.le_resolution606.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.le_resolution606,4,10,1,1)
        
        self.le_resolution775 = QtGui.QLineEdit('0.99',self)
        self.le_resolution775.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.le_resolution775,4,13,1,1)
        
        self.le_resolution160 = QtGui.QLineEdit('0.99',self)
        self.le_resolution160.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.le_resolution160,4,16,1,1)
        
        self.fill_window() #I call the function
        self.setLayout(self.grid)
        self.show()
        
        
    
    #functions that make the small hst images the main one after doing double click on their label
    def make_436smallHST_mainHST(self): 
        self.hst_data = self.hst_data_436
        self.hst_header = self.hst_header_436
        self.hst_level_max =self.hst_level_max436
        self.hst_level_min = self.hst_level_min436
        self.main_image_Widget()
        self.fill_window()     

    def make_606smallHST_mainHST(self): 
        self.hst_data = self.hst_data_606
        self.hst_header = self.hst_header_606
        self.hst_level_max =self.hst_level_max606
        self.hst_level_min = self.hst_level_min606
        self.main_image_Widget()
        self.fill_window()
 
    def make_775smallHST_mainHST(self): 
        self.hst_data = self.hst_data_775
        self.hst_header = self.hst_header_775
        self.hst_level_max =self.hst_level_max775
        self.hst_level_min = self.hst_level_min775
        self.main_image_Widget()
        self.fill_window()

    def make_160smallHST_mainHST(self): 
        self.hst_data = self.hst_data_160
        self.hst_header = self.hst_header_160
        self.hst_level_max =self.hst_level_max160
        self.hst_level_min = self.hst_level_min160
        self.main_image_Widget()
        self.fill_window()  
        
        
    #small HST images setting style   
    def image_Widget436(self,label,pos,data,header): #I pass the parameters defined in the function fill_window 
        here_label = QtGui.QLabel(label, self)        
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,0,pos,1,1) # row, column, spanning rows, spanning columns
        self.img436 = pg.GraphicsLayoutWidget(self)
        self.img436.setMinimumSize(180,180)
        self.img436.setMaximumSize(180,180)
        clickable(self.img436).connect(self.make_436smallHST_mainHST) #when we click twice on the image of the small hst image we call the function    
        self.grid.addWidget(self.img436,1,pos,3,3) 
        
        
        # cut HST image to relevant area
        hst_ps = pix_scale(header)  
        self.width_window_pix_hst = width_window[0] / hst_ps  
        x_hst,y_hst = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,header)  
        hst_data_cut =self.hst_data_436[int(y_hst-self.width_window_pix_hst):
                                    int(y_hst+self.width_window_pix_hst),
                                    int(x_hst-self.width_window_pix_hst):
                                    int(x_hst+self.width_window_pix_hst)]
        hst_data_cut=ndimage.gaussian_filter(hst_data_cut, sigma=.6, order=0) #smooth image with sigma=0.6
        self.img_hst436 = pg.ImageItem(border='k')
        self.img_hst436.setLookupTable(viridis)
        self.img_hst436.setImage(np.rot90(hst_data_cut)[::-1])
        self.img_hst436.setLevels([self.hst_level_min436, self.hst_level_max436]) 
        self.view_hst436 = self.img436.addViewBox()
        self.view_hst436.addItem(self.img_hst436)
        
        # zooming in a little bit in the object
        # giving the middle of the area and the width to be shown in x and y
        self.view_hst436.setRange(QtCore.QRectF(self.width_window_pix_hst/2., #left, top, width and height
                                             self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst,
                                             self.width_window_pix_hst))

        self.view_hst436.setAspectLocked(True)
        
            
    def image_Widget606(self,label,pos,data,header): #I pass the parameters defined in the function fill_window (2 functions below)
        here_label = QtGui.QLabel(self.tr(label), self)
        
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,0,pos,1,1) # row, column, spanning rows, spanning columns
        self.img606 = pg.GraphicsLayoutWidget(self)
        self.img606.setMinimumSize(180,180)
        self.img606.setMaximumSize(180,180)
        clickable(self.img606).connect(self.make_606smallHST_mainHST) #when we click twice on the image of the small hst image we call the function    
        self.grid.addWidget(self.img606,1,pos,3,3) 
        
        # cut HST image to relevant area
        hst_ps = pix_scale(header)
        self.width_window_pix_hst = width_window[0] / hst_ps
        x_hst,y_hst = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,header)
        
        hst_data_cut =self.hst_data_606[int(y_hst-self.width_window_pix_hst):
                                    int(y_hst+self.width_window_pix_hst),
                                    int(x_hst-self.width_window_pix_hst):
                                    int(x_hst+self.width_window_pix_hst)]
        hst_data_cut=ndimage.gaussian_filter(hst_data_cut, sigma=.6, order=0)
        img_hst606 = pg.ImageItem(border='k')
        img_hst606.setLookupTable(viridis)
        img_hst606.setImage(np.rot90(hst_data_cut)[::-1])
        img_hst606.setLevels([self.hst_level_min606, self.hst_level_max606]) 
        self.view_hst606 = self.img606.addViewBox()
        self.view_hst606.addItem(img_hst606)
        
        # zooming in a little bit in the object
        # giving the middle of the area and the width to be shown in x and y
        self.view_hst606.setRange(QtCore.QRectF(self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst,
                                             self.width_window_pix_hst))
        
        self.view_hst606.setAspectLocked(True)
        
    def image_Widget775(self,label,pos,data,header): 
        here_label = QtGui.QLabel(label, self)
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,0,pos,1,1) # row, column, spanning rows, spanning columns
        self.img775 = pg.GraphicsLayoutWidget(self)
        self.img775.setMinimumSize(180,180)
        self.img775.setMaximumSize(180,180)
        clickable(self.img775).connect(self.make_775smallHST_mainHST) #when we click twice on the small hst label we call the function    
        self.grid.addWidget(self.img775,1,pos,3,3) 
              
        # cut HST image to relevant area, otherwise GUI is too slow
        hst_ps = pix_scale(header)
        self.width_window_pix_hst = width_window[0] / hst_ps
        x_hst,y_hst = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,header)
        
        hst_data_cut =self.hst_data_775[int(y_hst-self.width_window_pix_hst):
                                    int(y_hst+self.width_window_pix_hst),
                                    int(x_hst-self.width_window_pix_hst):
                                    int(x_hst+self.width_window_pix_hst)]
        hst_data_cut=ndimage.gaussian_filter(hst_data_cut, sigma=.6, order=0)
        img_hst775 = pg.ImageItem(border='k')
        img_hst775.setLookupTable(viridis)
        img_hst775.setImage(np.rot90(hst_data_cut)[::-1])
        img_hst775.setLevels([self.hst_level_min775, self.hst_level_max775]) 
        self.view_hst775 = self.img775.addViewBox()
        self.view_hst775.addItem(img_hst775)
        
        # zooming in a little bit in the object
        # giving the middle of the area and the width to be shown in x and y
        self.view_hst775.setRange(QtCore.QRectF(self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst,
                                             self.width_window_pix_hst))
        
        self.view_hst775.setAspectLocked(True)

    def image_Widget160(self,label,pos,data,header):
        here_label = QtGui.QLabel(label, self)
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,0,pos,1,1) # row, column, spanning rows, spanning columns
        self.img160 = pg.GraphicsLayoutWidget(self)
        self.img160.setMinimumSize(180,180)
        self.img160.setMaximumSize(180,180)
        clickable(self.img160).connect(self.make_160smallHST_mainHST) #when we click twice on the image of the small hst image we call the function    
        self.grid.addWidget(self.img160,1,pos,3,3) 
        
        # cut HST image to relevant area, otherwise GUI is too slow
        hst_ps = pix_scale(header)  
        self.width_window_pix_hst160 = width_window[0]*2 / hst_ps #*2 to have the same size as the rest
        
        x_hst,y_hst = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,header)
        
        hst_data_cut =self.hst_data_160[int(y_hst-self.width_window_pix_hst160):
                                    int(y_hst+self.width_window_pix_hst160),
                                    int(x_hst-self.width_window_pix_hst160):
                                    int(x_hst+self.width_window_pix_hst160)]    
        hst_data_cut=ndimage.gaussian_filter(hst_data_cut, sigma=.6, order=0)
        img_hst160 = pg.ImageItem(border='k')
        img_hst160.setLookupTable(viridis)
        img_hst160.setImage(np.rot90(hst_data_cut)[::-1])
        img_hst160.setLevels([self.hst_level_min160, self.hst_level_max160]) 
        self.view_hst160 = self.img160.addViewBox()
        self.view_hst160.addItem(img_hst160)
        
        # zooming in a little bit in the object
        # giving the middle of the area and the width to be shown in x and y
        self.view_hst160.setRange(QtCore.QRectF(self.width_window_pix_hst160/2.,
                                             self.width_window_pix_hst160/2.,
                                             self.width_window_pix_hst160,
                                             self.width_window_pix_hst160))
        
        self.view_hst160.setAspectLocked(True)
        
                
    #MUSE image setting style 
    def image_Widget_narrowband_muse(self,label,pos,data,header): 
        here_label = QtGui.QLabel(label, self)
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,13,pos,1,1) # row, column, spanning rows, spanning columns
        self.imgM = pg.GraphicsLayoutWidget(self)
        self.imgM.setMinimumSize(360,360)
        self.imgM.setMaximumSize(360,360)
        self.grid.addWidget(self.imgM,14,pos,8,3) 
        
        # cut muse image to relevant area
        muse_ps = pix_scale(header)
        self.width_window_pix_muse = width_window[0] / muse_ps
        x_muse,y_muse = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,header) #ra, dec of gal in px
        muse_data_cut =self.muse_data #we dont need to cut muse data, because we have already cut a NB of NBwindow_size
        muse_data_cut=ndimage.gaussian_filter(muse_data_cut, sigma=.6, order=0)
              
        img_muse = pg.ImageItem(border='k')
        img_muse.setLookupTable(viridis)
        img_muse.setImage(np.rot90(muse_data_cut)[::-1])       
        img_muse.setLevels([self.muse_level_min, self.muse_level_max]) 
        self.view_muse = self.imgM.addViewBox()
        self.view_muse.addItem(img_muse)
        
        # zooming in a little bit in the object
        # giving the middle of the area and the width to be shown in x and y
        self.view_muse.setRange(QtCore.QRectF(self.width_window_pix_muse/200., #left, top, width and height
                                             self.width_window_pix_muse/200.,
                                             self.width_window_pix_muse,
                                             self.width_window_pix_muse))
        #CHANGE THE ZOOM FROM MUSE WITH THIS ABOVE
        self.view_muse.setAspectLocked(True)
      
    #Main HST image setting style, by default we choose the 436 band to be the main image   
    def main_image_Widget(self):
        # HST image with highest depth and clickable objects
        here_label = QtGui.QLabel('HST main image', self)
        here_label.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(here_label,0,0,1,1)
        self.img_hst = pg.GraphicsLayoutWidget(self)
        # make main HST image bigger
        self.img_hst.setMinimumSize(360,360)
        self.img_hst.setMaximumSize(360,360)
        self.grid.addWidget(self.img_hst,1,0,12,12) 
        
        self.x_hst,self.y_hst = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,self.hst_header)

        img_hst = pg.ImageItem(border='k')
        img_hst.setLookupTable(viridis)
        img_hst.setImage(np.rot90(ndimage.gaussian_filter(self.hst_data, sigma=.6, order=0))[::-1]) #change to self.hst_data_606 for example if I wanna display another band by default
        img_hst.setLevels([self.hst_level_min, self.hst_level_max]) #change to self.hst_level_max606 and self.hst_level_min606if I wanna display another band by default
        self.view_hst_main = self.img_hst.addViewBox()
        self.view_hst_main.addItem(img_hst)
        
        hst_ps = pix_scale(self.hst_header) #change to self.hst_header_606 if I wanna display another band by default
        self.width_window_pix_hst = width_window[0] / hst_ps
        
        self.view_hst_main.setRange(QtCore.QRectF(self.x_hst-self.width_window_pix_hst/2.,
                                             self.y_hst-self.width_window_pix_hst/2.,
                                             self.width_window_pix_hst,
                                             self.width_window_pix_hst))

        self.view_hst_main.setAspectLocked(True)
        
        #the function allows to click anywhere on the main hst image to obtain the position of the clicked place
        def clickToObtainPositionOnImage(event):
            position = event.pos()
            self.x_WithoutCounterpart = position[0]
            self.y_WithoutCounterpart = position[1]
            self.RA_WithoutCounterpart,self.Dec_WithoutCounterpart=pix_to_radec(self.x_WithoutCounterpart,self.y_WithoutCounterpart,self.hst_header,origin=0)
            return self.RA_WithoutCounterpart,self.Dec_WithoutCounterpart, self.x_WithoutCounterpart, self.y_WithoutCounterpart
        
        self.RA_WithoutCounterpart='-' #default if the hst img was not clicked, this is saved in the output catalog
        self.Dec_WithoutCounterpart='-'
        self.forced_photometry='-'        
        img_hst.mouseClickEvent = clickToObtainPositionOnImage #however, if the hst img is clicked we take the positions of the mouse click
      
    #we create all labels and buttons that will accompany the previous images    
    def fill_window(self):
     
        # create widgets/images with HST bands (436, 606, 775/816 and 160)
        self.image_Widget_narrowband_muse('MUSE narrowband',0,self.muse_data,self.header_muse)
        
        self.image_Widget436('HST band 435',7,self.hst_data_436,self.hst_header_436)
        self.image_Widget606(self.tr('HST band 606'),10,self.hst_data_606,self.hst_header_606)
        self.image_Widget775('HST band 775/814',13,self.hst_data_775,self.hst_header_775)
        self.image_Widget160('HST band 160',16,self.hst_data_160,self.hst_header_160)
        
        self.resolution_bt =  QtGui.QPushButton('Cuts:', self)
        self.resolution_bt.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.resolution_bt,4,7,1,1)
        
        self.resolution_bt.clicked.connect(self.small_hst_resolution)
               
        if self.le_resolution436.text() != '0.99':
            self.le_resolution436.setText(self.le_resolution436.text())
        
        #Button to reset the hst small images, we go back to the default cut
        self.btn_reset_hst = QtGui.QPushButton('Reset images', self)
        # trigger function reset_hst_images when clicking
        self.btn_reset_hst.clicked.connect(self.reset_hst_images)           
        self.btn_reset_hst.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.btn_reset_hst,4,17,1,1)
        
        #Button to clear the selected counterparts, all goes back to the default window without numbers
        self.btn_clear =  QtGui.QPushButton('Clear counterparts', self)
        self.btn_clear.clicked.connect(self.clear_counterparts) # trigger function clear_counterparts when clicking
        self.btn_clear.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.btn_clear,9,14,1,2)
        
        #confidence bullets to quantify how sure we are in the selection of the counterpart
        #0 means no sure at all and 3 means completely sure
        confidence_label =  QtGui.QLabel('Confidence:', self)
        confidence_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(confidence_label,7,11,1,1)
        #create the bullets of confidence      
        self.b1 = QRadioButton("0",self)
        self.b1.setChecked(True) #default marked button, confidence 0
        self.le_confidence = '0' #default to save in the output file if no confidence is clicked
        self.b1.toggled.connect(lambda:self.btnstate(self.b1)) #The lambda is returned from a function to prevent that modified values will be used later
        self.b1.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.b1,7,12,1,1)
		
        self.b2 = QRadioButton("1",self)
        self.b2.toggled.connect(lambda:self.btnstate(self.b2))
        self.b2.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.b2,7,13,1,1)
      
        self.b3 = QRadioButton("2",self)
        self.b3.toggled.connect(lambda:self.btnstate(self.b3))
        self.b3.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.b3,7,14,1,1)
      
        self.b4 = QRadioButton("3",self)
        self.b4.toggled.connect(lambda:self.btnstate(self.b4))
        self.b4.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.b4,7,15,1,1)
        
        
        #create a clickable square on the left of the labels that indicate the situation
        nomatch_checkbox =  QtGui.QCheckBox('No match to MUSE in photometry', self)
        self.le_nomatch = '-' #default to save in the output file if the checkbox is not clicked        
        nomatch_checkbox.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(nomatch_checkbox,8,11,1,2) 
        nomatch_checkbox.stateChanged.connect(self.nomatch_button)
        
        nocounter_checkbox =  QtGui.QCheckBox('Counterpart not in catalog', self)
        self.le_nocatalog = '-' #default to save in the output file if the checkbox is not clicked 
        nocounter_checkbox.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(nocounter_checkbox,9,11,1,2)
        nocounter_checkbox.stateChanged.connect(self.nocounter_button)
        
        more_checkbox =  QtGui.QCheckBox('More than one counterpart', self)
        self.le_more = '-' #default to save in the output file if the checkbox is not clicked 
        more_checkbox.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(more_checkbox,10,11,1,2)
        more_checkbox.stateChanged.connect(self.more_button)      
      
        # create the main HST image to click the possible counterparts
        self.main_image_Widget()
        # Get the photometric information
        self.get_pos_UV()  
        
        # display information on MUSE data on window (taken from the MUSE file that we pass on the first window of the GUI)
        MUSE_info = QtGui.QLabel('MUSE information', self)
        MUSE_info.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(MUSE_info,13,7,1,1)
        MUSE_ID_label = QtGui.QLabel('MUSE ID:\t'+str(self.MUSE_ID), self)
        MUSE_ID_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_ID_label,14,7,1,1)
        
        MUSE_RA_label = QtGui.QLabel('RA:\t'+str(np.round(self.MUSE_RA,3)), self)
        MUSE_RA_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_RA_label,15,7,1,1)
        MUSE_DEC_label = QtGui.QLabel('DEC:\t'+str(np.round(self.MUSE_Dec,3)), self)
        MUSE_DEC_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_DEC_label,16,7,1,1)
        MUSE_z_label = QtGui.QLabel('z:\t'+str(np.round(self.MUSE_Z,3)), self)
        MUSE_z_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_z_label,17,7,1,1)
        
        MUSE_confidence_label = QtGui.QLabel('Confidence:\t'+str(self.MUSE_confidence), self) #strongest line in the muse image
        MUSE_confidence_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_confidence_label,18,7,1,1)
        
        MUSE_line_label = QtGui.QLabel('Lead line:\t'+str(self.MUSE_lead_line), self)  
        MUSE_line_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_line_label,19,7,1,1)
        MUSE_SN_label = QtGui.QLabel('S/N:\t'+str(np.round(self.MUSE_SN,2)), self)
        MUSE_SN_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(MUSE_SN_label,20,7,1,1)
        

        # information on UV counterpart
        # blank by default until we click on any counterpart
        UV_info = QtGui.QLabel('Photometric information', self)
        UV_info.setFont(QtGui.QFont('SansSerif', 11, QtGui.QFont.Bold))
        self.grid.addWidget(UV_info,6,7,1,1)
        UV_ID_label =  QtGui.QLabel('UV ID:', self)
        UV_ID_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(UV_ID_label,7,7,1,1)
        self.le_UV_ID = QtGui.QLineEdit(self)
        self.le_UV_ID.setFont(QtGui.QFont('SansSerif', 11))
        self.le_UV_ID.setText('0')
        self.grid.addWidget(self.le_UV_ID,7,8,1,1)
        
        UV_RA_label =  QtGui.QLabel('<font>&Delta;</font> RA', self)
        UV_RA_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(UV_RA_label,8,7,1,1)
        self.le_UV_RA = QtGui.QLineEdit(self)
        self.le_UV_RA.setFont(QtGui.QFont('SansSerif', 11))
        self.le_UV_RA.setText('0')
        self.grid.addWidget(self.le_UV_RA,8,8,1,1)
        
        UV_DEC_label =  QtGui.QLabel('<font>&Delta;</font> Dec', self)
        UV_DEC_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(UV_DEC_label,9,7,1,1)
        self.le_UV_DEC = QtGui.QLineEdit(self)
        self.le_UV_DEC.setFont(QtGui.QFont('SansSerif', 11))
        self.le_UV_DEC.setText('0')
        self.grid.addWidget(self.le_UV_DEC,9,8,1,1)
        
        UV_z_label =  QtGui.QLabel('photo-z:', self)
        UV_z_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(UV_z_label,10,7,1,1)
        self.le_UV_z = QtGui.QLineEdit(self)
        self.le_UV_z.setFont(QtGui.QFont('SansSerif', 11))
        self.le_UV_z.setText('0')
        self.grid.addWidget(self.le_UV_z,10,8,1,1)
        
        UV_z_confidence_label =  QtGui.QLabel('z-confidence:', self)
        UV_z_confidence_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(UV_z_confidence_label,11,7,1,1)
        self.le_UV_z_confidence = QtGui.QLineEdit(self)
        self.le_UV_z_confidence.setFont(QtGui.QFont('SansSerif', 11))
        self.le_UV_z_confidence.setText('0')
        self.grid.addWidget(self.le_UV_z_confidence,11,8,1,1)
        
        # line edit for comments
        comments_label =  QtGui.QLabel('Comments:', self)
        comments_label.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(comments_label,11,11,1,1)
        self.le_comment = QtGui.QLineEdit(self)
        self.le_comment.setFont(QtGui.QFont('SansSerif', 11))
        self.le_comment.setText('')
        self.grid.addWidget(self.le_comment,11,12,1,5)
        
        # line edit to jump to another MUSE ID
        self.jump_bt =  QtGui.QPushButton('Jump to MUSE ID:', self)
        self.jump_bt.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.jump_bt,17,9,1,1)
        self.le_jump = QtGui.QLineEdit(self)
        self.le_jump.setFont(QtGui.QFont('SansSerif', 11))
        self.le_jump.setText('')
        self.grid.addWidget(self.le_jump,17,10,1,1)
        self.jump_bt.clicked.connect(self.find_object)
        
        # line edit to jump to the next unclassified object
        self.unclassi_bt =  QtGui.QPushButton('Jump to next unclassified ID', self)
        self.unclassi_bt.setFont(QtGui.QFont('SansSerif', 11))
        self.grid.addWidget(self.unclassi_bt,18,9,1,2)
        self.unclassi_bt.clicked.connect(self.unclassified_ID)
        
        # buttons to move forward or back
        self.btn_prev = QtGui.QPushButton('previous', self)
        self.btn_prev.setFont(QtGui.QFont('SansSerif', 11))
        self.btn_prev.clicked.connect(self.show_next) #we call the function
        self.grid.addWidget(self.btn_prev,20,15,1,1)
        
        self.btn_next = QtGui.QPushButton('next', self)
        self.btn_next.setFont(QtGui.QFont('SansSerif', 11))
        self.btn_next.clicked.connect(self.show_next) #we call the function
        self.grid.addWidget(self.btn_next,20,16,1,1)
    
    #function that deselect all counterparts and goes back to the default window, with the info from counterparts in blank
    def clear_counterparts(self):
        self.fill_window()
        
        
    #function that changes the cuts of the small HST images
    def small_hst_resolution(self ):
        if self.le_resolution436.text() != '0.99':
             print('printing', self.le_resolution436.text())
             self.hst_level_min436, self.hst_level_max436 = vminvmax(self.hst_data_436, scale=float(self.le_resolution436.text()), nans=0)
             
        if self.le_resolution606.text() != '0.99':
             self.hst_level_min606, self.hst_level_max606 = vminvmax(self.hst_data_606, scale=float(self.le_resolution606.text()), nans=0)
          
        if self.le_resolution775.text() != '0.99':
             self.hst_level_min775, self.hst_level_max775 = vminvmax(self.hst_data_775, scale=float(self.le_resolution775.text()), nans=0)
          
        if self.le_resolution160.text() != '0.99':
            self.hst_level_min160, self.hst_level_max160 = vminvmax(self.hst_data_160, scale=float(self.le_resolution160.text()), nans=0)
        
        self.fill_window()


   #function to reset the cuts in the small HST images, back to default scale=0.99
    def reset_hst_images(self,combo_resol):
          self.get_HST_images()
          self.hst_level_min436, self.hst_level_max436 = vminvmax(self.hst_data_436, scale=0.99, nans=0)
          self.hst_level_min606, self.hst_level_max606 = vminvmax(self.hst_data_606, scale=0.99, nans=0)
          self.hst_level_min775, self.hst_level_max775 = vminvmax(self.hst_data_775, scale=0.99, nans=0)
          self.hst_level_min160, self.hst_level_max160 = vminvmax(self.hst_data_160, scale=0.99, nans=0)
          self.le_resolution436.setText('0.99')
          self.fill_window()

    
    #function that quantifies which confidence circular button was pressed in order to save it in the output file
    def btnstate(self,b):
      if b.text() == "0":
         if b.isChecked() == True:
            self.le_confidence = b.text()
		
      if b.text() == "1":
         if b.isChecked() == True:
            self.le_confidence = b.text()

      if b.text() == "2":
         if b.isChecked() == True:
            self.le_confidence = b.text()

      if b.text() == "3":
         if b.isChecked() == True:
            self.le_confidence = b.text()
        
        
        
    #function that knows if the "no match to MUSE in photometry" checkbox was clicked in order to save it in the output file
    def nomatch_button(self,state):
        if state == QtCore.Qt.Checked:  
            self.le_nomatch = 'Yes'
        else: 
            self.le_nomatch = '-'   
       
        
        
    #function that knows if the "counterpart not in catalog" checkbox was clicked in order to save it in the output file    
    def nocounter_button(self,state):	
      if state == QtCore.Qt.Checked:
            self.le_nocatalog = 'Yes'
      else: 
            self.le_nocatalog = '-'



    #function that knows if the "more than one counterpart" checkbox was clicked in order to save it in the output file
    def more_button(self,state):	
      if state == QtCore.Qt.Checked:
            self.le_more = 'Yes'
      else: 
            self.le_more = '-'



    #function to jump to the MUSE ID that we specify   
    def find_object(self):
        new_object = self.le_jump.text()
        for i, item in enumerate(ids_all): #there are 2381 objects here, the number of gal in MW_60fields_main_table_v1.1.fits
            if item == new_object:
                # clear all Widgets for next object    
                for k in reversed(range(self.grid.count())): 
                    self.grid.itemAt(k).widget().setParent(None)
    
                self.j=i # move forward until the position of the new ID. I change the counter to where I am
                self.get_infos() # get info on new ID
                self.fill_window() # show main window again
                print('position',self.j)
                break

    #function to jump to the next unclassified ID  
    def unclassified_ID(self):
        output_filename = get_pkg_data_filename(dict_cats['output']) 
        hdul = fits.open(output_filename)
        output_data = hdul[1].data  
        hdul.close()
        for i, museid in enumerate(np.asarray(ids_all,dtype=str)): #I treat ids_all as int to be consistent with the output_data type
            if museid not in output_data['ID']:
                # clear all Widgets for next object    
                for k in reversed(range(self.grid.count())): 
                    self.grid.itemAt(k).widget().setParent(None)
           
                self.j=i # move forward until the position of the new ID. I change the counter to where I am
                print('position3',self.j)
                self.get_infos() # get info on new ID
                self.fill_window() # show main window again
                
                break #once I find the unclassified object I leave the for loop
                


    def get_HST_images(self): #list of filenames of the HST directory and MUSE
            
        #this should give me all images in a list if the have area-field in this writting format in the name
        #hst_images = [i for i in np.asarray(os.listdir('HST_images/')) if self.area+'-'+self.field in i]
        
        #the list 'hst_images' only gives you one image, because in the loop you ask if 'band' is in the 
        #file name (the file name in the loop is 'i'), so it will only give you the HST image for the 435 band in this case.
        band = '435'
        hst_image436 = [i for i in np.asarray(os.listdir('HST_images/')) 
                     if self.area+'-'+self.field in i and band in i]   # area is cdfs or cosmos     
        band = '606'
        hst_image606 = [i for i in np.asarray(os.listdir('HST_images/')) 
                      if self.area+'-'+self.field in i and band in i]
        band = '816'
        hst_image816 = [i for i in np.asarray(os.listdir('HST_images/')) 
                      if self.area+'-'+self.field in i and band in i]        
        band = '160'
        hst_image160 = [i for i in np.asarray(os.listdir('HST_images/')) 
                      if self.area+'-'+self.field in i and band in i]

        self.hst_data, self.hst_header = read_hdu(dict_cats['HST1'], 0, nans_to_value=True) #I set by default the hst data of the big hst image and muse narrow band as the 436 hst band
        
        self.hst_data_436, self.hst_header_436 = read_hdu(dict_cats['HST1'], 0, nans_to_value=True)
        self.hst_data_606, self.hst_header_606 = read_hdu(dict_cats['HST2'], 0, nans_to_value=True)
        self.hst_data_775, self.hst_header_775 = read_hdu(dict_cats['HST3'], 0, nans_to_value=True)
        self.hst_data_160, self.hst_header_160 = read_hdu(dict_cats['HST4'], 0, nans_to_value=True)

        self.hst_level_min, self.hst_level_max = vminvmax(self.hst_data) #set the cuts of the big images
        self.hst_level_min436, self.hst_level_max436 = vminvmax(self.hst_data_436) #set the cuts of the small hst images
        self.hst_level_min606, self.hst_level_max606 = vminvmax(self.hst_data_606)
        self.hst_level_min775, self.hst_level_max775 = vminvmax(self.hst_data_775)
        self.hst_level_min160, self.hst_level_max160 = vminvmax(self.hst_data_160)
        self.muse_level_min, self.muse_level_max = vminvmax(self.muse_data)
        
    #get information from the MUSE catalog to be displayed with the fill_window function    
    def get_infos(self):
        
        self.id_here = ids_all[self.j]   
      
        # read parameters from MUSE catalogue for current object
        self.MUSE_ID = info_MUSE_ids[self.id_here]['UNIQUE_ID']
        self.MUSE_Z = info_MUSE_ids[self.id_here]['Z']
        self.MUSE_confidence = info_MUSE_ids[self.id_here]['CONFIDENCE']
        self.MUSE_lead_line = info_MUSE_ids[self.id_here]['LEAD_LINE']
        self.MUSE_SN = info_MUSE_ids[self.id_here]['SN']
        self.MUSE_other_lines = info_MUSE_ids[self.id_here]['OTHER_LINES']
        self.MUSE_RA = info_MUSE_ids[self.id_here]['RA']
        self.MUSE_Dec = info_MUSE_ids[self.id_here]['DEC']
     
        #lambda, ra and dec of interest to extract the NB from the MUSE cube
        ra = self.MUSE_RA
        dec = self.MUSE_Dec
        z = self.MUSE_Z
        line = self.MUSE_lead_line
      
        if line == 'Lya':
            lam=(1+z)*1215.67 #air wavelength
        elif line == 'O3_2':
            lam=(1+z)*4958.91
        elif line == 'O2':
            lam=(1+z)*3727.09
        elif line == 'O3_1':
            lam=(1+z)*5006.84
        elif line == 'Ha':
            lam=(1+z)*6562.8
        elif line == 'Hb':
            lam=(1+z)*4861.36
        elif line == 'N2_2':
            lam=(1+z)*6548.05
        elif line == 'Ne3':
            lam=(1+z)*3869
        elif line == 'C4':
            lam=(1+z)*1549.48
        elif line == 'He2':
            lam=(1+z)*1640.4
        elif line == 'C3':
            lam=(1+z)*1908.73
        else:
            print('The line',line,'is not listed')
        lam1 = lam-400  
        lam2 = lam+400  
    
        #Open cube, get data and header
        cube_open = fits.open(MUSE_DATA_CUBE)
        cube_data = cube_open[1].data
        cube_header = cube_open[1].header 
        w = WCS(cube_header) 
        cube_open.close()    
      
        # read pixel size and wavelength range
        try:
            pix = cube_header['CD2_2'] #px size
        except:
            pix = cube_header['CDELT2']
        lam_bin = cube_header['CD3_3']
        lam_start = cube_header['CRVAL3']
    
        # window size in pixels
        if NBwindow_size != 'full':
            window_size_pix = NBwindow_size/(pix*60.*60.)
            window_size_pix_half = int(window_size_pix/2.)
    
        # position in pixel coordinates
        x,y,z = w.wcs_world2pix(ra,dec,lam1,0)   
    
        # wavelength in pixel
        lam1_pix = int((lam1-lam_start)/lam_bin)
        lam2_pix = int((lam2-lam_start)/lam_bin)
    
        # correct for small shift due to integer x,y positions
        ra_middle,dec_middle,lam_middle = w.wcs_pix2world(int(x),int(y),lam1,0)
    
        if NBwindow_size != 'full':
            x1,x2 = int(x-window_size_pix_half)+1,int(x+window_size_pix_half)+1
            y1,y2 = int(y-window_size_pix_half)+1,int(y+window_size_pix_half)+1 
            cutout = np.asarray(cube_data[lam1_pix:lam2_pix,y1:y2,x1:x2],dtype=float) #nb data creation
            self.narrowband = np.sum(cutout,axis=0) #I sum over wavelengths (axis 0) so I am left with a matrix of 2D (spatial)
            #to create a nb we need to do a plain unweighted summation, which is the default approach:
            #nb = np.sum( data[l0,j0,i0:l1,j1,i1], axis=0 )
            #where l0, l1 are the limits in wavelength direction and j0, j1, i0, i1 specifies 
            #the window in spatial direction. Recall that in the numpy convention, the first index 
            #is always included, the second index is not (so if l1 is the last spectral layer to be 
            #included, you have to extract up to l1+1).
    
        else:
            cutout = np.asarray(cube_data[lam1_pix:lam2_pix],dtype=float)
            self.narrowband = np.sum(cutout,axis=0)
     
        if np.size(self.narrowband)!=0: #if the narrowband array is not empty          
        
            # update fits header to show correct positions
        
            self.narrowband_header = copy.deepcopy(cube_header)
        
            if NBwindow_size != 'full':
                crpix1_new = int(window_size_pix_half)
                crpix2_new = int(window_size_pix_half)
                self.narrowband_header['CRPIX1'] = crpix1_new #ref px point along x
                self.narrowband_header['CRPIX2'] = crpix2_new #ref px point along y
            crval1_new = float(ra_middle)
            crval2_new = float(dec_middle)
        
            self.narrowband_header['CRVAL1'] = crval1_new #first ra value at ref px point along x
            self.narrowband_header['CRVAL2'] = crval2_new #first dec value at ref px point along y
        
            # remove third dimension
            self.narrowband_header['NAXIS'] = 2 #number of dimensions
            self.narrowband_header['WCSAXES'] = 2 
            del self.narrowband_header['CRVAL3']
            del self.narrowband_header['NAXIS3']
            del self.narrowband_header['CTYPE3']
            del self.narrowband_header['CUNIT3']
            del self.narrowband_header['CRPIX3']        
            del self.narrowband_header['CD1_3']
            del self.narrowband_header['CD2_3']
            del self.narrowband_header['CD3_1']
            del self.narrowband_header['CD3_2']
            del self.narrowband_header['CD3_3'] 
     
            #save the NB image in a fits file
            #hdu = fits.PrimaryHDU(self.narrowband,header=self.narrowband_header)
            #hdul = fits.HDUList([hdu])
            #hdul.writeto(dict_cats['MUSEnarrowband'],overwrite=True)  #if I want to save a fits file with the NB data, header
            print('NB created from',MUSE_DATA_CUBE)
            #hdul.close()
        
            '''
            #plot the NB to make sure we are selecting a galaxy
            w_cutout = WCS(self.narrowband_header)
            #DELETE UP TO SHOW ONCE I GET THE NB TO PROPERLY SHOW IN GUI
            fig, ax = plt.subplots(figsize=(9., 7.))
            ax = plt.subplot(projection=w_cutout)  #shows axes in h, min, s
            # get useful ranges for the brightness
            here_not_0 = (self.narrowband!=0)*(self.narrowband>-99999)
            perc = 99.5
            dat_perc = np.percentile(self.narrowband[here_not_0],perc)
            med_hst = np.median(self.narrowband[here_not_0])
            vmin,vmax = med_hst-dat_perc,med_hst+dat_perc
            im = ax.imshow(self.narrowband, origin='lower',vmin=vmin,vmax=vmax)
            x,y= w_cutout.wcs_world2pix(ra,dec,0)    
            ax.scatter(x,y, s=300, edgecolor='k', facecolor='none', linewidths=2) #black circle around x,y (galaxy position)
            ax.scatter(x,y, s=300, color='k', marker='+', linewidths=3)
            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')
            plt.show()   
            '''
            #The NB extraction from cube is finished
            
            self.muse_data, self.header_muse = self.narrowband, self.narrowband_header 
            
            # find out whether cdfs or cosmos from ID       
            area = (self.id_here)[0]  
            if area == '1':  
                self.area = 'cdfs'
            elif area == '2':
                self.area == 'cosmos'
            field = (self.id_here)[1:3] # find out the field number from the ID
            if field != self.field:
                self.field = field
                self.get_HST_images() #function defined above that passes the HST images
                
            self.UV_ID = '0' #0 is just to initialize it, after I click in each object I get the ID written in the little box
           
            #lets make a circle with an inner cross around the muse position of the gal
            self.MUSE_x,self.MUSE_y = radec_to_pix(self.MUSE_RA,self.MUSE_Dec,self.hst_header,origin=0) #hst_header or the one from muse nb???
        
        else: #if the narrowband is empty
            for i, item in enumerate(MUSE_DATA_CUBE):
                if item=='-':
                    print('Warning! The galaxy you are trying to access is at RA=',ra,'and DEC=',dec,'Is it outside the field',MUSE_DATA_CUBE[i+1:i+8],'? Then, write the correct name for MUSE_DATA_CUBE and HST images.')
                    self.close()
        
    # function to know which galaxies can be counterparts       
    # make possible counterparts clickable in the main hst big image and surround them in white by default until they are clicked
    def get_pos_UV(self):
                
        self.UV_RAs = info_UV['RA']
        self.UV_Decs = info_UV['DEC']
        self.UV_IDs = info_UV['ID']
        self.UV_x,self.UV_y = radec_to_pix(self.UV_RAs,self.UV_Decs,self.hst_header,origin=0)
        #I read the z from a different catalog than from where I take RA, Dec, ID
        #UV2_RAs = info_UV2['RA'] uncomment when we have the catalog
        #UV2_Decs = info_UV2['DEC'] uncomment when we have the catalog
        #self.UV2_IDs = info_UV2['ID'] uncomment when we have the catalog
        #self.UV_confidences = zconfidence #uncomment when we have the catalog
        #self.UV_zs = z #uncomment when we have the catalog
        self.UV_confidences = info_UV['RA'] #delete when we have the catalog!!!!!!!!!!!!!!!!!!
        self.UV_zs = info_UV['RA'] #delete when we have the catalog!!!!!!!!!!!!!!!!!!
        
        # find objects close to objects in other catalogue
        margin = 200 # how far away to look for counterparts in pixels to be displayed on the image 
        dist = np.sqrt((self.x_hst-self.UV_x)**2+(self.y_hst-self.UV_y)**2)
        
        close_bool = dist<margin        
        close_here_ID = self.UV_IDs[close_bool]
        close_here_x = self.UV_x[close_bool]
        close_here_y = self.UV_y[close_bool]

        # brush is transparent, only pen is shown as thick line
        self.s1 = pg.ScatterPlotItem(close_here_x,close_here_y,size=32,pxMode=False,pen=pg.mkPen('w',width=3),
                                     brush=pg.mkBrush(255,255, 255, 0),name=close_here_ID) #pen is the thick line and brush the filling (empty)
        
        self.s2 = pg.ScatterPlotItem([float(self.MUSE_x)],[float(self.MUSE_y)],symbol='o',pxMode=False,pen=pg.mkPen('k',width=2),size=15,brush=pg.mkBrush(255,255, 255, 0)) 
        self.s3 = pg.ScatterPlotItem([float(self.MUSE_x)],[float(self.MUSE_y)],symbol='+',pxMode=False,size=10) 
        
        self.view_hst_main.addItem(self.s1)
        self.view_hst_main.addItem(self.s2)
        self.view_hst_main.addItem(self.s3)
        # Make white points clickable
        self.s1.sigClicked.connect(self.clicked_counterpart)
       
        
    #get the information from the selected objects in the get_pos_UV function 
    #this info will be then displayed with the fill_window function   
    #make counterparts red once clicked
    def clicked_counterpart(self,plot,points):
        xy = points[0].pos()
        x,y = xy[0],xy[1]
        # retrieve the ID of the clicked object
        here_ID_UV = int(self.UV_IDs[np.argmin(abs(self.UV_x-x))])
        self.UV_ID = str(here_ID_UV)
        self.le_UV_ID.setText(self.UV_ID) #I write the ID to the little box with a 0 as a default
        
        #not displayed in gui but maybe useful to save in table?
        #here_ID_UV2 = int(self.UV2_IDs[np.argmin(abs(self.UV2_x-x))]) #uncomment when we have the catalog
        #self.UV2_ID = str(here_ID_UV2)#uncomment when we have the catalog
        
        here_DEC_UV = float("%.9f" % self.UV_Decs[np.argmin(abs(self.UV_x-x))])
        self.UV_DEC = here_DEC_UV #this is the DEC from the photometry
        delta = abs(float(self.MUSE_Dec)-float(self.UV_DEC))
        delta_DEC = str(delta)
        self.delta_DEC = '%.2E' % Decimal(delta_DEC)
        self.le_UV_DEC.setText(self.delta_DEC)
        
        here_RA_UV = float("%.9f" % self.UV_RAs[np.argmin(abs(self.UV_x-x))])
        self.UV_RA = here_RA_UV #this is the RA from the photometry
        delta = abs((float(self.MUSE_RA)-float(self.UV_RA))*math.cos(self.MUSE_Dec*math.pi/180))
        delta_RA = str(delta)
        self.delta_RA = '%.2E' % Decimal(delta_RA)
        self.le_UV_RA.setText(self.delta_RA)
        
        here_z_UV = float("%.3f" % self.UV_zs[np.argmin(abs(self.UV_x-x))])#delete when we have the catalog!!!!!!!!!!!!!!!!!!
        self.UV_z = str(here_z_UV)#delete when we have the catalog!!!!!!!!!!!!!!!!!!
        here_confidence_UV = int(self.UV_confidences[np.argmin(abs(self.UV_x-x))])#delete when we have the catalog!!!!!!!!!!!!!!!!!!
        self.UV_confidence = str(here_confidence_UV)#delete when we have the catalog!!!!!!!!!!!!!!!!!!
        #if isinstance(self.UV_zs[np.argmin(abs(self.UV2_x-x))], str) ==True:#uncomment when we have the catalog
        #   self.UV_z='-' #float(0) 
        #   self.le_UV_z.setText(self.UV_z)
        #   self.UV_confidence.setText(self.UV_z)
        #else:
        #   here_z_UV = float("%.3f" % self.UV_zs[np.argmin(abs(self.UV2_x-x))]) 
        #   self.UV_z = str(here_z_UV) 
        #   self.le_UV_z.setText(self.UV_z) 
           #here_confidence_UV = float("%.3f" % self.UV_zs[np.argmin(abs(self.UV2_x-x))])#uncomment when we have the catalog
           #self.UV_z = str(here_confidence_UV)#uncomment when we have the catalog
           #self.le_UV_confidence.setText(self.UV_confidence)#uncomment when we have the catalog
           
                
        # change colour of circle when clicked from white to red
        try:
            # change colour of previously clicked circle back
            self.last_clicked.setPen('w',width=3)
        except:
            pass
        points[0].setPen('r',width=3)
        
        self.last_clicked = points[0]




    #after we click in the button next or previous we call this function to save the clicked counterparts 
    #for the specific MUSE object as well as their info    
    def show_next(self):     
   
        # move to next or previous object
        where = self.sender().text()  # label of the button that was clicked
                
        #if we have clicked on any button of the GUI and on next only, we save data to output file
        #This saves the gal in random order           
        if where == 'next' and (int(self.le_UV_ID.text()) !=0 or type(self.RA_WithoutCounterpart)!=str or self.le_comment.text() != '' or self.le_confidence !='0' or self.le_nomatch != '-' or self.le_nocatalog != '-' or self.le_more != '-') :
            # write counterpart info in table            
            if self.le_UV_ID.text() != '0': #there is a counterpart
                separation = np.sqrt((float(self.MUSE_RA)-float(self.UV_RA))**2+(float(self.MUSE_Dec)-float(self.UV_DEC))**2) #separation between MUSE and the photometry object only if we have chosen a counterpart
                #muse id where we is given by int(self.id_here)
                counterparts[str(self.id_here)] = {'MUSE_RA':float(self.MUSE_RA),'MUSE_DEC':float(self.MUSE_Dec),'MUSE_z':float(self.MUSE_Z),
                                                       'UV_ID':str(self.le_UV_ID.text()),'UV_RA':float(self.UV_RA),'UV_DEC':float(self.UV_DEC),
                                                       'delta_RA':float(self.le_UV_RA.text()),'delta_DEC':float(self.le_UV_DEC.text()), 
                                                       'UV_z':float(self.le_UV_z.text()), 'Separation':separation, 
                                                       'Confidence':self.le_confidence,'No match':self.le_nomatch,'ra_noMatch': self.RA_WithoutCounterpart,
                                                       'dec_noMatch': self.Dec_WithoutCounterpart,'Photometry':self.forced_photometry,'Not in catalog':self.le_nocatalog, 
                                                       'More than one counterpart':self.le_more,                                           
                                                       'Comment':self.le_comment.text()}
            else: #if there is no counterpart
                if self.RA_WithoutCounterpart != '-': #if we have clicked somewhere on the HST image to force the photometry   
                    aperture=CircularAperture([(self.x_WithoutCounterpart, self.y_WithoutCounterpart)], r=20) #we perform the photometry in an aperture of radius of 20px=0.6arcsec
                    forced_photometry_table_instrumental=aperture_photometry(np.rot90(self.hst_data)[::-1], aperture)
                    forced_photometry=np.array(forced_photometry_table_instrumental['aperture_sum'])[0] #the phtometry is given in instrumental flux units (electrons/s)
            
                    #flux conversion to physical units
                    if self.hst_data.all() == self.hst_data_436.all(): #if the main HST is the band 435W
                        #to convert it to flux physical units: see https://stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
                        forced_photometry_inf=forced_photometry/0.941 #is the correction from .6 arcsec to inf of the filters (see first graph in https://stsci.edu/hst/instrumentation/acs/data-analysis/aperture-corrections)
                        q=acszpt.Query(date='2016-02-22',detector='WFC',filt='F435W') #date I take it from the HST headers
                        filter_zpt=q.fetch()
                        F_forced_photometry=forced_photometry_inf*filter_zpt['PHOTFLAM'].value[0] #flux density in physical units, erg/s/cm^2/Angstrom 
                    elif self.hst_data.all() == self.hst_data_606.all():
                        forced_photometry_inf=forced_photometry/0.947
                        q=acszpt.Query(date='2016-02-22',detector='WFC',filt='F606W') 
                        filter_zpt=q.fetch()
                        F_forced_photometry=forced_photometry_inf*filter_zpt['PHOTFLAM'].value[0]
                    elif self.hst_data.all() == self.hst_data_775.all(): 
                        forced_photometry_inf=forced_photometry/0.949
                        q=acszpt.Query(date='2016-02-22',detector='WFC',filt='F775W') 
                        filter_zpt=q.fetch()
                        F_forced_photometry=forced_photometry_inf*filter_zpt['PHOTFLAM'].value[0]
                    else: 
                        forced_photometry_inf=forced_photometry/0.940
                        q=acszpt.Query(date='2016-02-22',detector='WFC',filt='F435W') #there is no 160 band here
                        filter_zpt=q.fetch()
                        F_forced_photometry=forced_photometry_inf*filter_zpt['PHOTFLAM'].value[0]
                    counterparts[str(self.id_here)] = {'MUSE_RA':float(self.MUSE_RA),'MUSE_DEC':float(self.MUSE_Dec),'MUSE_z':float(self.MUSE_Z),
                                                       'UV_ID':int(0),'UV_RA':int(0),'UV_DEC':int(0),'delta_RA':int(0),'delta_DEC':int(0), 'UV_z':int(0), 'Separation':'-', 
                                                       'Confidence':self.le_confidence,'No match':self.le_nomatch,'ra_noMatch': float(self.RA_WithoutCounterpart),
                                                       'dec_noMatch': float(self.Dec_WithoutCounterpart),'Photometry':F_forced_photometry,'Not in catalog':self.le_nocatalog, 
                                                       'More than one counterpart':self.le_more,                                           
                                                       'Comment':self.le_comment.text()}
                else:
                    counterparts[str(self.id_here)] = {'MUSE_RA':float(self.MUSE_RA),'MUSE_DEC':float(self.MUSE_Dec),'MUSE_z':float(self.MUSE_Z),
                                                       'UV_ID':int(0),'UV_RA':int(0),'UV_DEC':int(0),'delta_RA':int(0),'delta_DEC':int(0), 'UV_z':int(0), 'Separation':'-', 
                                                       'Confidence':self.le_confidence,'No match':self.le_nomatch,'ra_noMatch': self.RA_WithoutCounterpart,
                                                       'dec_noMatch': self.Dec_WithoutCounterpart,'Photometry':self.forced_photometry,'Not in catalog':self.le_nocatalog, 
                                                       'More than one counterpart':self.le_more,                                           
                                                       'Comment':self.le_comment.text()}
                    
                    
            for i, item in enumerate(counterparts):
                position=i #number of objects we have in counterparts
            
            
            while True: #we need while to use break
                try: #if I cannot open the output file (it doesnt exist), then I create the table
                    hdu = fits.open(dict_cats['output'])
                    hdu.close() 
                    rwf.modify_output_table(counterparts,dict_cats['output'],1,position)
                    break
                except FileNotFoundError: #if trying to open the file we get this error (bc it doesnt exist)
                    #we create a new table from scratch
                    rwf.write_table_from_dict_in_dict(counterparts,dict_cats['output'])
                    print('Created output file',dict_cats['output'])
                    break     
        else: #if we only click next or previous, we save nothing        
            print('Nothing has been saved')
        
        # close if all objects were looked at
        max_len = len(ids_all)  
       
        if self.j >= max_len:
            print('Done! Good job ;)')
            self.close()
        else:
            # check if going to next object or previous
            if where == 'previous':
                direction = -1
            else: 
                direction = 1
                
            # clear all Widgets for next object    
            for i in reversed(range(self.grid.count())): 
                self.grid.itemAt(i).widget().setParent(None)            
            
            self.j+=direction # move forward. I change the counter to where I am
            print('position2',self.j)
            if self.j<0:
                print('There is no previous object! You tried to reach position',self.j)
                self.close()
            else:
                self.get_infos() # get info on new ID
                #self.get_HST_images()
                self.fill_window() # show main window again
                
        

def main():
    
    app = QtGui.QApplication(sys.argv)
    main = main_GUI()
        
    main.show()    
    app.exec_() # start the app

if __name__ == '__main__':

    main()


