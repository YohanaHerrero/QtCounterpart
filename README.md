# QtCounterpart

This README documents the necessary steps to get "QtCounterpart" running.

## **What is this repository for?**

QtCounterpart is a GUI that helps find HST counterparts for sources detected in integral field data cubes, i.e., three dimensional data with two spatial and one wavelength dimensions. QtCounterpart was build in Python for the HST counterpart assignation of detected objects in MUSE data cubes. The GUI takes each detected MUSE source and guesses which is the corresponding counterpart in HST data, if any. QtCounterpart can display a narrowband from the data cube and four HST bands for each spectroscopically detected object. In addition, the spectroscopic and photometric information of the object can also be displayed. As input you need a MUSE datacube, its corresponding spectroscopic catalogue of sources (at least with spatial and wavelenght/redshift data), HST band images (maximum four), and a photometric catalogue (at least with spatial information). 

## **Installing QtCounterpart**

QtCounterpart does not require any installation. Simply cloning the QtCounterpart GitHub repository and importing the scripts should be enough. Hence, QtCounterpart is "installed" by doing:

```
    git clone https://github.com/YohanaHerrero/QtCounterpart.git
```

QtCounterpart is written in Python and uses a range of default packages included in standard installations of Python. In order to run QtCounterpart you will need python version 3.7.x and the following packages:

### **Standard packages:**

- numpy  
- matplotlib
- math
- scipy 
- pylab

### **Special packages and input files:**

- acstools
- photutils
- astropy 
- pyqtgraph 
- PyQt5

Besides the python packages, you will also need an integral field data cube, HST bands of the same field, a catalogue of detected MUSE sources and a photometric catalog. Note that only the catalogs and one HST image are *must* arguments to the GUI, the rest is optional.
The MUSE and photometric catalogue must contain the following columns (they can also contain additional info such as ID, confidence, lead line, SN, flux/magnitude, photo_z):

- RA in degrees
- DEC in degrees
- Redshift or wavelength (if this info is no provided, no narrowband image from the data cube is displayed but the GUI runs)

ident_counterparts.py is the main programme for QtCounterpart. By doing ident_counterparts.py --help, the input files and parameters you should pass to the GUI are shown. It also shows the default names of the files (and parameters), their location, and what they exactly are. E.g:
```
usage: ident_counterparts.py [-h] [-cube INPUTDATA] [-c CATALOG] [-p PCATALOG]
                             [-hstMAIN HSTIMAGEMAIN] [-hst1 HSTIMAGE1]
                             [-hst2 HSTIMAGE2] [-hst4 HSTIMAGE4] [-o OUTPUT]
                             [-c_RA COLUMN_RA] [-c_DEC COLUMN_DEC]
                             [-c_lambda COLUMN_LAM] [-c_z COLUMN_Z]
                             [-c_ID COLUMN_ID] [-c_sn COLUMN_SN]
                             [-c_line COLUMN_LINE]
                             [-c_confidence COLUMN_CONFIDENCE]
                             [-c_pRA COLUMN_PRA] [-c_pDEC COLUMN_PDEC]
                             [-c_pID COLUMN_PID] [-c_photoz COLUMN_PHOTOZ]
                             [-c_pmag_flux COLUMN_PMAG_FLUX]

ident_counterparts.py - identify HST counterparts with a GUI!

optional arguments:
  -h, --help            show this help message and exit
  -cube INPUTDATA, --inputdata INPUTDATA
                        Input flux cube FITS file (e.g. DC subtracted
                        effnoised cube). (default: None)
  -c CATALOG, --catalog CATALOG
                        Input spectroscopic catalog FITS file (e.g. from
                        LSDCat, consolidated MUSE catalog, etc). (default: /Us
                        ers/yohana/Desktop/PhD_Yohana/Galaxy_counterpart_class
                        ification/Counterpart_GUI/Fields/lsdcat_candels-
                        cdfs-47.fits)
  -p PCATALOG, --pcatalog PCATALOG
                        Input photometric catalog FITS file. (default: /Users/
                        yohana/Desktop/PhD_Yohana/Galaxy_counterpart_classific
                        ation/Counterpart_GUI/catalog_photometry_candels-
                        cdfs-47.fits)
  -hstMAIN HSTIMAGEMAIN, --hstimageMAIN HSTIMAGEMAIN
                        HST main image for comparison. (default:
                        HST_images/acs_814w_candels-cdfs-47_cut_v1.0.fits)
  -hst1 HSTIMAGE1, --hstimage1 HSTIMAGE1
                        HST1 image for comparison. (default:
                        HST_images/acs_435w_candels-cdfs-47_cut_v1.0.fits)
  -hst2 HSTIMAGE2, --hstimage2 HSTIMAGE2
                        HST2 image for comparison. (default:
                        HST_images/acs_606w_candels-cdfs-47_cut_v1.0.fits)
  -hst4 HSTIMAGE4, --hstimage4 HSTIMAGE4
                        HST4 image for comparison. (default:
                        HST_images/wfc3_160w_candels-cdfs-47_cut_v1.0.fits)
  -o OUTPUT, --output OUTPUT
                        Output catalog FITS table. (default: /Users/yohana/Des
                        ktop/PhD_Yohana/Galaxy_counterpart_classification/Coun
                        terpart_GUI/output_test.fits)
  -c_RA COLUMN_RA, --column_RA COLUMN_RA
                        Column name in input FITS spectroscopic catalog for RA
                        (in degrees), default=RA. (default: RA)
  -c_DEC COLUMN_DEC, --column_DEC COLUMN_DEC
                        Column name in input FITS spectroscopic catalog for
                        DEC (in degrees), default=DEC. (default: DEC)
  -c_lambda COLUMN_LAM, --column_LAM COLUMN_LAM
                        Column name in input FITS spectroscopic catalog for
                        wavelength (in Angstrom), default=LAMBDA. (default:
                        LAMBDA)
  -c_z COLUMN_Z, --column_z COLUMN_Z
                        Column name in input FITS spectroscopic catalog for
                        redshift, default=Z. (default: Z)
  -c_ID COLUMN_ID, --column_ID COLUMN_ID
                        Column name in input FITS spectroscopic catalog for
                        ID, default=ID. (default: ID)
  -c_sn COLUMN_SN, --column_sn COLUMN_SN
                        Column name in input FITS spectroscopic catalog for
                        SN, default=S2N. (default: S2N)
  -c_line COLUMN_LINE, --column_LINE COLUMN_LINE
                        Column name in input FITS spectroscopic catalog for
                        lead line, default=LINE_ID. (default: LINE_ID)
  -c_confidence COLUMN_CONFIDENCE, --column_confidence COLUMN_CONFIDENCE
                        Column name in input FITS spectroscopic catalog for
                        confidence, default=CONFIDENCE. (default: CONFIDENCE)
  -c_pRA COLUMN_PRA, --column_pRA COLUMN_PRA
                        Column name in input FITS photometric catalog for RA
                        (in degrees), default=ra. (default: ra)
  -c_pDEC COLUMN_PDEC, --column_pDEC COLUMN_PDEC
                        Column name in input FITS photometric catalog for DEC
                        (in degrees), default=dec. (default: dec)
  -c_pID COLUMN_PID, --column_pID COLUMN_PID
                        Column name in input FITS photometric catalog for ID,
                        default=id. (default: id)
  -c_photoz COLUMN_PHOTOZ, --column_photoz COLUMN_PHOTOZ
                        Column name in input FITS photometric catalog for
                        photo_z, default=z_p. (default: z_p)
  -c_pmag_flux COLUMN_PMAG_FLUX, --column_pmag_flux COLUMN_PMAG_FLUX
                        Column name in input FITS photometric catalog for
                        pflux or magnitude (will be saved in output file),
                        default=F160W_LIMITING_MAGNITUDE. (default:
                        F160W_LIMITING_MAGNITUDE)
```

## **Running QtCounterpart**

Once you have gone through ident_counterparts.py --help and you know which files and parameters you are going to pass to the GUI, you can specify default directories and names of the necessary files in the first lines of code of ident_counterparts.py. In this way, you do not need to pass each single argument every time you run the GUI and the default names are used instead. E.g.:
```
###########################################################################################
############################## Change default file names here #############################

directory = '/Users/yohana/Desktop/PhD_Yohana/Galaxy_counterpart_classification/Counterpart_GUI/' #main directory (where scripts are saved)
field_dir = 'Fields/' #directory of MUSE cube and spectrocopic catalog
hst_dir = 'HST_images/' #directory of HST images
def_catalog = (directory+field_dir+'lsdcat_candels-cdfs-47.fits') #default spectroscopic catalog
def_pcatalog = (directory+'catalog_photometry_candels-cdfs-47.fits') #default photometric catalog
def_output = (directory+'output_test.fits') #default output table
#default HST images
def_hstMAIN = (hst_dir+'acs_814w_candels-cdfs-47_cut_v1.0.fits') 
def_hst1 = (hst_dir+'acs_435w_candels-cdfs-47_cut_v1.0.fits') 
def_hst2 = (hst_dir+'acs_606w_candels-cdfs-47_cut_v1.0.fits') 
def_hst4 = (hst_dir+'wfc3_160w_candels-cdfs-47_cut_v1.0.fits') 
```

After inserting your files, QtCounterpart is ready to run. An example of running the GUI with all possible files would be:
```
ident_counterparts() -cube /Users/yohana/Desktop/PhD_Yohana/Galaxy_counterpart_classification/Counterpart_GUI/Fields/DATACUBE_candels-cdfs-47_v2.0_dcsub_effnoised.old.fits -c_pRA RA -c_pDEC DEC -c_pID ID
```
In this case, only the cube is passed as a file argument, since I use the default files for the photometric/spectrocopic catalog, output table and the HST images (i.e. def_catalog, def_pcatalog, def_output, def_hstMAIN, def_hst1, def_hst2, def_hst4). The passed arguments for the photometric catalog -c_pRA RA -c_pDEC DEC -c_pID ID mean that the parameters in the photometric catalog are different from the default values. If they were the same, I would not have to pass them as arguments. Thus, if an argument is not passed is because you want to use the default value. If you do not want to use the file/parameter at all, you need to delete the part "default = X" of the lines below, where X if the default name of the files i.e. def_catalog, def_pcatalog, etc.
```
#Input cube, spectroscopy and photometric catalog, HST images, and output file
parser.add_argument("-cube","--inputdata", required=False, type=str,
                    help="Input flux cube FITS file (e.g. DC subtracted effnoised cube).")

parser.add_argument("-c","--catalog", required=False, type=str,
                    *default=def_catalog*, help="Input spectroscopic catalog FITS file (e.g. from LSDCat, consolidated MUSE catalog, etc).")

parser.add_argument("-p","--pcatalog", required=False, type=str,
                    *default=def_pcatalog*, help="Input photometric catalog FITS file.")
etc...
```

Once, the script is run, the GUI is displayed.

### **What do I see?**
The bottom left big panel shows the narrowband of the objects detected by MUSE ("MUSE narroband"), which are saved in the MUSE catalogue. The top left big panel shows a HST band cutout ("Main HST image") corresponding to the same region where the detected MUSE object is located. Similarly, the four top small panels (HST band 435, 606, 775/814, and 160) show different HST bands of the same main HST image cutout. Inside the main HST image, the black circle with the inner cross represents the MUSE position of the object, while the empty white circles show the position of the possible HST counterparts of the MUSE source. The numerical boxes below the small HST panels show the current resolution cut of the HST bands, representing the scaling factor used for the intensity of the HST images. The "Reset images" button offers the possibility to set all GUI images back to the original display. Spectroscopic data from MUSE is shown in the left half of the GUI. There, information such as MUSE ID, right ascension (RA), declination (Dec), wavelength (lambda), redshift (z), confidence of the source detection, the line that lead to the detection or the signal to noise (S/N) of the detection is found. Right above, photometric data from the phototometric catalogue is shown. The photometric UV ID, the difference in RA between the spectroscopic and the photometric values, same for Dec, and the photometric redshift are displayed. 

![GUI_screenshot](https://user-images.githubusercontent.com/81447012/207744202-4691e966-37ca-41a2-8af2-4f20ef12b32a.PNG)


### **What can I do?**

#### **General overview** 
In the images you can compare the source detections in MUSE with the corresponding ones observed in different HST bands. Besides, the possible counterparts to the MUSE object are marked with a white circle in the main HST image, based on spatial assumptions (max 6" apart from the MUSE source). You can select the most likely counterpart by clicking on those circles and click on the confidence buttons to record how sure you were about the decision. Then, the information of that counterpart will be recorded in the photometric information section. By clicking on the "next" button you will save all the displayed information in a fits output Table. 

#### **Script overview** 

There are several scripts that the main GUI script, ident_counterpart.py, uses:
- ```colormaps.py``` : matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt, and (in the case of viridis) Eric Firing
- ```colors.py```: get colors for the GUI images from colormaps.py
- ```read_write_files.py```: read, save and modify fits tables by creating dictionaries with information using the first column as identifiers. 

As already mentioned, the main script is ident_counterpart.py. Below, I explain the purpose of each of button and the general functionalities of the GUI. Also, by leaving the mouse on top of each button/box text a tip box is displayed to remind the user of the functionality of the button:
The contrast of the small HST images can be modified by writing a different cut number on the cut boxes and clicking afterwards on the "Cuts" button. The smaller that number is, the lower the resolution of the images. By using the "Reset images" button, the resolution of all HST images is set back to default (0.99).
Now, if one is more interested in one spectific HST band than in the rest of them, by double clicking on the wanted small HST image, you will convert it into the main HST image of the GUI. You can always zoom in or out all images in the GUI by using the mouse.

The main HST image already hints at possible counterparts (only when existing in the photometric catalogue) for the currently displayed MUSE source. You can select any of the white circles to choose the most likely counterpart. This circle will then turn into red and the photometric information of the selected object will be displayed in the photometric information section. In order to guide the eye during the counterpart selection, the main HST image also contains the position of the MUSE object according to the spectroscopic data (shown in the MUSE information section). Once you have selected the counterpart, you should tell the GUI how sure you are about your selection. For that, you can click on the "Confidence" bullet points, indicating "0" for the most unsure case and "3" for the most sure scenario. You have then found an HST counterpart for the MUSE source so you can click on "next" to save all information and find the counterpart for the next object in the MUSE catalogue. 

There might be of course cases for which the MUSE object has no counterpart in the HST data, the object is not registered in the photometric catalogue or there seems to be more than one counterpart. These cases can be recorded by clicking in the empty boxes below the "Confidence" region. You can also write a comment so it is easier to recall the case in future inspections. If you believe you made a mistake selecting the counterpart or choosing any of these possible scenarios, you only need to click on the "Clear counterparts" button to start classifying that MUSE source from scratch.

In case that the counterpart is not in the photometric catalogue, there will be no white circle around that galaxy in the main HST image. However, if you believe you see a counterpart in the main HST band, it is possible to perform a forced photemetry by cliking (in the main HST image) on the center of the galaxy that you believe is the counterpart. Once you have clicked on the area of the main HST image where that photometrically undetected counterpart is, the GUI measures the flux and coordinates of the object. Specify the case that you encountered (e.g. "Counterpart not in catalog"), leave a comment if you wish and click "next" to save the information.

Anytime it is possible to go back to the already classified objects by clicking on "previous" until you find the desired object. Remember that if you want to save or rewrite the classification of an object, you must click on "next" after you have clicked on any button/counterpart of the GUI. A more practical way of doing this is to write the MUSE ID of the desired object in the empty box right the "Jump to MUSE ID" button, then click on the button and your desired object will be shown in the interface.

Let's now say you do not want to classify more objects. Then, press the red cross on the top right corner of the display to exit the GUI. Once you want to classify more objects again, you only need to run the GUI and click on the "Jump to next unclassified ID" button to continue with the object where you left it.

The most important functions in the script are described below, but for further details, refer to the header of each routine.

- ```get_HST_images()```: extract HST and MUSE data and headers from the input fits files
- ```get_MUSE_info()```: create MUSE NB from the MUSE data cube and extract spectroscopic information to be displayed under "MUSE information"
- ```get_UV_info()```: extract photometric information to be displayed under "Photometric information" and locations of possible HST counterparts
- ```image_HST*i*()```, ```main_image_HST()```, ```image_narrowband_muse()```: cut the relevant area of the small HST band *i* or of the datacube and display it. Same for the main HST image
- ```clicked_counterpart()```: display information of the clicked counterpart
- ```clickToObtainPositionOnImage()```: allow to click anywhere on the main HST image to extract the position of clicked point (used in the forced photometry when the MUSE object is not recorded in the photometric catalogue)
- ```fill_window()```: display buttons and images on the screen
- ```show_next()```: only called when clicking on "next". Save displayed information and perform forced photometry (if needed)

## **QtCounterpart output**

The output file generated by QtCounterpart is a fits Table containing all the counterpart identification results. The table contains the following columns:

- ```'ID'```: ID of the MUSE object
- ```'UV_ID'```: photometric ID of the counterpart
- ```'MUSE_RA'```: spectroscopic RA of the MUSE object (in deg)
- ```'UV_RA'```: photometric RA of the counterpart (in deg)
- ```'delta_RA'```: difference between MUSE_RA and UV_RA
- ```'MUSE_DEC'```: spectroscopic DEC of the MUSE object (in deg)
- ```'UV_DEC'```: photometric DEC of the counterpart (in deg)
- ```'delta_DEC'```: difference between MUSE_DEC and UV_DEC
- ```'Separation'```: angular separation in rad between spectroscopic object and photometric location of the counterpart
- ```'MUSE_lambda'```: spectroscopic wavelength of the MUSE object
- ```'MUSE_z'```: spectroscopic redshift of the MUSE object
- ```'UV_z'```: photometric redshift of the counterpart
- ```'Photometry'```: flux or magnitude from photometric catalog
- ```'Confidence'```: confidence of your choice e.g. how sure are you of your counterpart (non-)selection?
- ```'Comment'```: comment written during the selection, if any
- ```'ra_noMatch'```: RA of the undetected counterpart in photometric catalogue, in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image
- ```'dec_noMatch'```: DEC of the undetected counterpart in photometric catalogue, in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image
- ```'Forced photometry'```: Forced photometry in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image (units: erg/s/cm^2/Angstrom)
- ```'More than one counterpart'```: 'Yes' if there is more than one counterpart and '-' if only one or none
- ```'No match'```: 'Yes' if there is no HST detection that matches the MUSE object and '-' if there is one/no counterpart 
- ```'Not in catalog'```: 'Yes' if the counterpart is not in the photometric catalogue but you see a likely counterpart without a white circle around, and '-' for other cases e.g. there is one/no counterpart


## **Tips for users**

When opening QtCounterpart, patience is needed. The datacubes are huge and it takes some time to load and extract the narrowband out of it. You need to make sure though that you used the right input and output catalogue file names, with the required column names and information, otherwise the GUI will crash.

Comments, issues and/or suggestions for improvement, new functions etc. are more than welcome. Please send these to my working email (yherreroalonso at aip.de) or add an 'issue' via GitHub.
