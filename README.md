# QtCounterpart

This README documents the necessary steps to get "QtCounterpart" running.

**What is this repository for?**

QtCounterpart is a GUI that helps find HST counterparts for sources detected in integral field data cubes, i.e., three dimensional data with spatial and wavelength dimensions. QtCounterpart was build in Python for the HST counterpart assignation of detected objects in MUSE data cubes. The GUI takes each detected MUSE source and guess which is the corresponding counterpart in HST data, if any. QtCounterpart displays a MUSE narrowband and four HST bands of the object. In addition, its spectroscopic and photometric information is also shown. As input you need a MUSE datacube, its corresponding catalogue of sources (with spatial and redshift data), four HST band images, and a photometric catalogue. 

**Installing QtCounterpart**

QtCounterpart does not require any installation. Simply cloning the QtCounterpart GitHub repository and importing the scripts should be enough. Hence, QtCounterpart is "installed" by doing:

    git clone https://github.com/YohanaHerrero/QtCounterpart.git

QtCounterpart is written in Python and uses a range of default packages included in standard installations of Python. In order to run QtCounterpart you will need python version 3.7.x and the following packages:

**Standard packages:**

numpy  
matplotlib
math
scipy 
pylab

**Special packages:**

acstools
photutils
astropy 
pyqtgraph 
PyQt5

After adding the QtCounterpart directory to the PYTHONPATH or changing location to the QtCounterpart directory, QtCounterpart can be imported in python with:

    import QtCounterpart

You will also need an integral field data cube, four HST bands of the same field, a catalogue of detected MUSE sources and a photometric catalog.
The MUSE and photometric catalogues will have to contain the following columns:

ID
RA in degrees
DEC in degrees
Z 

The MUSE catalogue should further include the columns:
Confidence
Lead_line
S/N

Now, in the first lines of code of ident_counterparts.py, which is the main programme for QtCounterpart, you can specify the directories and names of the necessary files. E.g.

MUSE_DATA_CUBE = 'DATACUBE_candels-cdfs-01_v2.0_dcsub_effnoised.fits' #MUSE cube

dict_cats = {
             'MUSE':'MW_60fields_main_table_v1.1.fits', # MUSE catalogue by default 
             'UV':'CANDELS.GOODSS.F160W.v1_1.photom.cat', # photometric catalogue 
             'output':'test_output.fits', # output table where counterparts info is saved
             'HST1':'HST_images/acs_435w_candels-cdfs-01_cut_v2.0.fits', # HST images inside HST directory
             'HST2':'HST_images/acs_606w_candels-cdfs-01_cut_v2.0.fits',  
             'HST3':'HST_images/acs_814w_candels-cdfs-01_cut.fits',  
             'HST4':'HST_images/wfc3_160w_candels-cdfs-01_cut_v2.0.fits'   
             }

After inserting your files QtCounterpart is ready to run. 

**Running QtCounterpart**

**What do I see?**
The lower left big panel shows the narrowband of the objects detected by MUSE (MUSE narroband), those saved in the MUSE catalogue. The top left big panel shows a HST band cutout (main HST image) corresponding to the same MUSE region where the detected object is located. Similarly, the four top small panels (HST band 435, 606, 775/814, and 160) show different HST bands of the same HST image cutout. Inside the main HST image, the black circle with the inner cross represents the MUSE position of the object, while the empty white circles show the position of the possible HST counterparts of the MUSE source. The numerical boxes below the small HST panels show the current resolution cut of the HST images, representing the scaling factor used for the intensity of the HST images. The "Reset images" button offers the possibility to set all GUI images back to the original display. Spectroscopic data from MUSE is displayed in the left half of the GUI. There, information such as MUSE ID, right ascension (RA), declination (Dec), redshift (z), confidence of the source detection, the line that lead to the detection or the signal to noise (S/N) of the detection is found. Right above this information, photometric data from the phtotometric catalogue is shown. The photometric UV ID, the difference in RA between the spectroscopic and the photometric values, same for Dec, the photometric redshift, and the confidence of the detection are displayed. 

![GUI screenshot](https://user-images.githubusercontent.com/81447012/154290161-9691d87a-1c94-4fda-974c-9a7f47626531.PNG)

**What can I do?**

**General overview** (inside what can I do section)
In the images you can compare the source detections in MUSE with the corresponding ones observed in different HST bands. Besides, the possible counterparts to the MUSE object are marked with a white circle in the main HST image, based on spatial assumptions. You can select the most likely counterpart by clicking on those circles and click on the confidence buttons to record how sure you were about the decision. Then, the information of that counterpart will be displayed below the photometric information section. By clicking on the "next" button you will save all the displayed information in a fits output Table. 

**Script overview** (inside what can I do section)

There are several scripts that the main GUI script, ident_counterpart.py, uses:
colormaps.py : New matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt, and (in the case of viridis) Eric Firing.
colors.py: get colors for the plots from colormaps.py
read_write_files.py: read, saves and modifies fits tables by creating dictionaries with information using the first column as identifiers. 

As already mentioned, the main script is ident_counterpart.py. Below, I explain the purpose of each of button and the general functionalities of the GUI:
The contrast of the small HST images can be modified by writing a different cut number on the cut boxes and clicking afterwards on the "Cuts" button. The smaller that number is, the lower the resolution of the images. By using the "Reset images" button, the resolution of all HST images is back to its default (0.99).
Now, if one is more interested in one spectific HST band than in the rest of them, by double clicking on the wanted HST band, you will convert it into the main HST image of the GUI.
The main HST image already hints at possible counterparts for the currently displayed MUSE source. You can select any of the white circles to choose the most likely counterpart. This circle will then turn into red and the photometric information of the selected object will be displayed below the Photometric information section. In order to guide the eye in the counterpart selection, the main HST image also contains the position of the MUSE object according to the spectroscopic data (shown below the MUSE information section). Once you have selected the counterpart, you must tell the GUI how sure you are about your selection. For that, you can click on the "Confidence" bullet points, indicating "0" for the most unsure case and "3" for the most sure scenario. Once you have done this, you have found an HST counterpart for the MUSE source so you can click on "Next" to save the displayed information and find the counterpart for the next object in the MUSE catalogue. Contrarily, if you believe you made a mistake selecting the counterpart, you only need to click on the "Clear counterparts" button to start classifying that MUSE source from the beginning.
There might be of course cases for which the MUSE object has no counterpart in the HST data, one because faint or bc doesnt exist


The most important functions in the script are described below, but for further details, refer to the header of each routine.
...


**QtCounterpart output**

The output file generated by QtCounterpart is a fits Table containing all the counterpart identification results. The table contains the following columns:

**Tips for users**

When opening QtCounterpart, patience is needed. The datacubes are huge and take some time to load and extract the narrowband out of it. You need to make sure though that you used the right input and output catalogue file names, with the required column names and information!

Comments, issues and/or suggestions for improvement, new functions etc. are more than welcome. Please send these to Yohana Herrero Alonso (yherreroalonso at aip.de) or add an 'issue' via GitHub.
