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
The main HST image already hints at possible counterparts (only when existing in the photometric catalogue) for the currently displayed MUSE source. You can select any of the white circles to choose the most likely counterpart. This circle will then turn into red and the photometric information of the selected object will be displayed below the Photometric information section. In order to guide the eye in the counterpart selection, the main HST image also contains the position of the MUSE object according to the spectroscopic data (shown below the MUSE information section). Once you have selected the counterpart, you must tell the GUI how sure you are about your selection. For that, you can click on the "Confidence" bullet points, indicating "0" for the most unsure case and "3" for the most sure scenario. Once you have done this, you have found an HST counterpart for the MUSE source so you can click on "Next" to save the displayed information and find the counterpart for the next object in the MUSE catalogue. 
There might be of course cases for which the MUSE object has no counterpart in the HST data, the object is not registered in the photometric catalogues or one MUSE object has more than one counterpart. These cases can be recorded by clicking in one of the empty boxes below the "Confidence" circles. You can also write a comment so it is easier to recall the case in future inspections. If you believe you made a mistake selecting the counterpart or any of these possible scenarios, you only need to click on the "Clear counterparts" button to start classifying that MUSE source from scratch.
In the case that the most likely counterpart does not appear in the photometric catalogue, there will be no white circle around that galaxy in the main HST image. However, it is possible to perform a forced photemetry by cliking (in the main HST image) on the center of the galaxy that you believe is the counterpart but has not been recorded in the photometric catalogues. Once you have clicked on the area of the main HST image where that undetected counterpart is, the GUI measures the flux and coordinates of the object. Specify the case that you encountered ("Counterpart not in catalog"), leave a comment if you wish and click "Next" to save the information.
Any time it is possible to go back to the already classified objects by clicking on "Previous" until you find the desired object. Remember that if you want to save or rewrite the classification of an object, you must click on "Next" after you have clicked on any button/counterpart of the GUI. A more practical way of doing this is to write the MUSE ID of the desired object in the empty box right the "Jump to MUSE ID" button, click afterward on the button and your desired object will be again shown in the display.
Finally, let's say you do not feel like classifying more objects. Then, press the red cross on the top right corner of the display to exit the GUI. Once you want to classify more objects again, you only need to run the GUI and click on the "Jump to next unclassified ID" button to continue with the object where you left it.

The most important functions in the script are described below, but for further details, refer to the header of each routine.
clickable(widget): detects double click on the small HST bands
get_HST_images(): extract HST and MUSE data and headers from the input fits files
get_infos(): create MUSE NB from the MUSE data cube and extract spectroscopic information to be displayed under "MUSE information"
get_pos_UV(): extract photometric information to be displayed under "Photometric information" and the locations of the possible HST counterparts
image_Widgeti(), main_image_Widget(), image_Widget_narrowband_muse(): cut the relevant area of the small HST band i and display it. Same for the main HST image or the MUSE NB image
small_hst_resolution(): change the resolution of the HST bands
make_ismallHST_mainHST(): convert the small HST band i into the main HST image
clicked_counterpart(): display information of the clicked counterpart on the main HST image
clickToObtainPositionOnImage(): allow to click on main HST image to extract the position of clicked area (used in the forced photometry when the MUSE object is not recorded in the photometric catalogues)
fill_window(): display buttons and images on the screen
show_next(): called when clicked on "Next". Save displayed information and perform forced photometry (if needed)


**QtCounterpart output**

The output file generated by QtCounterpart is a fits Table containing all the counterpart identification results. The table contains the following columns:
'ID': ID of the MUSE object
'UV_ID': photometric ID of the counterpart
'MUSE_RA': spectroscopic RA from MUSE (in deg)
'UV_RA': photometric RA of the counterpart (in deg)
'delta_RA': difference between MUSE_RA and UV_RA
'MUSE_DEC': spectroscopic DED from MUSE (in deg)
'UV_DEC': photometric DEC of the counterpart (in deg)
'delta_DEC': difference between MUSE_DEC and UV_DEC
'MUSE_z': spectroscopic redshift
'UV_z': photometric redshift of the counterpart
'Separation': spatial separation between the spectroscopic object and photometric location of the counterpart (np.sqrt((MUSE_RA-UV_RA)**2+(MUSE_Dec-UV_DEC)**2))
'Confidence': confidence of our selection
'Comment': comment written during the selection
'More than one counterpart': 'Yes' if there is more than one counterpart and '-' if only one or none
'No match': 'Yes' if there is no match between the spectroscopic and photometric catalogue for that object and '-' if there is no or one counterpart 
'Not in catalog': 'Yes' if the counterpart is not in the photometric catalogue and '-' if it is or if there is no counterpart
'Photometry': Forced photometry in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image (in erg/s/cm^2/Angstrom)
'dec_noMatch': DEC of the undetected counterpart in photometric catalogues, in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image
'ra_noMatch': RA of the undetected counterpart in photometric catalogues, in case 'Not in catalog' = 'Yes' and one has clicked somewhere on the main HST image

**Tips for users**

When opening QtCounterpart, patience is needed. The datacubes are huge and take some time to load and extract the narrowband out of it. You need to make sure though that you used the right input and output catalogue file names, with the required column names and information!

Comments, issues and/or suggestions for improvement, new functions etc. are more than welcome. Please send these to Yohana Herrero Alonso (yherreroalonso at aip.de) or add an 'issue' via GitHub.
