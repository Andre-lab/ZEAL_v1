# User guide

## Loading structures
Import the two structures you want to align in the *import-tab*: the **fixed" structure** is the reference structure and the **rotating" structure** is the structure whose (rotational) orientation is going to be changed after the alignment. ZEAL allows structures to be loaded from pdb files, or to be downloaded directly from the PDB server - use the load-fetch switch to change between the modes. If a particular chain should be used, enter this before clicking import.

## Shape-alignment 
In the *alignment tab* there are three modes available for generating a shape alignment:

| Button | Description |
|--|--|
| <img src="../images/ZEAL/icons/the_flash_100px.png" height="30px"> | **Global search**  Automatic search for the best alignment using a machine learning algorithm called [surrogate optimization](howItWorks.md). The search progress is shown in a new window, showing the best correlation found for each function evaluation. The search can be stopped by clicking stop, otherwise the search will proceed until the stopping criteria is met (default is 300 function evaluations). Clicking on <img src="../images/ZEAL/icons/settings_100px.png" height="15px"> opens a settings window where the stopping criteria can be changed, and the number of random samples used in the surrogate creation phase.  |
| <img src="../images/ZEAL/icons/gyro_100px.png" height="30px"> | **Manual alignment** Interactive alignment by using the mouse to change the orientation of the rotating structure. The new orientation will be the starting point for the Global-search algorithm. |
| <img src="../images/ZEAL/icons/dice_filled_60px.png" height="30px"> | **Randomize** Apply a random rotation to the rotating structure. The maximmum change (in Euler angles) can be set by clicking on <img src="../images/ZEAL/icons/settings_100px.png" height="15px"> which will open a settings window. |


After selecting a mode, clicking on <img src="../images/ZEAL/icons/next_96px.png" height="20px"> will activate it. 

### Region of interest
Sometimes a structure contains a region that you want to ignore for the shape alignment. For such cases, ZEAL allows interactive selection of a region-of-interest (ROI) using JSmol. Defining ROIs for each structure is done in the "Define Region of Interest" tab. 

As an ***example***, let's consider the fixed structure ***6fi8-H*** and the rotating structure ***5mn2-D***. This is how the shape alignment looks like after running the search

<img src="../images/ZEAL/ROI/ZealRF_noROI_duoWin.png"/>

The fixed structure has an alpha-helical tail which we wish do leave out to get a different alignment. We do the following steps to do it:

<img src="../images/ZEAL/ROI/ZealRF_ROIsetup.png" class="callout" height="200px"/>

1. CLick on the *Define Region of Interest* tab and select the fixed structure by clicking on ***fixed***. 

2. Click on ***define*** and a new window with JSmol will be shown. 

3. Define your ROI using the mouse: ***click-and-drag while holding the left SHIFT key*** to select/deselect atoms within the selection box. Holding left SHIFT+ALT (SHIFT+OPTION in Mac) will deselect all atoms within the selection-box. Selected atoms will have yellow halos around them. You can invert the selection by clicking on invert. When you are done, click ***Use selection*** in the JSmol window and you will see how the resulting structure is defined in ZEAL, with the non-ROI colored in white. 
<img src="../images/ZEAL/ROI/ZealRF_ROIsel.png"/>

4. Enable ROI-based shape alignment by clicking on the ROI-checkbox in the *Setup tab* and then run the search. 
<img src="../images/ZEAL/ROI/ZealRF_ROI_import.png"/>


This is how the new shape alignment looks like when running the search with ROI-constraint.

<img src="../images/ZEAL/ROI/ZealRF_ROI_duoWin.png"/>


## Graphics and visualizations
The *visualization tab* contains options that control the apperance and visualization of the structures, including the open source molecular viewer JSmol.

<img src="../images/ZEAL/win/ZealRF_vizTab.png"/>

| Button | Description|
|--|--|
| <img src="../images/ZEAL/icons/settings_100px.png" height="20px">  |Opens a new window with settings controlling scene lights and surface smoothing for shapes shown in ZEAL. Some parameters for JSmol can also be changed here.  |
| <img src="../images/ZEAL/icons/paint_palette_100px.png" height="20px"> | Opens a color-picking window to set another color for the structure. This will also be the default color used by JSmol. |
| <img src="../images/ZEAL/icons/hide_100px.png" height="20px"> |	Hides the structure from view. |
|  <img src="../images/ZEAL/icons/opacity_100px.png" height="20px">	 | Makes the surface transparent. |


### JSmol
ZEAL has JSmol with webGL built in to allow more sophisticated visualizations of the structures after performing a shape alignment. JSmol can be started by clicking on the JSmol-logo in the visualization tab.   

## Settings 
The *Setup* tab contains a number of options that effect the ZC-decomposition of the structures and thus the resulting shape alignment.

### Zernike-Canterakis moments
For computational reasons, only ZC moments of maximum order 1-20, 25 and 30 can be selected. However, order 20 has been shown to capture enough salient features of 3D shapes. See [shape reconstructions](howItWorks.md) to see how much shape information a set of ZC moments contain using the available moment orders. 

### Shape representation
There are three representations to choose frome:

1. The shell defined by the solvent-excluded or solvent-accessible surface. The probe-radius can be changed as well as 
the thickness (in grid units) of the resulting "surface shell".
2. A density analogous to the electron density. The atoms are "smeared" out on the grid, and the relative size of the radii of the smeared atoms is controlled by the "Smear factor" (the fraction of the grid to smear out over). The iso-surface to be shown is set by the relative iso-value (0-1). 




