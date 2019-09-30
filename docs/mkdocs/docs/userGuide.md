# User guide

## Loading structures
Import the two structures you want to align by clicking the green load buttons in the import-tab. The "fixed" structure is the reference structure and 	the "rotating" structure is the structure whose rotation is going to be changed after the alignment. ZEAL allows structures to be loaded from pdb files, or to be downloaded directly from the PDB server - use the load-fetch switch to change between the modes.

## Shape-alignment 
In the alignment tab there are three "strategies" available for generating a shape alignment:

  > <img src="../images/ZEAL/icons/the_flash_100px.png" height="30px"> Global search

  > <img src="../images/ZEAL/icons/gyro_100px.png" height="30px"> Manual alignment

  > <img src="../images/ZEAL/icons/dice_filled_60px.png" height="30px"> Randomize

Clicking on one of them will select the given strategy, and clicking on the <span style="color:green"> &#9654; </span> (*Play*) button will apply it.
In summary, the Global search option will automatically find the optimal alignment, wherease the other can be used to generate starting orientations for search algorithm. 

### Global search
This will automatically search for the best alignment. The search progress is shown in a new window, showing the best correlation found for each function evaluation. The search can be stopped by clicking on the stop button, otherwise the search will proceed until a stopping criteria is met; these can be changed in the settings - click on the &#9881; (gear) button. 

#### Region of interest
Sometimes a structure contains a region that you want to ignore for the shape alignment. For such cases, ZEAL allows interactive selection of a region-of-interest (ROI) using JSmol. Defining ROIs for each structure is done in the "Define Region of Interest" tab. Follow these steps to enable shape alignment with a ROI:

1. Select a structure by clicking on "rotating/fixed". 

2. Click on "Select". You will first be prompted to save the structure to a pdb-file, and then prompted to load it in JSmol. Once the structure is loaded, click on the button to setup JSmol selection. (In the case that you do not see a structure in JSmol after loading it, please redo close the JSmol window and redo this step). 

3. Select your ROI using the mouse: click-and-drag while holding the left SHIFT key to define select/deselect atoms within the selection box. Holding SHIFT+ALT (SHIFT+OPTION in Mac) will deselect all atoms within the selection-box. Selected atoms will have yellow halos around them. You can invert the selection by clicking on invert. When you are done, click save and you will be prompted to save a file. 

4. Import the ROI by clicking on import. You will prompted to load the file you saved in step (3).

5. Enable ROI-based shape alignment by clicking on the relevent checkbox in the Setup tab. 

##### Example


### Manual alignment
This will allow interactive alignment by using the mouse to change the orientation of the rotating structure.

### Randomize 
The will apply a random rotation to the rotating structure. The amount can be changed in the settings - click on the &#9881; (gear) button. 

## Settings


### Zernike-Canterakis moments

### Shape representation
The protein-shape representation can be defined in the setup tab. There are three representations to choose frome:

1. The shell defined by the solvent-excluded or solvent-accessible surface. The thickness (in grid units) of these shells can be changed as well.
2. A density analogous to the electron density. The atoms are "smeared" out on the grid, and the relative size of the smeared atoms is controlled by the "Smear factor" (DEFINITION!). The iso-surface to be shown is set by the relative iso-value (0-1)

### Graphics
Settings controlling scene lights and surface smoothing can be controlled in the graphics settings by clicking on 🌄 (*Sunrise Over Mountains*) in the "View structures" tab. 

#### Light

#### Smoothing 

#### Colors 








