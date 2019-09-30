# . <img src="./images/zealLogo.png" height="100px"> 

## Shape-based molecular alignment

---

## Overview
ZEAL is an app that allows two molecules to be aligned based on their shape resemblance alone. This is done by projecting the shape onto 3D Zernike-Canterakis functions up to a specific order; the higher the order, the higher is the level of detail captured by the moments associated with the functions. ZEAL will find the optimalshape alignment by searching for the rotation that maximizes the correlation of these moments. 


## Installation

1. Download the ZEAL installer from [github.com/Andre-lab/ZEAL](https://github.com/Andre-lab/ZEAL/tree/master/installation/ZEALstandalone).
   The installer will download the [MATLAB-runtime](https://se.mathworks.com/products/compiler/matlab-runtime.html) enviroment necessary for ZEAL to run. 

2. To install the standalone application, double-click the **ZEALinstaller_web.app** executable.
	
3. To complete installation, follow the instructions on the user interface.

	??? note "Proxy server"
		If you want to connect to the Internet using a proxy server, click Connection Settings. Enter the proxy server settings in the provided window. Click OK.
	

	!!! Note
		On Linux and Mac OS X, you do not have the option of adding a desktop shortcut.




### MatlabApp
If you have MATLAB installed, you can install ZEAL as a [Matlab app](https://se.mathworks.com/discovery/matlab-apps.html) directly. 

1. Download the mlapp-installer file ´ZEAL.mlappinstall´  at [github.com/Andre-lab/ZEAL/](https://github.com/Andre-lab/ZEAL/tree/master/installation/MatlabApp)
2. Double-click the **ZEAL.mlappinstall** executable and ZEAL will be installed as an app within MATLAB.

Start ZEAL via *The Apps tab* of the MATLAB Toolstrip, which shows you the apps that you currently have installed.

## Launching the app
To run your standalone application:

1. Open a terminal window.

	??? summary "Windows"
				1. Click the "Start" button to open the Start menu.
				2. Open the "All Programs" menu, followed by the "Accessories" option.
				3. Select the "Command Prompt" option from the "Accessories" menu.

	??? summary "Mac OS X"		
				1. Open a Finder window.
				2. Select Applications from the left side.
				3. Click the arrow to expand the Utilities folder.
				4. Double-click Terminal.

2. Navigate to the folder into which you installed the application. If you accepted the default settings, you can find the folder in one of the following locations by executing the command in the terminal:

	|  Platform                       | Path                         | 
	| ----------------------------- |:------------------------------:| 
	| Windows                       | `cd C:\Program Files\ZEAL`         | 
	| Mac OS X                      | `cd /Applications/ZEAL`             | 
	| Linux                         | `cd /usr/ZEAL`                      | 

3. Run the application using one of the following commands:

	|  Platform | Path                           | 
	| ----------|:------------------------------:| 
	| Windows   | ```application\ZEAL  ```        | 
	| Mac OS X  | ``` ./application/run_ZEAL.sh /Applications/MATLAB/MATLAB_Runtime/v97 ``` |
	| Linux     | ``` ./ZEAL ```                     | 

	??? note "Mac OS X"
				- The shell script `run_ZEAL.sh` will temporarily set the necessary environment variables and start 
			 	the application)

			   	- to run the shell script, type
			   
			       		./run_ZEAL.sh <mcr_directory> <argument_list>
			       
			    at Linux or Mac command prompt. `<mcr_directory>` is the directory 
			    where version 9.7 of the MATLAB Runtime is installed or the directory where 
			    MATLAB is installed on the machine. 

			    For example, if you have version 9.7 of the MATLAB Runtime installed in 
			    /mathworks/home/application/v97, run the shell script as:
			    
			       ./run_ZEAL.sh /mathworks/home/application/v97
			       
			    If you have MATLAB installed in /mathworks/devel/application/matlab, 
			    run the shell script as:
		    
		       ./run_ZEAL.sh /mathworks/devel/application/matlab

## Getting started
Once ZEAL is running, performing a shape alignment is easy and the app design should make most operations feel intuitive. 

1. Import the two structures you want to align by clicking the green load buttons in the import-tab. The "fixed" structure is the reference structure and 	the "rotating" structure is the structure whose orientation is going to be changed after the alignment. ZEAL allows structures to be loaded from pdb files, or to be downloaded directly from the PDB server - use the load-fetch switch to change between the modes. <img src="./images/ZEAL/processed/ZEAL_start.png" height="500px"> 

2. After importing the two structures you will see a low-resolution representation of the solvent-accessible surface of the structures - this can be changed in the setup tab. Note that the structures might have different relative sizes due to scaling (this is to achieve [scale invariance](howItWorks.md) when comparing the shapes). <img src="./images/ZEAL/processed/ZEAL_loaded.png" height="500px"> 

3. Click the *Play button* (*Global search*) to start the search for optimal shape alignment. The progress of the search will be shown in a new window; the alignment with the best correlation coefficient found after a certain number of function evaluations. The search can be stopped by clicking on the stop button, otherwise it will stop when the [stopping criteria](userGuide.md) is met. 
<img src="./images/ZEAL/processed/ZEAL_Search.png" height="500px"> 

4. To inspect the alignment in more detail you can view the structures in JSmol if you like. Click the export tab if you want to save the aligned structures to a pdb file. 
<img src="./images/ZEAL/processed/ZEAL_JSmol.png" height="500px"> 
