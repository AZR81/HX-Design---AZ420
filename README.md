# HX-Design---AZ420
The code used for the heat exchanger design project.
The prerequisite libraries are Numpy, Matplotlib and open-cv.

The PNG images were taken from Coulson and Richardson's book.  
Each file performs the following tasks:

iterator.py:  
-Performs all design calculations.
-Saves the results of the calculations in a text file in the Output folder for future processing. The folder is needed for this file to run.
  
constants.py:  
-Stores all the data and the functions needed for processing it.
  
image_to_poly.py  
-Loads grayscale images and fits a polynomial or returns the scaled data so it can be used for linear interpolation.
  
mechanical_design.py  
-Performs calculations for the mechanical design section.

log_analyser.py  
-Re-plots graphs. Useful for changing the font and colours without having to recalculate the data.  
-Converts the chat output of the iterator file into a more understandable format.

page_scraper.py  
-Loads the heat exchanger specifications from a set of given URLs (model_links.txt). Data is taken from https://www.shell-tube.com/Manufacturers/americanindustrial.html. Saves the data to model_output.txt.

data_processor_1.py  
-Gives a list of materials used by the heat exchangers found in model_output.txt.
