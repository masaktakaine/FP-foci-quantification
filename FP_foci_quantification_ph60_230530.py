# FP foci quantification, jython version
# developed by Masak Takaine

# This FIJI macro allows automatic detection and analysis of intracellular fluorescent foci in yeast cells.
# "Foci" are defined as intracellular fluorescent structures that have stronger intensitie than fluorescent granules.
# As input, two channel image file that contains a pair of fluorescence (Ch #1) and phase contrast (Ch #2) microscopic images is required.
# This macro is optimized for images accuired by using a 60x phase-contrast objective lens.

# Cell outlines are extracted from phase-contrast images, fluorescent foci inside the cell boundary are detected using the FindMaxima function.
# You can specify "prominence" (former "Noise Tolerance") of the FindMaxima function.

from __future__ import division
from ij import IJ, ImagePlus, Prefs
from ij.process import ByteProcessor
from ij.process import ImageStatistics as IS
options = IS.MEAN | IS.AREA | IS.STD_DEV  # many others
from ij.gui import Roi

from ij.measure import ResultsTable
from ij.plugin import ImageCalculator
from ij.plugin import ChannelSplitter as CS
from ij.plugin import RGBStackMerge, RGBStackConverter
from ij.plugin import Duplicator
from ij.plugin import ZProjector as ZP
from ij.plugin.frame import RoiManager as RM
from ij.plugin.filter import GaussianBlur
from ij.plugin.filter import MaximumFinder
from ij.plugin.filter import BackgroundSubtracter as BS
createBackground = False
lightBackground = False
useParaboloid = False
doPresmooth = False
correctCorners = False

import os
from os import path
import random

#@ String(label="Date of experiments, e.g., 2022-02-05") edate1
#@ int(label="Prominence for FindMaxima") prom
#@ File (label="Choose source Folder", style="directory") dirS0 
#@ File (label="Choose destination Folder", style="directory") dirD0

file_name = [] 	# An array to store filenames
total_cell_number = []	# An array to store the totall cell number counted in a file
foci_cell = [] 	# An array to store the number of cells having fluorescent foci 
pct_foci_cell = [] 	# An array to store the percent of cells having fluorescent foci
prominence = []
defocus = []

sno = 0 # A variable to count observation： Serial No. of observation, equals total number of rows in table "foci_data"
snf = 0 # A variable to count the total number of fluorescent foci: Serial No. of Foci

date = []			#An array to store date of experiments
image_file = []  	#An array to store filenames
cell_number =[]		#An array to store the serial number of cell
observation = []		#An array to store the serial number of observation
foci_serial = []		#An array to store the serial number of fluorescent granule
foci_meanints =[]		#An array to store mean fluorescence intensity of the granule
foci_x = []			#An array to store the x coordinate of the granule
foci_y = []			#An array to store the y coordinate of the granule
cell_area = []		#An array to store the area of cell
prominences = []		#An array to store the current prominence
withfoci = []		#An array to check existance of granule in the cell

# Save the Result Table in csv format
def save_result_table(directory, filename, result_table):
    resultfile = os.path.join(directory, filename + ".csv") 
    result_table.saveAs(resultfile)

# Save the image file in tiff format
def save_image_as_tif(directory, filename, image):
    outputfile = os.path.join(directory, filename + ".tif")
    IJ.saveAs(image, "TIFF", outputfile) # 保存する画像、tiff形式、パスを指定

def FP_foci_quantification_findmaxima(current_file_path):
	global edate, sno, snf
	imp0 = IJ.openImage(current_file_path)
	width = imp0.getWidth()
	height = imp0.getHeight()
	filename = imp0.getTitle().split(".")[0]
	
	# Max intensity projection
	imp0_mip = ZP.run(imp0, "max all")
	
	cs = CS()
	image_list = cs.split(imp0_mip)
	ch1 = image_list[0] # Channel #1：fluorescence image
	ch2 = image_list[1] # Channel #2：phase-contrast image
	
	ch2_mask = ch2.duplicate()
	ch2_mask.setTitle("mask")
	
	# Gaussian filter
	gb = GaussianBlur()
	sigma = 1
	accuracy = 0.01
	gb.blurGaussian(ch2_mask.getProcessor(), sigma, sigma, accuracy)  # x & y sigmas
	
	# BackGround subtraction, 
	radius = 25
	bs = BS()
	bs.rollingBallBackground(ch2_mask.getProcessor(), radius, createBackground, lightBackground, useParaboloid, doPresmooth, correctCorners)
	
	IJ.setAutoThreshold(ch2_mask, "Li dark")
	Prefs.blackBackground = True
	IJ.run(ch2_mask, "Convert to Mask", "")
	IJ.run(ch2_mask, "Dilate", "")
	IJ.run(ch2_mask, "Skeletonize", "")
	IJ.run(ch2_mask, "Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune_0")
	IJ.selectWindow("Tagged skeleton")
	ts = IJ.getImage()
	ts.hide()
	IJ.run(ts, "Convert to Mask", "")
	IJ.run(ts, "Dilate", "")
	IJ.run(ts, "Close-", "")
	
	ts2 = ts.duplicate()
	IJ.run(ts2, "Fill Holes", "")
	ts2.hide()
	
	cell_ol = ImageCalculator().run("Subtract create", ts2, ts)
	#imp3.show()
	IJ.run("Set Measurements...", "area redirect=None decimal=3")
	IJ.run(cell_ol, "Analyze Particles...", "size=400-3000 pixel show=Outlines display exclude clear add")  # size unit is µm
	IJ.selectWindow("Drawing of Result of DUP_Tagged skeleton")
	drawing = IJ.getImage()
	drawing.hide()
	rm = RM.getRoiManager()
	nRois = rm.getCount() # ROIの総数を取得
	cwf = 0 # A variable to count the number of cells showing fluorescent foci: Cell With Foci
	
	rt = ResultsTable.getResultsTable()
	rt.reset()
	for k in range(0,nRois):
	    rt.reset()
	    
	    # Assign the k-th ROI on ch1
	    rm.select(ch1, k)
	    area = IJ.getValue(ch1, "Area")
	    MaximumFinder().findMaxima(ch1.getProcessor(), prom, MaximumFinder.COUNT, True)
	    nf_k = rt.getValue("Count", 0)
	    rt.reset()
	    # If the k-th cell showed some foci
	    if nf_k >= 1:
			cwf = cwf + 1
			MaximumFinder().findMaxima(ch1.getProcessor(), prom, MaximumFinder.LIST, True)
			for j in range(0, int(nf_k)):
				observation.append(sno+1)
				date.append(edate)
				foci_serial.append(snf+1)
				sno = sno + 1
				snf = snf + 1
				cell_number.append(k+1)
				image_file.append(filename)	
				x = rt.getValue("X", j)
				y = rt.getValue("Y", j)
				foci_x.append(x)
				foci_y.append(y)
				cell_area.append(area)
				prominences.append(prom)
				IJ.run(ch1, "Specify...", "width=2 height=2 x="+str(x)+" y="+str(y)+" oval centered");
				foci_meanints.append(IJ.getValue(ch1, "Mean"))
				withfoci.append("TRUE")
				
		# nf_k == 0, i.e., if the k-th cell showed no foci			    	
	    else: 
			observation.append(sno+1)
			sno = sno + 1
			date.append(edate)
			foci_serial.append("NaN") 			
			cell_number.append(k+1)
			image_file.append(filename)	
			foci_x.append("NaN")
			foci_y.append("NaN")		
			cell_area.append(area)
			prominences.append(prom)
			foci_meanints.append("NaN")
			withfoci.append("FALSE")
			
	total_cell_number = nRois
	foci_cell = cwf
	pct_foci_cell = 100*cwf/nRois
	prominence = prom
	
	# Create a merge image
	IJ.run(ch1, "Grays", "")
	merge_cm = RGBStackMerge.mergeChannels([None, None, None, None, ch2, ch1, None], True)
	RGBStackConverter.convertToRGB(merge_cm)
	
	# Create a mask image of foci detected
	bp = ByteProcessor(width, height)            # create ImageProcessor
	title_foci_mask = str(filename) + "_foci_location"
	foci_mask = ImagePlus(title_foci_mask, bp)    # create ImagePlus with specific title and ImageProcessor object
	for x, y in zip(foci_x, foci_y):
		if not x == "NaN":
			bp.setColor(255)
			bp.fillOval(int(x)-2, int(y)-2, 4, 4)
	
	# Collect parameters or lists of parameters into a dictionary
	params = {"date":date,"image_file":image_file, "observation":observation, "cell_number":cell_number, "withfoci":withfoci, "foci_serial":foci_serial,\
	"foci_meanints":foci_meanints, "foci_x":foci_x, "foci_y":foci_y, "cell_area":cell_area, "prominences":prominences}
	params2 = {"filename":filename, "total_cell_number":total_cell_number, "foci_cell":foci_cell, "pct_foci_cell":pct_foci_cell, "prominence":prominence}
	
	return ch1, ch2, drawing, merge_cm, params, params2, foci_mask

#### Main code
# Insert a blank to prevent automatic modification on Excel.
edate = " "+edate1
# Make directories
dirD = os.path.join(str(dirD0), edate1 + "_output")
if not os.path.exists(dirD):
	os.mkdir(dirD)
dirBF = os.path.join(str(dirD), "BF")
if not os.path.exists(dirBF):
	os.mkdir(dirBF)                           
#Create a folder for mask images and ROI data
dirDR = os.path.join(str(dirD), "Drawings")
if not os.path.exists(dirDR):
	os.mkdir(dirDR)
dirGreen = os.path.join(str(dirD), "Green")
if not os.path.exists(dirGreen):
	os.mkdir(dirGreen)
dirMerge = os.path.join(str(dirD), "merge")
if not os.path.exists(dirMerge):
	os.mkdir(dirMerge)

# Acquire a list of files in the directory
filelist = os.listdir(str(dirS0))

# List comprehension, extract nd2 files.
nd2_files = [f for f in filelist if f.split(".")[-1] == "nd2"]
#filenames = [f.split(".")[0] for f in filelist]
nd2_files = sorted(nd2_files)

# Create a table that summarises the averages of particle parameters in a image file
foci_stat = ResultsTable()

for nd2_file in reversed(nd2_files):  # reversed() generates a revered iterator
    current_file_path = os.path.join(str(dirS0), nd2_file) 
    results = FP_foci_quantification_findmaxima(str(current_file_path))
      
    params = results[4]
    params2 = results[5]
    filename = params2["filename"]
    
    save_image_as_tif(str(dirGreen), filename, results[0])
    save_image_as_tif(str(dirBF), filename, results[1])
    save_image_as_tif(str(dirDR), filename, results[2])
    save_image_as_tif(str(dirMerge), filename, results[3])
    save_image_as_tif(str(dirDR), results[6].getTitle(), results[6])
    
    foci_stat.addValue("date", edate)
    foci_stat.addValue("file_name", filename)
    foci_stat.addValue("total_cell_number", params2["total_cell_number"])
    foci_stat.addValue("foci_cell", params2["foci_cell"])
    foci_stat.addValue("pct_foci_cell", params2["pct_foci_cell"])
    foci_stat.addValue("prominence", params2["prominence"])
    foci_stat.addRow()

#　Remove the last empty row
foci_stat.deleteRow(foci_stat.size() - 1)
save_result_table(str(dirD), edate1+"_foci_stat", foci_stat)

# Create a table that summarises all observations
foci_data = ResultsTable(sno) # Construct a table with "sno" rows

# Extract each list from the dictionary, and add the values as a column in the table
for k, v in params.items():
		for j in range(0, len(v)):
			foci_data.setValue(str(k), j, v[j])
		
save_result_table(str(dirD), edate1+"_foci_data", foci_data)
#print "date", len(params["date"])
#print "image_file", len(params["image_file"])
#print "observation", len(params["observation"])
#print "cell_number", len(params["cell_number"])
#print "withfoci", len(params["withfoci"])
#print "foci_serial", len(params["foci_serial"])
#print "foci_meanints", len(params["foci_meanints"])
#print "foci_x", len(params["foci_x"])
#print "sno", sno
#print "snf", snf

print "Done. \n"
IJ.run("Clear Results")
rm = RM.getRoiManager()
rm.reset()
IJ.run("Close All")