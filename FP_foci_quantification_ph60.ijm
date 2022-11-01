// FP foci quantification
// author: Masak Takaine

// This FIJI macro allows automatic detection and analysis of intracellular fluorescent foci in yeast cells.
// "Foci" are defined as intracellular fluorescent structures that have stronger intensitie than fluorescent granules.
// As input, two channel image file that contains a pair of fluorescence (Ch #1) and phase contrast (Ch #2) microscopic images is required.
// This macro is optimized for images accuired by using a 60x phase-contrast objective lens.

// Cell outlines are extracted from phase-contrast images, fluorescent foci inside the cell boundary are detected using the FindMaxima function.
// You can specify "prominence" (former "Noise Tolerance") of the FindMaxima function.

macro "FP_foci_quantification" {

#@ String(label="Date of experiments, e.g., 2022-02-05") edate1
// Input a prominence(= noise tolerance) used in the FindMaxima function
#@ int(label="Prominence for FindMaxima") prom 
#@ File (label="Choose source Folder", style="directory") dirS0
#@ File (label="Choose destination Folder", style="directory") dirD0
#@ String(label="Hide/Show the active image? The Show slows the analysis.", choices={"hide","show"}, style="radioButtonHorizontal") sbm

setBatchMode(sbm); // hides the active image, required ImageJ 1.48h or later
dirS = dirS0 + File.separator; // "File.separator" returns the file name separator character depending on the OS system used.
dirD1 = dirD0 + File.separator;

//setOption("ExpandableArrays", true);  // Enables/disables support for auto-expanding arrays, In ImageJ 1.53g or newer, arrays automatically expand in size as needed.
edate = " "+edate1; // Insert a blank to prevent automatic modification on Excel.

imagefilelist = getFileList(dirS);
file_name = newArray; 	// An array to store filenames
total_cell_number = newArray;	// An array to store the totall cell number counted in a file
foci_cell = newArray; 	// An array to store the number of cells having fluorescent foci 
pct_foci_cell = newArray; 	// An array to store the percent of cells having fluorescent foci
prominence = newArray;
defocus = newArray;  // An array to check images that are not suitable for analysis

dirD = dirD1 +  "/"+ "prom_" + prom;
File.makeDirectory(dirD);

dirBF = dirD + "/BF/";				//Create a folder for phase-contrast or biright-field images
File.makeDirectory(dirBF);				//

dirGreen = dirD + "/green/";				//Create a folder for fluorescent images
File.makeDirectory(dirGreen);				//

dirMer = dirD + "/merge/";				//Create a folder for a merge image (fluorescent + phase-contrast)
File.makeDirectory(dirMer);	

dirDR = dirD + "/Drawings/";				//Create a folder for mask images and ROI data
File.makeDirectory(dirDR);			//

dirCSV = dirD +"/" + edate1 + "_prom" + prom + "_csv/";
File.makeDirectory(dirCSV);		// Create a folder for csv files

for (i = 0; i < imagefilelist.length; i++) {
  currFile = dirS + imagefilelist[i];
    if((endsWith(currFile, ".nd2"))||(endsWith(currFile, ".oib"))||(endsWith(currFile, ".zvi"))) { // process if files ending with .oib or .nd2, or .zvi
		run("Bio-Formats Macro Extensions"); 
		Ext.openImagePlus(currFile)}
	else if ((endsWith(currFile, ".tif"))||(endsWith(currFile, ".tiff"))) {// process if files ending with .tif or .tiff (hyperstack files)
			open(currFile); 
		}

// Create arraies for table "foci_data"
date = newArray;			//An array to store date of experiments
image_file = newArray;  	//An array to store filenames
cell_number =newArray;		//An array to store the serial number of cell
obserbation = newArray;		//An array to store the serial number of observation
foci_serial = newArray;		//An array to store the serial number of fluorescent granule
foci_meanint =newArray;		//An array to store mean fluorescence intensity of the granule
foci_x = newArray;			//An array to store the x coordinate of the granule
foci_y = newArray;			//An array to store the y coordinate of the granule
cell_area = newArray;		//An array to store the area of cell
prominences = newArray;		//An array to store the current prominence
withfoci = newArray;		//An array to check existance of granule in the cell
defocuses = newArray;  // An array to check cells that are not suitable for analysis

print("\\Clear"); 							
title = getTitle();
// Remove the extension from the filename
title_s = replace(title, "\\.nd2", ""); 
title_s = replace(title_s, "\\.tif", "");
title_s = replace(title_s, "\\.tiff", "");

imagewidth = getWidth();   // Acquire size of the image
imageheight= getHeight();

run("Z Project...", "projection=[Max Intensity]");
title2 = getTitle();
run("Split Channels");

c1 = "C1-" + title2; // Channel #1：A z-stack fluorescence image
c2 = "C2-" + title2; // Channel #2：x60 phase-contrast image

selectWindow(c2);
BFID = getImageID();
run("Duplicate...", " ");
run("Gaussian Blur...", "sigma=1");   // sigma = 2 for x100 phase-contrast image 
run("Subtract Background...", "rolling=25");  // rolling = 42 for x100 phase-contrast image 

setAutoThreshold("Li dark");
setOption("BlackBackground", true);
run("Convert to Mask");

run("Dilate");
//run("Dilate");  //repeat if x100 phase-contrast image 

run("Skeletonize");
run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune_0");
run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune_0");  // Repeat twice to ensure generation of a "Tagged skeleton" image
selectWindow("Tagged skeleton");
run("Convert to Mask");
rename("debranch_" + title2);

run("Dilate");
run("Close-");

c3 = "debranch_" + title2;

run("Duplicate...", " ");
run("Fill Holes");
rename("filled_" + title2);

c4 = "filled_" + title2;

imageCalculator("Subtract create", c4,c3);
selectWindow("Result of " + c4);
rename("segmented_" + title2);
run("Set Measurements...", " area redirect=None decimal=3");
run("Analyze Particles...", "size=400-3000 pixel show=Outlines display exclude clear add"); // size = 667-5000 for x100 phase-contrast image 

sno = 0; // A variable to count observation： Serial No. of observation, equals total number of rows in table "foci_data"
cwf = 0; // A variable to count the number of cells showing fluorescent foci: Cell With Foci
snf = 0; // A variable to count the total number of fluorescent foci: Serial No. of Foci

selectWindow(c1);							// Superimpose cell boundaries detected (ROIs) on a fluorescence image
//run("Subtract Background...", "rolling=10");  // Important step to detect fine fluorescent condensates (granules) 
roiManager("Set Color", "white");
roiManager("Show None");
roiManager("Show All");
roiManager("Save", dirDR + "ROIs-" + title2 + ".zip"); 
run("Clear Results");

for(k=0; k<roiManager("count"); k++) {   // Activate and analyse a ROI one by one
	run("Clear Results");
	roiManager("select", k);
	Roi.setStrokeColor(255*random,255*random,255*random);
	area = getValue("Area");		// Store the area of the ROI in a variable		
	run("Find Maxima...", "prominence=" + prom + " exclude output=Count");// prominence = noise tolerance
	
	nf_k = getResult("Count", 0);	//Store the number of maxima = the number of foci in a variable
									// number of foci in k-th cell
	run("Clear Results");			
	if(nf_k >= 1){						// If the k-th cell showed some foci
		cwf = cwf + 1;
		
		run("Find Maxima...", "prominence=" + prom + " exclude output=List");//Detect maxima (fluorescent foci), and show their x-y coordinates in the Result table
		for (l=0; l< nf_k; l++){						// Analyze {nf_k} foci in the k-th cell
			obserbation[sno] = sno+1;
			foci_serial[sno] = snf+1;					
			cell_number[sno] = k+1;						
			image_file[sno] = title_s;					
			
			foci_x[sno] = getResult("X", l);			// The x coordinate of the snf-th granule
			foci_y[sno] = getResult("Y", l);			// The y coordinate of the snf-th granule
			cell_area[sno] = area;
			prominences[sno] = prom;
			date[sno] = edate;
			withfoci[sno]="TRUE";

			run("Specify...", "width=2 height=2 x="+foci_x[sno]+" y="+foci_y[sno]+" slice=1 oval centered");
			foci_meanint[sno] = getValue("Mean");		// Measure mean intensity of the ROI
	
			Table.update;
			sno = sno+1;
			snf= snf+1;		
		}
	}else {											// nf_k == 0, i.e., if the k-th cell showed no foci
		obserbation[sno] = sno+1;
		foci_serial[sno] = "NaN"; 			
		cell_number[sno] = k+1;
		image_file[sno] = title_s;		
		foci_x[sno] = "NaN";		
		foci_y[sno] = "NaN";		
		cell_area[sno] = area;
		foci_meanint[sno] = "NaN";
		prominences[sno] = prom;
		date[sno] = edate;
		withfoci[sno]="FALSE";
		Table.update;
		sno=sno+1;
}
}
file_name[i] = title_s;									
total_cell_number[i] = k;								
foci_cell[i] = cwf;									
pct_foci_cell[i] = 100*cwf/k;	// Calculate and store the percentage of cell showing foci in the i-th iamge file
prominence[i] = prom;

saveAs("Tiff", dirGreen+getTitle());

c1t = replace(c1, "\\.nd2", "\\.tif");
					
// Create a merge image
run("Merge Channels...", "c4=["+c2+"] c6=["+c1t+"] keep ignore"); // Sandwich variables with {"+} when using in option
selectWindow("RGB");
rename("BFM-" + title);
RGID = getImageID();
saveAs("Tiff", dirMer+getTitle());
close();

selectImage(BFID);
saveAs("Tiff", dirBF+getTitle());
selectWindow("Drawing of segmented_"+ title2);
saveAs("Tiff", dirDR+getTitle());
close();

// Create a mask image of foci detected
newImage(title_s+"_foci_location", "8-bit", imagewidth, imageheight, 1);
setColor(0,0,0);
fill();
//setLineWidth(1);
for (m = 0; m < foci_x.length; m++) {
	setColor(255,255,255);
	fillOval(foci_x[m]-2, foci_y[m]-2, 4, 4);
}
saveAs("Tiff", dirDR+getTitle());
close();


// Create a table summarizing data of individual granules
Array.show("foci_data", date,image_file,cell_number, obserbation,withfoci, foci_serial, foci_meanint, foci_x, foci_y, cell_area, prominences, defocuses);
selectWindow("foci_data");
saveAs("Results", dirCSV + edate1 + "_"+ title_s + "_prom" + prom + ".csv"); 
run("Close");
run("Close All");

}
    
// Create a table summarizing data of individual image files
Array.show("foci_stat(row numbers)", file_name, total_cell_number, foci_cell, pct_foci_cell, prominence, defocus);														
    selectWindow("foci_stat");
    saveAs("Results", dirD +"/"+ edate1 + "_prom" + prom + "_foci_stat.csv"); 
   
	run("Close");	
    run("Clear Results");						// Reset Results window
	print("\\Clear"); 							// Reset Log window
	roiManager("reset");						// Reset ROI manager

  showMessage(" ", "<html>"+"<font size=+2>Process completed<br>");
}
