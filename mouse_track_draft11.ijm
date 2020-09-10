//Value for mean filtering speed, angle, and dTheta values
meanFilter = 4;
tailLength = 0;
houghThresh = 0.6;
knownRadius = 15; //Known radius of retroreflective tracking dots
useKnownRadius = true; //Whether to use known or measured radii to set area for finding orientation dot

run("Bio-Formats Macro Extensions");
dir = getDirectory("Choose a Directory ");
setBatchMode(true);
count = 0;
countFiles(dir);
n = 0;
processFiles(dir);
setBatchMode(true);

function countFiles(dir) {
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/")) countFiles(""+dir+list[i]);

		//Check to make sure a file is an 8-bit tiff with more than 100 slices
		else if (endsWith(list[i], ".tif")){
			Ext.setId(dir+list[i]);
			Ext.getSizeZ(sizeZ);
			Ext.getSizeT(sizeT);
			Ext.getSizeC(sizeC);
			Ext.getPixelType(pixelType);
			if(sizeZ*sizeT > 100  && matches(pixelType, "uint8") && sizeC == 1) count++;
		}
	}
}

function processFiles(dir) {
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/")) processFiles(""+dir+list[i]);
		else if (endsWith(list[i], ".tif")){
			Ext.setId(dir+list[i]);
			Ext.getSizeZ(sizeZ);
			Ext.getSizeT(sizeT);
			Ext.getSizeC(sizeC);
			Ext.getPixelType(pixelType);

			//Sometimes time domain is reported as Z domain, so this ensures that the macro works either way
			if(sizeZ*sizeT > 100 && matches(pixelType, "uint8") && sizeC == 1){
				showProgress(n++, count);
				print("Processing file " + n + " of " + count + ".");
				path = dir+list[i];
				processFile(path);
			}
		}
	}
}

function processFile(file) {
	//Close images and results
	close("*");
	if (isOpen("Results")) { 
		selectWindow("Results"); 
		run("Close"); 
	}
	
	open(file);
	directory = File.directory();
	title = getTitle();
	getDimensions(width, height, channels, slices, frames);
	run("Set Measurements...", "center redirect=None decimal=3");
	
	//Create binary mask of edges in movie
	run("Duplicate...", "title=mask duplicate");
	selectWindow("mask");
	run("Find Edges", "stack");
	nBefore = nImages; //Variable to keep track of when Hough transform finishes
	setAutoThreshold("Triangle dark");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Moments background=Dark calculate black");
	
	//NOTE: 0.1 allows for local tracking of background when there is no circle, and it fails to recover (i.e. not stingent enough).
	run("Hough Circle Transform","minRadius=10, maxRadius=24, inc=1, minCircles=1, maxCircles=1, threshold=" + houghThresh + ", resolution=100, ratio=1.0, bandwidth=5, local_radius=20,  reduce local_search show_centroids show_mask results_table");
	
	//wait until hough transform finishes
	while(nImages == nBefore){
		wait(100);
	}

	close("mask");
	selectWindow(title);
	run("Duplicate...", "title=Dot duplicate");
	run("Median...", "radius=2 stack");
	run("Invert", "stack");
	
	//Isolate the position of the dot
	selectWindow("Dot");

	for(a=0; a<nResults; a++){
		setSlice(getResult("Frame (slice #)",a));
		x = getResult("X (pixels)", a);		
		y = getResult("Y (pixels)", a);
		if(useKnownRadius) r = knownRadius;
		else r = getResult("Radius (pixels)", a);
		makeOval(x-r, y-r, 2*r, 2*r);
		getStatistics(dummy, dummy, min, max);
		backtrack = false; //Record if prominence search is backtracking - gaurds against infinite loop if two points have identical prominence
		for(b=1; b<256; b=b){ //Look for maximum point in selection
			makeOval(x-r, y-r, 2*r, 2*r);
			run("Find Maxima...", "noise="+b+" output=[Point Selection]");
			getSelectionCoordinates(xpoints, ypoints);
			if(xpoints.length > 1 && !backtrack) b *= 2;  //Search by power to speed search for max point
			else if (xpoints.length < 1){
				b--;
				backtrack = true;
			}
			else b = 257;	
			if(b==256) b = 255; //If power of 2 reaches 256, check back to 255.
		}
		if(xpoints.length > 1){ //If there are still multiple points, find the one closest to the prior point
			minIndex = -1;
			minDist = 99999999 * (knownRadius+1);
			xT = getResult("xT", a-1); //Get previous tip position
			yT = getResult("yT", a-1);
			for(b=0; b<xpoints.length; b++){
				dist = (xT-xpoints[b])*(xT-xpoints[b]) + (yT-ypoints[b])*(yT-ypoints[b]);
				if(dist < minDist){
					minDist = dist;
					x = xpoints[b];
					y = ypoints[b];
				}
			}
			print(a);
		}
		else getSelectionBounds(x, y, dummy, dummy);
		setResult("xT", a, x);
		setResult("yT", a, y);
		run("Select None");
	}
	updateResults();
	close("Dot");
	selectWindow(title);
	
	//Calculate theta and velocity
	dThetaArray = newArray(nResults);
	angleArray = newArray(nResults);
	speedArray = newArray(nResults);
	xCarray = newArray(nResults);
	yCarray = newArray(nResults);
	xTarray = newArray(nResults);
	yTarray = newArray(nResults);
	
	setResult("Speed",0,0);
	xAngle = -1;
	yAngle = 0.1;
	hyp = sqrt(xAngle*xAngle + yAngle*yAngle);
	if(xAngle<0 && yAngle>=0) radAngle = PI - asin(yAngle/hyp);
	else if(xAngle<0 && yAngle<0) radAngle = (-1)*PI - asin(yAngle/hyp);
	else radAngle = asin(yAngle/hyp);
	degAngle = radAngle*180/PI;
	setResult("Theta",0,degAngle);
	setResult("dTheta",0,0);
	dThetaArray[0] = 0;
	angleArray[0] = degAngle;
	speedArray[0] = 0;
	xCarray[0] = getResult("X (pixels)",0);
	yCarray[0] = getResult("Y (pixels)",0);
	xTarray[0] = getResult("xT",0);
	yTarray[0] = getResult("yT",0);
	
	
	for(a=1; a<nResults; a++){
		//Measure speed
		dTime = getResult("Frame (slice #)", a) - getResult("Frame (slice #)", a-1);
		dX = getResult("X (pixels)", a) - getResult("X (pixels)", a-1);
		dY = getResult("Y (pixels)", a) - getResult("Y (pixels)", a-1);
		speed = sqrt(dX*dX + dY*dY)/dTime;
		setResult("Speed",a,speed);
	
		//Measure angle
		xAngle = getResult("xT", a) - getResult("X (pixels)", a);
		yAngle = getResult("yT", a) - getResult("Y (pixels)", a);
		hyp = sqrt(xAngle*xAngle + yAngle*yAngle);
		if(xAngle<0 && yAngle>=0) radAngle = PI - asin(yAngle/hyp);
		else if(xAngle<0 && yAngle<0) radAngle = (-1)*PI - asin(yAngle/hyp);
		else radAngle = asin(yAngle/hyp);
		degAngle = radAngle*180/PI;
		setResult("Theta",a,degAngle);
	
		//Measure change in angle
		dTheta = ((sin(degAngle) - sin(getResult("Theta",a-1)))*180/PI)/dTime;
		setResult("dTheta",a, asin(dTheta));
	
		//Update arrays
		dThetaArray[a] = dTheta;
		angleArray[a] = degAngle;
		speedArray[a] = speed;
		xCarray[a] = getResult("X (pixels)",a);
		yCarray[a] = getResult("Y (pixels)",a);
		xTarray[a] = getResult("xT",a);
		yTarray[a] = getResult("yT",a);
	}
	updateResults();
	
	//Run mean filter on arrays to smooth out pixel jitter while preserving high-freq components
	newImage("dTheta", "32-bit black", dThetaArray.length, 1, 1);
	newImage("angle", "32-bit black", dThetaArray.length, 1, 1);
	newImage("speed", "32-bit black", dThetaArray.length, 1, 1);
	newImage("xC", "32-bit black", dThetaArray.length, 1, 1);
	newImage("yC", "32-bit black", dThetaArray.length, 1, 1);
	newImage("xT", "32-bit black", dThetaArray.length, 1, 1);
	newImage("yT", "32-bit black", dThetaArray.length, 1, 1);
	for(a=0; a<7; a++){
		if(a==0) selectWindow("speed");
		else if (a==1) selectWindow("angle");
		else if (a==2) selectWindow("dTheta");
		else if (a==3) selectWindow("xC");
		else if (a==4) selectWindow("yC");
		else if (a==5) selectWindow("xT");
		else if (a==6) selectWindow("yT");
				
		for(b=0; b<dThetaArray.length; b++){
			if(a==0) setPixel(b,0,speedArray[b]);
			else if (a==1) setPixel(b,0,angleArray[b]);
			else if (a==2) setPixel(b,0,dThetaArray[b]);
			else if (a==3) setPixel(b,0,xCarray[b]);
			else if (a==4) setPixel(b,0,yCarray[b]);
			else if (a==5) setPixel(b,0,xTarray[b]);
			else if (a==6) setPixel(b,0,yTarray[b]);
		}
	
		run("Mean...", "radius=" + meanFilter);
	
		for(b=0; b<dThetaArray.length; b++){
			if(a==0) setResult("speed (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==1) setResult("Theta (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==2) setResult("dTheta (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==3) setResult("xC (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==4) setResult("yC (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==5) setResult("xT (mean " + meanFilter + ")",b,getPixel(b,0));
			else if (a==6) setResult("yT (mean " + meanFilter + ")",b,getPixel(b,0));
		}
	}
	updateResults();
	
	//calculate 8-bit LUT for speed 
	selectWindow("speed");
	getStatistics(dummy, dummy, min, max);
	setMinAndMax(min,max);
	run("8-bit");
	run("Add...", "value=1");
	dLUT = newArray(dThetaArray.length);
	for(b=0; b<dLUT.length; b++){
		dLUT[b] = getPixel(b,0);
	}
	close("speed");
	close("dTheta");
	close("angle");
	close("xC");
	close("yC");
	close("xT");
	close("yT");
	
	
	//Generate angle masks
	//Draw track arrows with dragon tails
	newImage("Arrow track", "8-bit black", width, height, slices);
	for(a=0; a<nResults; a++){
		setSlice(getResult("Frame (slice #)",a));
		for(b=a-tailLength; b<=a; b++){
			if(b>=0){
				makeArrow(getResult("X (pixels)", a-0), getResult("Y (pixels)", a-0), getResult("xT", a-0), getResult("yT", a-0), "filled");
				Roi.setStrokeWidth(2);
				Roi.setStrokeColor("white");
				setForegroundColor(dLUT[b], dLUT[b], dLUT[b]);
				run("Draw", "slice");
				run("Select None");
			}
		}
	}
	
	//Apply physics LUT with black background
	run("physics");
	getLut(reds, greens, blues);
	reds[0] = 0;
	greens[0] = 0;
	blues[0] = 0;
	setLut(reds, greens, blues);
	selectWindow("Arrow track");
	run("RGB Color");
	selectWindow(title);
	run("RGB Color");
	imageCalculator("Transparent-zero create stack", title,"Arrow track");
	selectWindow("Result of " + title);
	saveAs("tiff", directory + "Arrow track for " + title);
	close(title);
	close("Arrow track");
	saveAs("Results", directory + "Results for " + replace(title, ".tif$", ".csv"));
}

