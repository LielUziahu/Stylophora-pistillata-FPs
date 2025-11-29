/*
 * SOP SCRIPT: AUTOMATED ALGAE COUNTER (Variable Replicates)
 * Lab Protocol: 11/2025
 * * PARAMETERS:
 * - Crop Size: 512x512 (Calculated for 10X Binning 3)
 * - Min Size: 50 pixels (Filters dust <14px)
 * - Threshold: Triangle (Optimized for dim fluorescence)
 */

macro "SOP Algae Variable [F2]" {
    // --- CONFIGURATION ---
    minSize = 50;          
    boxSize = 512; 
    threshMethod = "Triangle"; 
    
    // --- SETUP ---
    run("Close All");
    print("\\Clear");
    roiManager("Reset");
    
    // Folder Selection
    inputDir = getDirectory("Select the INPUT Folder (Images)");
    outputDir = getDirectory("Select the OUTPUT Folder (Results)");
    
    list = getFileList(inputDir);
    Array.sort(list); 
    
    // Excel Header
    print("Sample No\tConc (Cells/mL)\tAverage\tCount 1\tCount 2\tCount 3\t(etc...)");

    currentSample = "";
    counts = newArray(0); 
    
    // Interactive Mode (BatchMode OFF to allow erasing)
    setBatchMode(false); 

    for (i = 0; i < list.length; i++) {
        filename = list[i];
        if (endsWith(filename, ".tif") || endsWith(filename, ".jpg") || endsWith(filename, ".png")) {
            
            // 1. Parse Sample ID
            dotIndex = lastIndexOf(filename, ".");
            nameNoExt = substring(filename, 0, dotIndex);
            parts = split(nameNoExt, "_");
            sampleID = parts[0]; 

            // 2. Data Grouping
            if (sampleID != currentSample) {
                if (currentSample != "") { printRow(currentSample, counts); }
                currentSample = sampleID;
                counts = newArray(0); 
            }
            
            // 3. Image Prep
            open(inputDir + filename);
            originalTitle = getTitle();
            
            // Auto-Crop to 0.1 uL Volume
            makeRectangle(getWidth()/2 - boxSize/2, getHeight()/2 - boxSize/2, boxSize, boxSize);
            run("Crop");
            
            // 4. USER INTERACTION: Scale Bar Eraser
            setTool("rectangle");
            waitForUser("SOP Action", "Draw a box over the SCALE BAR text.\nThen click OK.");
            
            if (selectionType() != -1) {
                run("Set...", "value=0"); // Paint black
                run("Select None");
            }

            // 5. Analysis
            run("Duplicate...", "title=Detection_Mask");
            run("8-bit");
            run("Subtract Background...", "rolling=50");
            setAutoThreshold(threshMethod + " dark");
            run("Convert to Mask");
            run("Watershed");
            
            // Count
            roiManager("Reset");
            run("Analyze Particles...", "size=" + minSize + "-Infinity circularity=0.30-1.00 exclude add");
            
            thisCount = roiManager("count");
            counts = Array.concat(counts, thisCount);
            
            // 6. Save Evidence
            selectWindow(originalTitle);
            run("8-bit");
            roiManager("Show All with labels");
            run("Labels...", "color=Cyan font=14 show use draw");
            run("Flatten"); 
            saveAs("Jpeg", outputDir + "Checked_" + filename);
            
            close("*"); 
            roiManager("Reset");
        }
    }

    // Print Final Row
    if (currentSample != "") { printRow(currentSample, counts); }
    
    showMessage("SOP Complete", "All images processed.\n\n1. Click inside the LOG window.\n2. Press Ctrl+A (Select All).\n3. Press Ctrl+C (Copy).\n4. Paste into Excel.");
}

function printRow(sample, cArr) {
    sum = 0; n = cArr.length;
    
    // Calculate Stats
    for (k=0; k<n; k++) { sum = sum + cArr[k]; }
    avg = 0; if (n > 0) { avg = sum / n; }
    conc = avg * 10000; 
    
    // Stats First Layout
    line = sample + "\t" + conc + "\t" + avg;
    
    // Append Variable Counts
    for (k=0; k<n; k++) {
        line = line + "\t" + cArr[k];
    }
    print(line);
}
