/*
 * SOP SCRIPT: AUTOMATED ALGAE COUNTER (Nikon .nd2 Version)
 */

macro "SOP Algae ND2 [F2]" {
    // --- CONFIGURATION ---
    minSize = 50;          
    boxSize = 512; 
    threshMethod = "Triangle"; 
    
    // --- SETUP ---
    run("Close All");
    if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
    print("\\Clear");
    
    inputDir = getDirectory("Select the INPUT Folder");
    outputDir = getDirectory("Select the OUTPUT Folder");
    
    list = getFileList(inputDir);
    Array.sort(list); 

    // Create CSV and write header
    csvPath = outputDir + "Algae_Results_ND2.csv";
    File.saveString("Sample No,Conc (Cells/mL),Average,Counts...\n", csvPath);

    currentSample = "";
    counts = newArray(0); 
    imageCount = 0;

    // Batch mode can be tricky with Bio-Formats; we'll use it but keep an eye on it
    setBatchMode(true); 

    for (i = 0; i < list.length; i++) {
        filename = list[i];
        
        // --- FILTER FOR ND2 FILES ---
        if (endsWith(toLowerCase(filename), ".nd2")) {
            imageCount++;
            showStatus("Opening ND2: " + filename);

            // 1. Parse Sample ID
            nameNoExt = replace(filename, ".nd2", "");
            if (indexOf(nameNoExt, "_") >= 0) {
                parts = split(nameNoExt, "_");
                sampleID = parts[0]; 
            } else {
                sampleID = nameNoExt;
            }

            // 2. Data Grouping
            if (sampleID != currentSample) {
                if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
                currentSample = sampleID;
                counts = newArray(0); 
            }
            
            // 3. Image Prep using Bio-Formats
            // "color_mode=Default" ensures it opens as intended by the microscope
            path = inputDir + filename;
            run("Bio-Formats Importer", "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
            originalTitle = getTitle();
            
            // Auto-Crop
            makeRectangle(getWidth()/2 - boxSize/2, getHeight()/2 - boxSize/2, boxSize, boxSize);
            run("Crop");

            // 4. Analysis
            run("Duplicate...", "title=Detection_Mask");
            run("8-bit");
            run("Subtract Background...", "rolling=50");
            setAutoThreshold(threshMethod + " dark");
            run("Convert to Mask");
            run("Watershed");
            
            // Count Particles
            roiManager("Reset");
            run("Analyze Particles...", "size=" + minSize + "-Infinity circularity=0.30-1.00 exclude add");
            
            thisCount = roiManager("count");
            counts = Array.concat(counts, thisCount);
            
            // 5. Save Evidence (Flattened JPEG)
            selectWindow(originalTitle);
            run("8-bit");
            roiManager("Show All with labels");
            run("Labels...", "color=Cyan font=14 show use draw");
            run("Flatten"); 
            saveAs("Jpeg", outputDir + "Checked_" + nameNoExt + ".jpg");
            
            // 6. Cleanup
            run("Close All"); 
        }
    }

    // Save final sample
    if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
    
    setBatchMode(false);
    
    if (imageCount == 0) {
        showMessage("Error", "No .nd2 files found. Ensure the files end in .nd2 (lowercase or uppercase).");
    } else {
        showMessage("Success!", imageCount + " .nd2 images processed.\nResults saved to: " + csvPath);
    }
}

function saveRowToCSV(sample, cArr, path) {
    sum = 0; n = cArr.length;
    for (k=0; k<n; k++) { sum = sum + cArr[k]; }
    avg = 0; if (n > 0) { avg = sum / n; }
    conc = avg * 10000; 
    line = sample + "," + conc + "," + avg;
    for (k=0; k<n; k++) { line = line + "," + cArr[k]; }
    File.append(line + "\n", path);
}