// Draw and segment image ROIs for downstream FLIM analysis

run("Close All"); roiManager("Reset"); run("Collect Garbage");

// Select and open reference image for deliniating ROIs; set up seperate directory for ROI and analysis files
file_path = File.openDialog("Select reference file");
dir_path = File.getParent(file_path) + File.separator; 
dir_parts = split(dir_path,File.separator);
dir_name = Array.slice(dir_parts,dir_parts.length-1); // Array.show(dir_name);
dir_name = String.join(dir_name); // print(dir_name); 
anlys_parentDir = dir_path + dir_name + "_analysis"; // print(anlys_parentDir);
File.makeDirectory(anlys_parentDir);  

// open(file_path);
run("Bio-Formats Importer","open=[" + file_path + "]");
filename = getTitle();
rootname = File.nameWithoutExtension;
anlys_dir = anlys_parentDir + File.separator + rootname; // print(anlys_dir);
File.makeDirectory(anlys_dir);
run("Z Project...", "projection=[Max Intensity]");

//==================Selection of ROIs for analysis================================
select_roi = 1; // Boolean for continuing to select ROIs
roi_num = 0; // selected ROI number
setTool("rectangle");

while (select_roi == 1){

    waitForUser("Draw the ROI for analysis and click OK");

    // record ROI in manager and rename
    roiManager("Add");
    roi_ID = roiManager("count") - 1; // ROI manager starts at zero
    roiManager("Select", roi_ID);
    roi_num = roi_num + 1;
    roiManager("Rename", "roi_" + roi_num);

    // create cropped subroi composite image and save
    selectImage(filename);
    roiManager("Select", roi_ID);
    roi_name = Roi.getName;
    run("Duplicate...", "duplicate");
    rename(rootname + "_" + roi_name + ".tif"); crop_name = getTitle();
    crop_path = anlys_dir + File.separator + crop_name;
    save(crop_path); close(crop_name);

    roiManager("Deselect");		    
    select_roi = getBoolean("Do you want to select more ROIs?");

}

close(filename);
close("Max_" + filename);

roiManager("Deselect");
roi_path = anlys_dir + File.separator + rootname + "_ROIset.zip"; 
roiManager("Save","" + roi_path);
roiManager("Reset");