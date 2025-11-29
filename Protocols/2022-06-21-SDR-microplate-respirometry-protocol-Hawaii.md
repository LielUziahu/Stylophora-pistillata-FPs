---
layout: post
title: SDR microplate respirometry protocol for coral planulae
category: [Respiration, Photosynthesis, Protocol ]
tags: [ Coral, Planulae, Physiology, pH, Temperature,Montipora capitata ]
---

This protocol details the step-by-step use of the PreSens/Loligo Microplate respiration system for measuring respiration and photosynthesis of coral planulae.  

[Original protocol](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Early-Life-History-Respirometry-Protocol/) written by Ariana Huffmyer and modified by [Jill Ashey](https://github.com/JillAshey) and me during our time at the Hawaiʻi Institute of Marine Biology in June-July 2022, conducting temperature and pH stress experiments with planulae and spats of the corals _Pocillopora acuta_ and _Montipora capitata_.


### Prepare equipment  
-	[PreSens Oxy-10 Sensor Dish Reader (SDR)](https://www.presens.de/products/detail/sdr-sensordish-reader-basic-set) in combination with [Loligo 24-well microplate](https://loligosystems.com/24-well-glass-microplate-80-ul) (80 µL each well)
-	Glass coverslips
-	PreSens software downloaded with calibrations entered as provided by manufacturer
-	Aquarium lighting [Prime 16HD AI](https://www.aquaillumination.com/products/prime)
-	Digital Incubator [MyTemp mini](https://www.benchmarkscientific.com/rpproducts/mytemp-mini-digital-incubator/)
-	Filtered seawater (1 µm)
-	MilliQ filtered water
-	Bubbler and aeration
-	Glass Transfer Pipets
-	PC computer
-	Kimwipes

### Prepare planulae
-	Take planulae out of the experimental system and rinse them in filtered seawater to clean any debris off of them 
-	Based on previous respirometry trials by Ariana Huffmyer, 5 _P. acuta_ larvae per well or 6 _M. capitata_ larvae per well will give good signals
-	Use a disposable pipette to transfer planulae into a 6 or 24 well plate and separate them out to easily get the desired number of planulae 

### Calibrate oxygen spots 
-	Start by hydrating sensors in the wells ~20 minutes before using the plate with filtered seawater
-	After hydration, fill the wells with 100% air saturated filtered seawater, which was aerated for about 20-30 minutes with the bubbler 
-	Place the 24 well plate on the SDR inside the incubator ![incubator_resp](https://github.com/fscucchia/FScucchia_Lab_Notebook-Mass_Lab/blob/master/images/incubator_resp.jpg?raw=true "incubator_resp") ![incubator_plate](https://github.com/fscucchia/FScucchia_Lab_Notebook-Mass_Lab/blob/master/images/incubator_plate.jpg?raw=true "incubator_plate")
-	Open the PreSens software
-	Make sure the SDR is connected to the computer
-	Select “Connect SDR” and you should see SDR 710 and SDR 641 show up
-	Select one of the SDR tabs
-	Select “Oxygen” from the drop down menu. Select PSt5-1624-01_25C under Batch No. 
-	For each SDR (either 641 or 710), enter the correct Meas. Temp. depending on what incubator it is in and what temperature the incubator is currently at 
-	Enter 37 for the salinity for both SDRs 
-	Select pO2 (%air sat) under Oxygen Unit  
-	After these data are entered, select one SDR tab
-	Click the calibration menu and select One-Point Adjustment
-	Enter 100% at the bottom and click “OK”  
-	Perform this calibration at the start of each measuring day  

### Prepare plates and metadata  
-	Before preparing planulae, turn on the incubators and lights to the desired temperature and level to allow the incubator to come to temperature. Set the lights at 30% intensity on all channels with 0% on the UV channel. Measure the light at the level of the SDR plate and record 
-	Place a cover over the doors of the incubator to block ambient light 
-	*Note*: Only use kimwipes and be very gentle and careful with the plates. Be careful not to touch or scratch the oxygen sensor spots in the bottom of each well 
-	Before starting trials, set up a plate map and a metadata file   
-	Include 2-3 wells to serve as blanks, loaded with filtered seawater. Record the location of blanks on the plate map and record NA on the plate map for any wells that are missing sensors. Change the location of blanks in each run 
-	Well sensors need to be hydrated 20 minutes before loading the planulae as specified above  
-	Put the planulae in the wells and add a drop of filtered seawater to fill the well so that a convex meniscus forms and the water tension holds a tiny bulge slightly above the top of the plate. Ensure the plate is very dry to allow for this tension to form. If the plate is too wet a bubble will not form. If you need to refill or reset a well make sure to dry it with kimwipes. 
![plate1](https://github.com/fscucchia/FScucchia_Lab_Notebook-Mass_Lab/blob/master/images/plate1.jpg?raw=true "plate1")
-	Check the well for bubbles and pipette up and down gently to remove if needed. 
-	Carefully drop a micro coverslip over each well to seal it. It’s best to drop the coverslip gently, and at an angle to the well plate. 
![plate2](https://github.com/fscucchia/FScucchia_Lab_Notebook-Mass_Lab/blob/master/images/plate2.jpg?raw=true "plate2")
-	Check to make sure there are no bubbles in the wells. 
-	Put the plate on the SDR in the incubator. Make sure the oxygen sensors are right on top of the readers on the SDR. Ensure that A1 is aligned between the plate and the SDR.
-	Check the temperature in the incubator with the digital thermometer regularly and launch hobo loggers to log temperature during the run.  

### Photosynthesis and respiration measurements 
-	For both SDR tabs, press “Log measurements”. The software will ask to save the file before starting the run. Save the file in the Hawaii 2022 folder on the PC. Name the file with date, run number, temperature, and initials. 
-	Make sure the lights are on and the temperature is correct. 
-	Once a measurement has been logged for both SDR’s, you can click START.
-	Start the run with the light on for photosynthesis measurements (~30-40 minutes). 
-	Then turn off the light for light-enhanced respiration measurements (~20-25 minutes). 
-	Record start and stop times and the time the light was turned off. Record any notes on wells that have bubbles or other issues. 
-	You should see a general increase in photosynthesis rates and a decrease in respiration rates. The respiration rates will be much more stable than the photosynthesis rates due to the effect of light. 

### Close and save data  
-	At the end of the run, press “Export measurements” to save the data. Do not press close measurement before exporting!
-	Save using the same naming as you did to log the measurements. 
-	Upload data sheets to the google drive or save onto a USB.  

### Clean equipment  
-	Use a glass pipette to remove the water and planulae from the wells GENTLY. 
-	Rinse the microplate with milliQ water and dry using a kimwipe. Roll the kimwipe up and gently dry each individual well being careful not to damage the spot. 
-	Because the oxygen spots are sensitive to light, put it in the bubble wrap and put it in the dark bag. Open the bag slightly so the plate can fully dry. 




