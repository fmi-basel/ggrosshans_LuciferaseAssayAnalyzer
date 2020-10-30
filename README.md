Luciferase Assay Analyzer 

# Overview

1. Description
2. Installation
3. Input
4. Execution
5. Output
6. References/citation


# (1) Description

This tool is a standalone application for analyzing C. elegans developmental progression using luciferase assay time course data.  
The Luciferase Assay Analyzer is a Graphical User Interface (GUI) based on MATLAB code.  
It performs automatic detection of hatch and molts in the luminescence trace, but also allows the user to manually annotate hatch and molts.


# (2) Installation

## Basic installation of the application for end users
To install the Luciferase Assay Analyzer desktop application, download the Application from github:
click on **tag** (next to branches), click on **releases**. 
Choose the appropriate MyAppInstaller, depending on whether you have windows or mac.
For windows: unzip the executable and follow the instructions in the installation menu
For MAC: unzip the application and follow the instruction in the installation menu. 
If you encounter problems with the installation, please have a look under **Issues** -> **Issue with mac installation**, to set up your proxy server correctly. 
FMI MAC users have experienced installation failure at 99% of the installation. Please retry several times, network issues may be the cause.

### Dependencies
None.

## Developer installation
To convert the MATLAB code into a standalone application, use the MATLAB Application Compiler App. 
The Application Compiler packages MATLAB programs into applications that can run outside of MATLAB. 
Detailed information on the Application Compiler App can be found on the [mathworks website](https://nl.mathworks.com/help/compiler/applicationcompiler-app.html)
1. Open the Application Compiler App: enter `applicationCompiler` in the MATLAB command prompt
![ApplicationCompiler](/Images_README/ApplicationCompiler.png)
2. Add main file: gui_select_molts.m
3. Change the **title** to LuciferaseAssayAnalyzer_windows or LuciferaseAssayAnalyzer_mac and choose the appropriate **version**. 
Note: for creating the Luciferase Assay Analyzer executable for windows the compiling needs to be done on a windows computer, while for creating the Luciferase Assay Analyzer application for mac, the compiling needs to be done on a mac computer. 
4. Once you have selected the main file, files required for the application to run will automatically be uploaded as well. Those are gui_select_molts.fig and inputsdlg.m
5. Fill in the remaining fields (**author name, email, company, summary, description**), **additional installer options** can remain unchanged.
6. Select **custom splash screen**: LuciferaseAssayAnalyzer_logo.jpg
7. Click on the little icon next to the title to select the icon that will appear on the desktop: LuciferaseAssayAnalyzer_logo.jpg, and select use mask. When you click save and use, a new window will pop up which allows you to save the all the settings in the compiler in a so-called project. This allows you to load the project and do the compiling at a later stage.
8. In the menu bar, choose: **Runtime downloaded from web**
9. To change the output directories of the testing files, end user files and packaged installers,  choose: **Settings**. 
10. To start the compiling process, click on: **Package**.
11. The executable for end users will appear in the folder: for_redistribution
12. These executables should be zipped manually
13. A new release can be drafted under: **tag** (next to branches) -> **releases** -> **draft a new release**

### Dependencies
MATLAB R2016a

### Run tests
The folder *tests* contains two example data sets (.txt files) and the corresponding experimental designs (pd.xlsx). 
These datasets can be used to test the application.

# (3) Input

The Luciferase Assay Analyzer is generated for the data coming from the Berthold Centro XS3 LB960 machine with the MikroWin software version 5.19. 
The input files for the Luciferase Assay Analyzer are: 
a. a text file with the data from the luminometer
b. an excel file with the experimental design of the 384-wells plate
Below you will find a description of both of these files, and examples can be found in the folder *tests*.
It is important that the format of both files matches exactly the description below.

## (a) Luminometer data file
For a data file to be compatible with the Luciferase Assay Analyzer, it is important that the data and the export setting are correct.
- For running a time course experiment on the luminometer, change the settings: 
**settings** -> **measurement** -> **repeated**: 
	- Total time (total duration assay)
	- Counting time (how long each measurement takes)
	- Cycle time in seconds (after how much time the measurement is repeated) 

- Once the experiment is finished, the data can be exported. 
Please use exactly the following format below.
To set up this format: 
**installation** -> **driver** -> **export** -> **RawData Export Driver** (double click)
	- data layout: by well
	- kinetic layout: Time/Position
	- time format: ss
	- operation mode: no box ticked
	- Statistics export: general statistics
	- format: text file
	- directory: C:\MikroWinData\Transfer

- To export your Luciferase Assay data: 
**file** -> **export** -> **Active Export Driver** (RawData Export Driver)
Choose format: .txt

## (b) Experimental plate design file
For each experiment a plate-design (_pd.xlsx) file should be created. This file tells the Luciferase Assay Analyzer which sample is in each well. 
This file should have the same filename as the .txt file, supplemented with _pd.xlsx. 
For example, for the file 20201231_exp1.txt, the plate design name should be:  20201231_exp1_pd.xlsx. 
See the folder *tests* for examples.
The format of the _pd.xlsx file is the following:
- A1: empty cell 
- A2, A3, A4, ... : A, B, C, ...(alphabet of nr of rows of plate)
- B1, C1, D1, ... : 1, 2, 3, ...(nr of columns of plate)
Fill from B2 with the name of the condition in each well.
Please only use the following characters: ‘a-z’ ‘0-9 ‘-’
Other characters (such as : or / or () or []) are NOT supported and will prevent the Luciferase Assay Analyzer from running
Wells without any sample can be left empty in the plate design file


# (4) Execution
To analyze the developmental progression of C. elegans based on luciferase assay data, please follow the manual below. 
There are critical steps and optional steps. 
The critical steps are required to make proper boxplots and heatmaps. 
If one of the following criteria is not met, the correct output cannot be guaranteed: annotate the hatch of ALL animals, annotate FOUR molts for ALL animals, remove empty wells, CHECK ANNOTATION. These steps are explained in more detail below.

Analysis manual:
1. Open the Luciferase Assay Analyzer application
2. To load your data, click: **Load data**. When you load your dataset for the first time, load the .txt file. When you load your dataset subsequent times, either the .txt file or the .mat file can be chosen
3. Once your data is loaded you will see an individual trace of a single animal. 
In this trace the pre-hatch is already annotated (red dots). To annotate the molts, click: **Detect molts**. Green dots will appear.
4. To manually change the hatch, click: **Annotate Hatch**. A horizontal and vertical line will appear. The vertical line (at any horizontal position) can be used to indicate the hatch. I do not recommend to use the middle of the cross, because the lines disappear at that position.
![LuciferaseAssayAnalyzer](/Images_README/AnnotateHatch.png)
5. If you pipetted L1 animals instead of eggs into the plate, the pre-hatch will not be present. You can check the box **Set all hatch to 0** to remove the red dots and set the hatch of all the animals to zero. If you wish to undo this, unclick **Set all hatch to 0**. The red dots will now appear again.
6. To manually change the molts, use one of the following:
- Case 1: four regions are already detected as molts, but the molt entry or molt exit of these molts are not annotated properly. Unclick the corresponding **Mn valid** (with n=1,2,3,4). The green dots of that molt will disappear and the cross for annotation will appear. First click on the molt entry, next click on the molt exit. New green dots will appear.
- Case 2: four regions are already detected as molts, but one of the 4 real molts is missing and an additional region is detected. Regardless of where this wrong region is located, ALWAYS unclick **M4 valid**, and use the cross to annotate the molt entry and the molt exit of the molt that the algorithm missed. The software will internally rearrange the molt, such that the first detected region is molt 1, the second detected region is molt 2, and so on.
- Case 3: less than 4 molts are detected. Regardless of which molt is missing, ALWAYS unclick **M4 valid**, and use the cross to annotate the molt entry and the molt exit of the molt that the algorithm missed. The software will internally rearrange the molt, such that the first detected region is molt 1, the second detected region is molt 2, and so on.
- Case 4: my strain is a mutant and makes less than, or more than 4 molts. The Luciferase Assay Analyzer is not able to process data with less than 4 molts yet. Therefore, 4 molts ALWAYS have to be annotated. For example, if your mutant strain makes 3 molts, annotate a fake molt at the end of the time course and manually remove the fake 4th molt from the boxplots. The white-black heatmap can be used without any problems, but in this case it is better to not use the blue-yellow molt annotation plot. Another example, if your mutant makes 5 molts, leave the first molt unannotated. Note that the boxplot will now give you the results of molt 2, 3, 4, and 5 (despite them being names 1,2,3, and 4 in the boxplot). Also in this case the white-black heatmap can be used without any problems, but it is better to not use the blue-yellow molt annotation plot.
- Case 5: the well is empty, i.e. the luminescence levels stay low and no molts can be observed, unclick **Valid sample**. Unvalid samples are excluded from heatmaps and boxplots. I would not recommend to set luminescence traces that look weird to you as unvalid. This is the variability in your data and you don’t know whether this reflects true biology or not. 
7. The box **y max** allows to manually change the maximum value of the y-axis. This can be done by typing the number in this box.
8. To switch to the next sample, use **Next sample**. To switch to the previous sample, use **Previous sample**.
9. Alternatively, the box **Go to sample** allows to jump to a specific sample number. This can be done by typing the sample number in this box.
10. Once the hatch is annotated for ALL animals, FOUR molts are annotated for ALL animals and empty wells are removed, the annotation can be checked using: **Check annotation**. Note that the analysis does not work and a correct output cannot be guaranteed if 1 or more molts are not annotated, and if the molts are not in the correct order (which is not visible for the user). Check annotation goes through the annotations of all the animals and check exactly that. When the annotation is not correct for one animal, the application will jump to that animals and display its luminescence trace. Start the annotation from scratch for this trace, by unclicking **Valid sample** (annotations will disappear), re-clicking **Valid sample** (pre-hatch annotation will appear), and clicking **Detect molts** (molt annotation will appear). Now, adjust the hatch and molt annotation appropriately following the instructions above. 
Click again on **Check annotation** until it does not jump to a specific animal any more. Now, you are ready for generating the output files.
11. To export your annotations, use the button **Annotation to Excel**. A csv file will be exported to the same folder as where you data is located with the same name as your .txt file extended with _molt_detection_results.csv. This is tab separated data which can be converted to columns by choosing: **Data** -> **Text to Columns** -> **Delimited** -> **Next** -> **Tab** -> **Next**. This file contains the columns: condition and the time points of the hatch, M1_entry, M1_exit, M2_entry, M2_exit, M3_entry, M3_exit, M4_entry, M4_exit in hours. Unvalid animals are excluded.
12. **Make boxplots** (for stage durations and luminescence values during the molt): a pop-up menu will appear in which you can select the conditions you would like to plot together in one boxplot. A new pop-up menu will appear in which you can select the order of the conditions. Select a unique number for each condition. The number refers to the position in the box plot. The condition with number 1 will be plotted as the first condition, the condition with number 2 will be plotted as the second condition, and so on. Three boxplot will be generated: molts, intermolts, larval stages. Boxplots can be saved: **File** -> **save**. The luminescence during the molt is the log10 of the median luminescence during the molt. The median prevents high values due to measurement errors or mis-annotation to skew the value of the relatively short molt which would otherwise occur with the mean. The log10 is to bring the different conditions closer together for visualization purposes.
![DropDownMenu](/Images_README/dropdownmenu.png)
![Boxplot](/Images_README/BoxplotLS.png)
13. **Make heatmaps** (for stage durations): a pop-up menu will appear in which you can indicate whether you want the pre-hatch to be included in the heatmap, choose yes or no. A new pop-up menu will appear in which you can indicate whether you want to plot all strains together in one heatmap or plot strains in separate heatmaps. Plotting strains separately is generally recommended. A new pop-up menu will appear in which you can indicate on which molt the data should be sorted. For example, if you choose molt 1, the animal that enters molt 1 the earliest will appear on top of the heatmap and the animal that enters molt 1 the latest will appear on the bottom of the heatmap. Two different plots are generated. The blue-yellow plot shows the annotation of the molts in yellow and intermolts in blue. The white-grey-black heatmap shows the trend corrected luminescence values, with values ranging from low (black) to high (white) corresponding to the molts and intermolts respectively. 
Plots can be saved: **File** -> **save**
![Heatmaps](/Images_README/Heatmaps.png) 
14. **Statistics to Excel** (for stage durations and luminescence values during the molt): In the drop-down menu Str1. and Str2. select the strains that you would like to test against one another. When clicking on Statistics to Excel, a csv file will be exported to the same folder as where you data is located with the same name as your .txt file extended with _log10medLumiMolt_stats (for luminescence during the molt) or _DevStages_stats_(for stage durations), including the strains tested against one another. This is tab separated data which can be converted to columns by choosing: **Data** -> **Text to Columns** -> **Delimited** -> **Next** -> **Tab** -> **Next**. These files contain t-test statistics and variance test statistics. 
NOTE: The results presented in these csv files can be used to get a first indication about the statistical differences. However, it is very important that you think for yourself which test is appropriate for your data (parametric vs non-parametric test). Different statistical tests have different assumptions. Therefore, please do NOT use the default tests provided here without thinking about it. 


# (5) Output

The output files of this application include:
- Csv file with molt annotation in hours
- Csv file with statistics of stage durations (to get a first impression of the significance)
- Csv file with statistics of the log10 of the median luminescence during the molt (to get a first impression of the significance)
- Boxplot of stage durations and the log 10 of the median luminescence during the molt
- Heatmap of trend-corrected luminescence sorted by molt 
- Plot of molt annotation over time sorted by molt


# (6) References/Citation

[Olmedo, et al. Genetics, 2015](https://www.genetics.org/content/201/2/443)  
[Meeuse, et al. Molecular Systems Biology, 2020](https://www.embopress.org/doi/full/10.15252/msb.20209498)
