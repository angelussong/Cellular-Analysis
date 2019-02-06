Instructions for running a neuron model optimization on the Neuroscience Gateway
================================================================================

This process goes from having a 'blank' neuron model with no tuned parameters to a model with optimized ion channels in 3 stages.

Stage 0) Preprocessing
    - Morphology and electrophysiology files need to be set up. The name of the neuron needs to be consistent in several places:
        - the 'CELL' parameter in 'settings.hoc'
        - a folder in the 'ephys' folder, that then contains ephys .dat files which are named 'CELL_protocol.dat'
        - a .hoc file in the 'morph' folder
        - files for errorfunctions, generators and parameters, in the appropriate subfolders of the 'setup' folder

1) Subthreshold
    a) ensure that the errorfunctions, generators and parameters files for the current cell are the subthreshold versions. Usually the suprathreshold versions of these files are present in the setup directories, and will need renaming, and then subthreshold versions of these files will need adding. They need to have the original name, which should look like: CELL_errfuncs.dat, CELL_generators.dat and CELL_params.dat
    b) Check through the settings.hoc file. Most settings can be changed, like the type of morphology, axon and ion channels to load. The following should be set to:
        - MODE = 3, to run DE
        - MULOBJ = 0, to prevent NSDE from running for subthreshold
        - max_time = 1440 (4 hours in seconds)
        - pop_size = 256 (or however many cores will be used. 256 is 16 nodes, 16 cores/node on Trestles)
        - LOAD_GENERATION = 0 (generate a new population)
        - CELL = "cellname"

    c) zip the folder up to be uploaded to NSG.
        - init.hoc should be in the root directory of the zip, so you should be one level above new_optimize
        - I usually use the command 'zip -r inputfile.zip new_optimize'
        - Recommended: copy the entire 'optimize' directory to a new location where results are going to be stored
             - something like 'cellname/subthreshold/'
             - then in that new folder, delete ephys (and morph) data for cells that aren't used here - saves a lot of size of zip file
        - Remove x86_64 folder before zipping

    d) log into NSG at nsgportal.org/portal2/home.action
        - Select the folder to run the simulation in (or Create New Folder)
        - In that folder, select Data
        - Then Upload/Enter Data
        - Type a name for the subthreshold simulation in 'Label'
        - Click 'Choose File', select inputfile.zip
        - Click 'Save' to upload the file
        
    e) Select 'Tasks', then 'Create new task'
        - Enter a description (I use the same as I named the data file)
        - Under Input, click Select Input Data and select the Data just uploaded
        - Under Tool, select 'Neuron7.3 on Comet' (or Neuron7.3 on Stampede) [actually Neuron7.4 is up there now, so could try that]
        - Under Input Parameters, select '4 parameters set'
            - Maxium Hours to Run = 4 (or whatever was set in settings.hoc)
            - Main Input Filename = init.hoc (leave it alone)
            - sub-directory name = blank (leave it alone)
            - command line options = blank (leave it alone)
            - Enter Number of Nodes = 16 (if using 256 cores on Comet)
            - Enter Number of Cores per Node = 16 (if using 256 cores on Comet) [these 2 settings need changing for different machines and population sizes]
            Click 'Save Parameters'
        - Click 'Save and Run Task', then confirm if it asks you to

    f) Wait for the run to terminate, then log back onto NSG and navigate to the correct Tasks page
        - On the task that just finished, click View Output
        - Click 'Download' on the outputfile output.tar.gz
        - Move output.tar.gz to wherever you want the results to go (usually the folder that was made prior to zipping up the inputfile), and extract
        - tar -zxf output.tar.gz

    g) Look at the output files, found in the output/ folder in optimize/
        - curr_population contains the final population of 256 individuals
        - fvevo.dat catalogues every generation of the population, with error functions values
        - de_log.dat records the best population member at every generation
            - For the subthreshold, we usually see complete convergence to one point in parameter space
            - So the final line of de_log.dat shows the best parameters found by the optimization, along with the total error associated with them
    
    h) The final step is to set the subthreshold parameters in the code ahead of running the LHS and suprathreshold optimizations
        - I have been doing this by keeping a spreadsheet of subthreshold results, then manually entering the parameter values in the .hoc file where the ion channels are setup
        - This would be in traub_channels.hoc, in the setup_parameters() method, around line 160
        - The dummy parameters are the ones that need to have their values set to the values we just found.
            - dg_pas, dr_pas, dcm are all in the first block, lines 165 ish
            - Then dgar is line 168, and dtmar and dvsar are further down at the end of the method, lines 275 and 285
        - When a model is initiated, set_params() is used which inserts these dummy parameters throughout the model, so the dummies overwrite the traub default parameters above every time the model is instantiated.
        - I usually make another copy of the whole original folder (with only the ephys for this cell), and then change the subthreshold parameters, and then I have the base folder I can use to make all of the LHS and suprathreshold runs I'm going to do with this cell - so these subthreshold parameters are never written into the main new_optimize folder, but of course this could be done potentially with a method to look at the cellname and set dummy passive parameters to different values after the cell is initialized - that would allow one main folder with all passive parameters for different neurons, while keeping the traub_channels file intact - might be better.


2) LHS
    a) Now the suprathreshold versions of the cell-specific setup files need to be replaced in the subfolders of the setup directory:
        - errorfunctions, generators, and parameters files
    
    b) The matlab scripts for performing the LHS are found in new_optimize/analysis/preprocess/
        - The setupLHS.m script will look at the number of parameters and their boundaries from the setup/parameters/CELL_params.dat file
        - One parameter in setupLHS.m that may need changing is nPerParam = 100 at the top of that file
        - run the script ( I use matlab -nodisplay -nojvm -r setupLHS from linux command line)
        - This can take a few minutes (and quite a bit longer with more parameters) due to the 10,000 samplings, although this number could probably be lowered
        - The result of the script will be a file output/curr_population, which contains all of the parameter combinations

    c) We then want to run this on the NSG: 
        - In settings.hoc, change the MODE to 4
        - zip the folder, and upload to NSGportal, as for subthreshold
        - when setting the run parameters, I still use 256 cores (16x16 on Comet), but more could be used.
        - max_time I set to 2 hours, but 1 hour is also fine and may help with queueing faster
    
    d) When the run completes, download and extract the output.tar.gz to a new directory - this will be where the suprathreshold is run from
        - This time the full result is in output/lhsresults.dat
        - go back to new_optimize/analysis/preprocess/ and check file getWeightsFromLHS.m - line 37 needs to have the correct number for top 5% to set 'ffsTop20' ( so the line should be ffsTop20 = ffsSorted( 1 : [top5%num], : ) ;
        - run the script: matlab -nodisplay -nojvm -r getWeightsFromLHS
        - efweights.dat file is written to setup/ directory
        - order of weights listed in that file is the same as the order of error functions in setup/errorfunctions.dat
        - weights from efweights.dat need inserting into the correct place in errorfunctions.dat
            - error functions have three numbers at the end - usually 1s or 0s - for the most part one of these 1s needs replacing with the weight
            - to confirm which slot the weight should be in, consult setup/generatorMould.hoc

3) Suprathreshold DE run
    Now we have error functions with weights, and we are close to running the main optimization
    
    a) alter various settings to make sure suprathreshold DE is run:
        - check settings.hoc, and change MODE to 3
        - also in settings.hoc, set DE parameters max_time and pop_size (max_time = 144000 for 40 hour DE, pop_size = 256 for 256 population and cores)
        - LOAD_GENERATION should be 0 if this is the first DE run and there is no existing curr_population
        - move output/lhsresults.dat to another folder (I move it to where I originally set up the LHS. This isn't essential, but slightly reduces the size of the zip file to be uploaded)
        
    b) check multi-objective settings:
        - in settings.hoc, set MULOBJ=1 to use NSDE, and set NSDE_COMB_FFS=0 to use every error function as a separate objective.
        - Leave MULOBJ=1, but set NSDE_COMB_FFS=1 to combine some error functions into fewer objectives
            - The number of objectives is determined by the file setup/nsde_ff_combs.dat
            - This file is in order of error functions in errorfunctions.dat, and is a list of numbers indicating which objective each error function should be included
            - So any error functions with the same number will be combined into the same objective for the NSDE

    c) empty the output folder - just delete everything in the output folder if this is the start of the DE
        - if an optimization is being continued make sure to leave curr_population in output/ (and also select LOAD_GENERATION=1 in settings.hoc)
    
    d) zip the new_optimize folder, upload to the NSG, and run the same as for subthreshold, but select correct number of hours of course..

4) Results
    a) download output.tar.gz from NSG site and extract - I extract to the same folder I set up the suprathreshold optimization in.
    b) can run ./postprocess or ./analysis/megafig/makeMegaFig.sh to generate plots
        - these will make lots of output files in the data/ folder, and put plots into a plots/ folder
        - requires octave and gnuplot v4.7 at least, or things aren't going to work
        - also, things probably won't work anyway, as the scripts are a bit wonky... so probably new analysis/plotting scripts are required...


