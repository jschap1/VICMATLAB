# Input file for plotting the results of two VIC simulations against one another

# DOES NOT WORK for my setup. The issue is in func.read.flux.files. 
# This function successfully reads in the data, but it does not output anything to compare.results()

# This works because the directory ../vis_Stehekin points at the 
#location where the function "compare results" is saved
setwd("/Users/jschapMac/Desktop/vis_Stehekin")
source(".Src/compare.results.s") 

flux_file_1 = "/Users/jschapMac/Desktop/results_VIC/fluxes_48.1875_-120.6875"
snow_file_1 = "/Users/jschapMac/Desktop/results_VIC/snow_48.1875_-120.6875"
fdepth_file_1 = ""
model_ver_1 = "4.1.0" # The version is actually 4.2.0, but the code only allows 4.0.x or 4.1.x
run_type_1  = "wb"
run_title_1 = "run_1"

flux_file_2 = "/Users/jschapMac/Desktop/results_VIC/fluxes_48.1875_-120.8125"
snow_file_2 = "/Users/jschapMac/Desktop/results_VIC/snow_48.1875_-120.8125"
fdepth_file_2 = ""
model_ver_2 = "4.1.0" 
run_type_2  = "wb"
run_title_2 = "run_2"

num_layers = 3
num_fronts = 1
# rec_interval = 24
rec_interval = "daily"
start_rec  = 1
end_rec    = 30

plot_title = "run 1 vs. run 2"
output_file_name = "/Users/jschapMac/Desktop/vis_Stehekin.ps"

compare.results(flux_file_1,snow_file_1,fdepth_file_1,model_ver_1,
                run_type_1,run_title_1,flux_file_2,snow_file_2,
                fdepth_file_2,model_ver_2,run_type_2,run_title_2,
                num_layers,num_fronts,rec_interval,start_rec,end_rec,plot_title,output_file_name)
