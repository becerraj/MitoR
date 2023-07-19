library(mitor)
message("Welcome to mitor!")
path_dir <- rstudioapi::selectDirectory("Select the patient you want to analyze: ")
mitor_analysis(path_dir)
 
