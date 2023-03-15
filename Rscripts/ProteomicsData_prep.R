
## 1) Attach UniprotID (and Defline) to the Qspec Reults file

source("Rscripts/FindUniprotID.R")


# FindUniprotID("InputSourceFile.txt","Output_file_name.csv" )
# Use the qspec output file as an input file (.txt)
# Example:
FindUniprotID("Output/OOvsON_pooled_corals.txt_qspec_paired_fdr.txt","Output/Olowalu_uniprot.csv")
FindUniprotID("Output/WNvsWO_corals.txt_qspec_paired_fdr.txt","Output/Wahikuli_uniprot.csv")
FindUniprotID("Output/qspec_results_NgvsUong.txt","Output/Palau_NgvsUl_uniprot.csv")


#run this function for all qspec result files
f<-list.files("Output/Olowalu2/", pattern= "paired_fdr$")
for (i in 1: length(f)){
        fname<-gsub(".txt_qspec_paired_fdr.txt","",f[i])
        FindUniprotID(paste0("Output/Olowalu2/",f[i]),paste0("Output/Olowalu2/",fname,".uniprot.csv"))
}


### for .csv data file ###
source("Rscripts/FindUniprotID2.R")

FindUniprotID2("Data/Plobata_tungsten.pooled.csv","Output/Plobata_tungsten_pooled_uniprot.csv")


FindUniprotID2("Data/Honokowai_QspecResults.csv","Output/BDSP/Honokowai2017_uniprot.csv")
FindUniprotID2("~/Dropbox/BDSP/Data/Hono2013_2.csv","~/Dropbox/BDSP/Data/Hono2013_2_uniprot.csv")

FindUniprotID2("Output/Rat/RAT_FSW.vs.HBP_qspec_paired_fdr.csv", "Output/Rat/RAT_FSW.vs.HBP_uniprot.csv")
FindUniprotID2("Output/Rat/HIPvsHBP.txt_qspec_paired_fdr.csv", "Output/Rat/HIP.vs.HBP_uniprot.csv")
FindUniprotID2("Output/Rat/LIPvsLBP.txt_qspec_paired_fdr.csv", "Output/Rat/LIP.vs.LBP_uniprot.csv")


## 2) Create iPath input files
source("Rscripts/iPATH_Prep.R")

# 1.Input file = files you want to run with iPath 3.
# 2.sample or location name to be used in Output files
# 3. Color for Site 1 (or Sample "0" in Qspec output)  #F68D09=Orange
# 4. Color for Site 2 (or Sample "1" in Qspec output)  #0A41E9=Blue
# 5. Width of the lines

# Example
iPATH_prep("Output/Olowalu_uniprot.csv", "Olowalu", "#F68D09","#0A41E9",15)
iPATH_prep("Output/Wahikuli_uniprot.csv", "Wahikuli", "#F68D09","#0A41E9",15)
iPATH_prep("Output/Polanui_uniprot.csv", "Polanui", "#F68D09","#0A41E9",15)
iPATH_prep("Output/Palau_NgvsUl_uniprot.csv", "Palau3", "#296318","#0000FF",15)


#for rat position data with 2 columns modify the function
iPATH_prep2("Output/Rat/RAT_FSW.vs.HBP_uniprot.csv", "FSWvsHBP", "#0000FF","#F68D09",15)
iPATH_prep2("Output/Rat/HIP.vs.HBP_uniprot.csv", "HIPvsHBP", "#0000FF","#F68D09",15)

source("Rscripts/CompGo_Prep.R")
CompGo_Prep2("Output/Rat/RAT_FSW.vs.HBP_uniprot.csv", "FSWvsHBP")
CompGo_Prep2("Output/Rat/HIP.vs.HBP_uniprot.csv", "HIPvsHBP")


