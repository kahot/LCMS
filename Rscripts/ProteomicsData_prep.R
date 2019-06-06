
## 1) Attach UniprotID to the Qspec Reults file

source("Rscripts/FindUniprotID.R")


# FindUniprotID("InputSourceFile.csv","Output_file_name.csv" )
# Input file should be csv file
# Example:
FindUniprotID("Data/Honokowai_HSvsHD.csv","Output/Honokowai_uniprot.csv")



## 2) Create iPath input files
source("Rscripts/iPATH_Prep.R")

# 1.Input file = files you want to run in iPath 3.
# 2.sample or location name to be used in Output files
# 3. Color for Site 1 (or Sample "0" in Qspec output)
# 4. Color for Site 2 (or Sample "1" in Qspec output)
# 5. Width of the lines

# Example
iPATH_prep("Output/Honokowai_uniprot.csv", "Honokowai", "#F68D09","#0A41E9",15)

