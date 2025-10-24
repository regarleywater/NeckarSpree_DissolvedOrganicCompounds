### Code to assess and compare organic compounds' targeted concentrations and non-targeted peak areas along and between the Neckar and Spree Rivers in Germany ###
### Code replicates main Hitchhikers script for calling files, blank removal and data set selection for data analysis, which is linked below
### #https://github.com/Functional-Metabolomics-Lab/FBMN-STATS/blob/main/R/.ipynb_checkpoints/Stats_Untargeted_Metabolomics-checkpoint.ipynb

# The main steps are to (1) load the libraries and (2) call the data files to be processed. 
# Subsequently (3) the non-targeted and (4) targeted data area analyzed, which entails producing
# an overview of data in the two rivers, how peak areas or concentrations vary along the rivers,  
# and how peak areas or concentrations correlate with the rivers' site characteristics. 

##### 1. LOAD LIBRARIES #####

# Global settings for plot size in the output cell:
options(repr.plot.width=10, repr.plot.height=10, res=600) # the parameters: width, height & resolution can be changed

library("tidyverse")
library("factoextra")
library("KODAMA")
library("vegan")
library("IRdisplay")
library("svglite")
library("matrixStats")
library("ggsci")
library("FSA")
library("cowplot")
library("BiocManager")
library("ComplexHeatmap")
library("dendextend")
library("NbClust")
library("xgboost")
library("caret")
library("doParallel")
library("ggrepel")
library("rfPermute")
library("stringr")
library("VennDiagram")
library("ggplot2")
library("ggvenn")
library("ggtext")
library("gridExtra")
library("viridis")
library("ggbreak")
library("kml")
library("broom")
library("data.table")
library(dplyr)
library(stringr)
library(rstatix)
library(rlang)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggridges)
library(forcats)
library(ggforce)

##### 2. CALL FILES and SET UP #####

#set directory
Directory <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Neckar and Spree Rivers Organic Pollutants Data and Code"
setwd(Directory)

# Create folders for output if not already there
dir.create(file.path(Directory, "02_Output_DataProducts"))
dir.create(file.path(Directory, "03_Output_Figures"))

# Call file with supplemental data that can go in scatter plots 
fileName_anthro <- file.path(Directory, "01_Data", "FeatureGroups_Assignment.csv")
df_csv_anthro <- read.csv(fileName_anthro,  head=TRUE, sep=",")
colnames(df_csv_anthro)

# Call file with MN results
fileName_edges <- file.path(Directory, "01_Data", "FeatureEdges_MolecularNetwork.csv")
df_csv_edges <- read.csv(fileName_edges,  head=TRUE, sep=",")
dim(df_csv_edges)

# call in csv file with Targeted Values info
fileName_TA <- file.path(Directory, "01_Data", "Targeted_Concentrations.csv")
df_csv_TAorig <- read.csv(fileName_TA,  head=TRUE, sep=",")
dim(df_csv_TAorig)

#load files from directory
file_names <- list.files('.') #list all the files in the working directory (mentioned by 'dot symbol')
file_names

#identify names of quant and metadata tables in directory
input_str <- "NonTargeted_PeakAreas_Raw.csv,Sample_Metadata.txt,FeatureAnnotations_GNPSoutput.tsv"

#read the files indicated by the indecies 
input <- (strsplit(input_str, ",")[[1]])

# reading feature table
first_line <- readLines(file.path(Directory, "01_Data", input[1]), n = 1)
if (length(strsplit(first_line, ';')[[1]]) > 1) {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F, sep = ';') # in case, ';' is the separator
} else {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F)
}

### Metadata ###
#reading metadata
md <- read.csv(file.path(Directory, "01_Data", input[2]), header = T, check.names = F, sep = '\t') # mention seperator as "/t"(tab-separated) in case of txt or tsv fi

### Annotations ###
#GNPS Annotations:
an_gnps <- read.csv(file.path(Directory, "01_Data", input[3]), header = T, check.names = F, sep = '\t') 

#explore table dimensions
dim(ft) # gives the number of rows and  columns
dim(md)
dim(an_gnps) 



#function to summarize metadata
InsideLevels <- function(metatable){
  LEVELS <- c() #creating empty vector to store information
  typ <-c()
  COUNT <- c()
  for(i in 1:ncol(metatable)){ # for each metadata column
    temp <- as.data.frame(table(metatable[,i])) #table function gives the category in each column and the count of each category
    x <- temp$Var1 #getting the name of each category in every column
    if(is.double(metatable[,i])==T){ # for numeric columns in metadata table, round the category values
      x=round(as.double(x),2)
    } 
    LEVELS <- rbind(LEVELS,toString(x)) # adding all the category values in a row
    COUNT <- rbind(COUNT,toString(temp$Freq)) # getting the frequency of each level in every column
    typ <- rbind(typ,class(metatable[,i])) # getting the class of each column
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)), #creating an output dataframe with 1st column as INDEX
                    ATTRIBUTES=colnames(metatable), #2nd column ATTRIBUTES will be the column name of metadata table
                    LEVELS, #3rd column LEVELS will give the different categories in each ATTRIBUTE
                    COUNT, #4th column COUNT will give the number of files present with each category
                    'ATTRIBUTE_CLASS'=typ, row.names=NULL) #Final column indicating the Class or datatype of each ATTTRIBUTE
  return(out)
}

#check out metadata
ncol(md) #number of columnns of metadata
InsideLevels(md[, 2:ncol(md)]) #excluding 1st filename

#check compatibility of annotations with feature table ids by checking their class
identical(class(ft$`row ID`),class(an_gnps$`#Scan#`))

#arrange annotations by Scan number
an_final <- an_gnps %>% arrange(`#Scan#`) #arranging by scan ID

#compare GNPS and analog annotations
#Function to compare between gnps_compound_name and its library match (analog)
combine_names <- function(compound_name) {
  return(paste(compound_name))
}

# POS Consolidate multiple annotations for a single '#Scan#' into one combined name
an_final_single <- an_final %>%
  group_by(`#Scan#`) %>%
  summarise(Combined_Name = combine_names(Compound_Name[1])) %>%
  ungroup() %>%
  as.data.frame()

#merge annotations with feature table
ft_an <- merge(ft, an_final_single, by.x="row ID", by.y="#Scan#", all.x= TRUE) 
head(ft_an, 2)
dim(ft_an) 


#Arrange feature table and metadata in same order
#create duplicate (working) files
new_ft <- ft #storing the files under different names to preserve the original files
new_md <- md

#clean the new files
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft)) #removing Peak area extensions from the column names of ft
new_ft <- new_ft[order(new_ft$`row ID`),,drop=F] #arranging the rows of ft file by  by ascending order of row ID
new_ft <- new_ft[,colSums(is.na(new_ft))<nrow(new_ft)] #removing if any NA columns present in the ft file,
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #removing if any NA columns present in the md file,
new_md <- new_md[apply(new_md != "", 1, any), ] # Removing rows that are completely filled with empty strings,
new_md <- new_md[, apply(new_md != "", 2, any)] # Removing columns that are completely filled with empty strings

#remove the (front & tail) spaces, if any present, from the filenames of md,
new_md$filename <- trimws(new_md$filename, which = c("both"))
rownames(new_md) <- new_md$filename
new_md <- new_md[, -which(names(new_md) == "filename")]

#update row names of feature table
if(exists("ft_an")){identical(ft_an$`row ID`,new_ft$`row ID`)} #should return TRUE if you have annotation file

#Changing the row names of the files into the combined name as "XID_mz_RT":
rownames(new_ft) <- paste(paste0("X",new_ft$`row ID`),
                          round(new_ft$`row m/z`,digits = 3),
                          round(new_ft$`row retention time`,digits = 3),
                          if(exists("ft_an")){ft_an$Combined_Name}, 
                          sep = '_') 

rownames(new_ft) <- sub("_$", "", rownames(new_ft)) #to remove the trailing underscore at rownames


#in the feature table, identify which columns correspond to samples
# Check if columns contain 'mzXML' or 'mzML' extensions
if (any(grepl('.mzML', colnames(new_ft)))) {
  # Picking only the files with column names containing 'mzXML' or 'mzML'
  new_ft <- new_ft[, grepl('.mzML', colnames(new_ft))]
  
  # Message if both .mzXML and .mzML files are present
  if (any(grepl('.mzXML$', colnames(new_ft))) && any(grepl('.mzML$', colnames(new_ft)))) {
    print("Both .mzXML and .mzML file types are present in the data")
  }
} else {
  # Ask the user for the extension if neither 'mzXML' nor 'mzML' is found
  
  your_extension <- readline('If your file extension is not .mzML or .mzXML, enter your extension (e.g., ".txt"): ')
  new_ft <- new_ft[, grepl(your_extension, colnames(new_ft))]
}

# Checking the files again to see if the above changes have been made
#head(new_ft,2)
dim(new_ft)

#head(new_md,2)
dim(new_md)

#check overlap between the feature table and metadata
new_ft <- new_ft[,order(colnames(new_ft)), drop=F] #ordering the ft by its column names
new_md <- new_md[order(rownames(new_md)),, drop=F] #ordering the md by the 1st column filename

#how many files in the metadata are also present in the feature table
table(rownames(new_md) %in% colnames(new_ft))

#are the sample names the same
identical(rownames(new_md), colnames(new_ft))

# which file names in the metadata are not in the feature table?
setdiff(rownames(new_md),colnames(new_ft))
# print(colnames(new_ft)) # uncomment to check the column names of new_ft

#checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

# Remove rows with only 0 peak area values
new_ft <- new_ft[rowSums(new_ft != 0) > 0,]
dim(new_ft)

cat("The number of rows and non-zero columns in our new ft is:",dim(new_ft))

##### 3 PROCESSING For NON-TARGETED FEATURES #####
##### 3.1 BLANK REMOVAL and SUBTRACTION #####

#transpose and merge feature table and metadata
ft_t <- as.data.frame(t(new_ft)) #transposing the ft
ft_t <- ft_t %>% mutate_all(as.numeric)  #converting all values to numeric
identical(rownames(new_md),rownames(ft_t)) #should return TRUE now

#check out the results
head(ft_t,3)

#get the index levels in your data
InsideLevels(new_md)

#enter indecies
#sample_attribute <- as.numeric(readline('Enter the index number of the attribute containing sample and blanks information: '))
sample_attribute <- as.numeric("7")
unique_sampletypes <- unique(new_md[, sample_attribute])

# Display the unique sample types along with their index
cat("Available sample types:/n")
print(unique_sampletypes)
#display(data.frame(INDEX = 1:length(unique_sampletypes), LEVELS = unique_sampletypes))

blank_ID_str <- readline('Enter the index number(s) of the blanks: ')
#blank_ID <- as.numeric(strsplit(blank_ID_str, ",")[[1]])
blank_ID <- as.numeric("1")

sample_ID_str <- readline('Enter the index number(s) of the samples: ')
# sample_ID <- as.numeric(strsplit(sample_ID_str, ",")[[1]])
sample_ID <- as.numeric("2")

# Filtering the rows from metadata with the condition = blank and sample
md_Blank <- new_md[new_md[, sample_attribute] %in% unique_sampletypes[blank_ID],]
md_Samples <- new_md[new_md[, sample_attribute] %in% unique_sampletypes[sample_ID],]

# Getting the corresponding rows from ft_t
Blank <- ft_t[which(rownames(ft_t) %in% (rownames(md_Blank))), , drop=F]
Samples <- ft_t[which(rownames(ft_t) %in% (rownames(md_Samples))), , drop=F]
InsideLevels(new_md)

#check out results
#head(Blank,n=2)
dim(Blank) 
#head(Samples, n=2)
dim(Samples)

#set Blank cutoff
#When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
#Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
Cutoff <- as.numeric(0.1)

#Getting mean for every feature in blank and Samples in a data frame named 'Avg_ft'
Avg_ft <- data.frame(Avg_blank=colMeans(Blank, na.rm= F)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
# Avg_ft$Max_blank <- apply(Blank, 2, max, na.rm = F) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ft$Avg_samples <- colMeans(Samples, na.rm= F) # adding another column 'Avg_samples' for feature means of samples

#Getting the ratio of blank vs Sample
Avg_ft$Ratio_blank_Sample <- (Avg_ft$Avg_blank+1)/(Avg_ft$Avg_samples+1)

# Creating a bin with 1s when the ratio>Cutoff, else put 0s
Avg_ft$Bg_bin <- ifelse(Avg_ft$Ratio_blank_Sample > Cutoff, 1, 0 )

#Calculating the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ft)))
print(paste("No.of Background or noise features:",sum(Avg_ft$`Bg_bin` ==1,na.rm = T)))
print(paste("No.of features after excluding noise:",(ncol(Samples) - sum(Avg_ft$`Bg_bin` ==1,na.rm = T))))

blk_rem_1 <- merge(as.data.frame(t(Samples)), Avg_ft, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) %>% #removing the last 4 columns
  column_to_rownames(var="Row.names") 
  blk_rem <- as.data.frame(t(blk_rem_1))
  
        # Blank subtraction step (using max value observed in the blanks) ##
        # Add Max_blank column to Avg_ft
        Avg_ft$Max_blank <- apply(t(Blank), 1, max, na.rm = TRUE)
        
        # Correctly match feature names using row names from Avg_ft
        Blank_max_filtered <- Avg_ft[names(blk_rem), "Max_blank", drop = FALSE]$Max_blank
        
        # Subtract Max blank values from Samples, ensuring no negative values
        blk_rem_sub <- sweep(blk_rem, 2, Blank_max_filtered, FUN = "-")
        blk_rem_sub[blk_rem_sub < 0] <- 0  # Set negative values to zero
        
        # Determine the limit of detection (LOD)
        Cutoff_LOD <- round(min(blk_rem[blk_rem > 0], na.rm = TRUE))
        print(paste0("The limit of detection (LOD) is: ", Cutoff_LOD))
        
        # Apply LOD correction
        blk_rem_sub[blk_rem_sub > 0 & blk_rem_sub < Cutoff_LOD] <- 0
        
        # View the final processed dataset
        print(dim(blk_rem_sub))
        
        # Print results
        # Define output file path
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_BlanksRemoved_and_MaxSubtracted_with_cutoff_',Cutoff,'.csv'))
        # Write the CSV file
        write.csv(t(blk_rem_sub), output_path, row.names = TRUE)

        
# Ensure the final dataframe is formatted correctly
blk_rem_sub <- as.data.frame(blk_rem_sub)

# review dimensions of feature table wihtout blanks
dim(blk_rem_sub)

# review dimensions of metadata without the blanks info 
dim(md_Samples)


##### 3.2 PEAK AREA REPLICATE AVERAGING and SCALING #####

# make copies of data tables
blk_rem_use <- blk_rem_sub
md_Samples_use <- md_Samples

# Create a function to calculate the average with the specified condition
average_with_zero_condition <- function(values) {
  if(any(values == 0)) {
    return(0)
  } else {
    return(mean(values))
  }
}

# Create a new column 'sample_id' by removing the replicate information from row names
blk_rem_new <- blk_rem_use %>%
  rownames_to_column("sample") %>%
  mutate(sample_id = gsub("(_both_[12]).*", "_both_1.mzML", sample)) 

# Group by 'sample_id' and calculate the average for each feature
blk_avg_pa <- blk_rem_new %>%
  group_by(sample_id) %>%
  summarise(across(!starts_with("sample"), average_with_zero_condition)) %>%
  column_to_rownames("sample_id")
dim(blk_avg_pa)

# remove zero-variance columns
blk_avg_pa <- blk_avg_pa[, apply(blk_avg_pa, 2, sd) != 0]

dim(blk_avg_pa)

# Define output file path
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged.csv'))

# Write the CSV file
write.csv(t(blk_avg_pa), output_path, row.names = TRUE)


## Scaling
# feature max norm function
normalize <- function(x) {
  return(x / max(x))
}

# Applying the normalization function to each column
blk_avg_fm <- as.data.frame(apply(blk_avg_pa, 2, normalize))

# Define output file path
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged_Scaled.csv'))

# Write the CSV file
write.csv(t(blk_avg_fm), output_path, row.names = TRUE)

# Apply normalization to all columns except "Erpe_both_1.mzML" and drop that column
blk_avg_fm_noErpe <- as.data.frame(
  apply(blk_avg_pa[rownames(blk_avg_pa) != c("Erpe_both_1.mzML"),], 2, normalize))

# Define output file path
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged_Scaled_noErpe.csv'))

# Write the CSV file
write.csv(t(blk_avg_fm_noErpe), output_path, row.names = TRUE)


##### 3.3 SUBSET ANTHROPOGENIC FEATURES USING ANNOTATIONS and MOLECULAR NETWORKS#####

library(dplyr)
library(igraph)

# Create a copy of the original feature group table
features_use <- df_csv_anthro

# Define the group columns
group_cols <- paste0("Group", 0:2)

# Find the rows where cosine is less than 0.8 and the MZerror is greater than 10 - the annotations will not be used from these features
rows_to_blank <- which(!is.na(features_use$Cosine) & features_use$Cosine < 0.8 | !is.na(features_use$Mzerror) & features_use$Mzerror >= 10)
print(rows_to_blank)

# Set the group column values to blank if a feature's cosine or Mzerror values do not meet the constraints
features_use[rows_to_blank, group_cols] <- ""

# View result
dim(features_use)
names(features_use)
str(features_use[group_cols])
length(rows_to_blank)

# Filter the features based on if they have an annotation category
features_filtered <- features_use %>%
  filter(Group0 %in% c("all anthropogenic"))

# Extract FeatureIDs present in the feature table after blank removal and subtraction, replicate averaging and scaling
blk_feature_ids <- gsub(".*X([0-9]+).*$", "\\1", colnames(blk_avg_fm))

# Filter df_csv_edges to keep only rows where BOTH nodes are in blk_rem
df_csv_edges_use <- df_csv_edges %>%
  filter(node1 %in% blk_feature_ids & node2 %in% blk_feature_ids)
dim(df_csv_edges_use)

# Initialize all features df
all_features_df <- data.frame(
  FeatureID = character(),
  Parent = character(),
  Source = character(),
  stringsAsFactors = FALSE
)

# Add filtered features with their subsets
all_features_df <- features_filtered %>%
  select(FeatureID, Group0) %>%
  rename(Source = Group0) %>%
  mutate(Parent = FeatureID)

# For each anthro feature, get first neighbors
for (feature in features_filtered$FeatureID) {

  connected_node2 <- df_csv_edges_use %>%
    filter(node1 == feature) %>%
    select(node2) %>%
    rename(FeatureID = node2) %>%
    mutate(Source = "MN", Parent = feature)

  connected_node1 <- df_csv_edges_use %>%
    filter(node2 == feature) %>%
    select(node1) %>%
    rename(FeatureID = node1) %>%
    mutate(Source = "MN", Parent = feature)

  connected_features <- bind_rows(connected_node1, connected_node2)
  all_features_df <- bind_rows(all_features_df, connected_features)
}

# Collect all IDs found so far
combined_feature_ids <- all_features_df$FeatureID

# Create new df to use
all_features_df_1 <- all_features_df 
# %>%
#   merge(features_use[, c("FeatureID", "Group0", "Group1", "Group2")], by.x = "Parent", by.y = "FeatureID", all.x = TRUE)

### Use igraph to add nth MN and cluster info ###

# Ensure only 2 columns are passed to define edges
edge_df <- df_csv_edges_use[, c("node1", "node2")]

# Build the graph cleanly
g <- graph_from_data_frame(edge_df, directed = FALSE)

# Check vertex names
head(V(g)$name)

# Get anthro features from the combined df
anthro_ids <- all_features_df_1 %>%
  filter(Source %in% c("all anthropogenic")) %>%
  pull(FeatureID) %>%
  intersect(V(g)$name)  # make sure they exist in graph

# Get component membership
components_info <- components(g)
node_to_component <- components_info$membership

# Add cluster ID to all_features_df_1
all_features_df_1$ClusterID <- node_to_component[match(all_features_df_1$FeatureID, names(node_to_component))]

# Identify which clusters contain anthro features
anthro_clusters <- unique(all_features_df_1$ClusterID[all_features_df_1$Source == "all anthropogenic"])

# Now find nodes in those anthro clusters
nodes_in_anthro_clusters <- names(node_to_component)[node_to_component %in% anthro_clusters]

# Find nth MN: in anthro clusters, but not anthro or 1st neighbors
nth_mn <- setdiff(nodes_in_anthro_clusters, all_features_df_1$FeatureID[all_features_df_1$Source %in% c("all anthropogenic", "MN")])

# Add these to your dataframe
nth_mn_df <- data.frame(
  FeatureID = nth_mn,
  Source = "nthMN",
  Parent = nth_mn
)

# Make sure the FeatureID columns are the same type.
all_features_df_1$FeatureID <- as.numeric(all_features_df_1$FeatureID)
all_features_df_1$Parent <- as.numeric(all_features_df_1$Parent)
features_use$FeatureID <- as.numeric(features_use$FeatureID)
nth_mn_df$FeatureID <- as.numeric(nth_mn_df$FeatureID)
nth_mn_df$Parent <- as.numeric(nth_mn_df$Parent)

# Merge all features together and merge with group info.
all_features_df_2 <- bind_rows(all_features_df_1, nth_mn_df)

# Collect all features into one df with their group info
combined_df <- full_join(all_features_df_2,
                         features_use[, c("FeatureID", "Feature", "Group0", "Group1", "Group2")],
                         by = "FeatureID") 

# Create a lookup table for parent feature group info
parent_groups <- features_use %>%
  select(FeatureID, Group0, Group1, Group2) %>%
  rename(Parent = FeatureID,
         Group0_parent = Group0,
         Group1_parent = Group1,
         Group2_parent = Group2)

# Join parent info to combined_df
combined_df_joined <- combined_df %>%
  left_join(parent_groups, by = "Parent")

# Replace Group0–2 values with parent values if Source == "MN"
all_features_df_use <- combined_df_joined %>%
  mutate(
    Group0 = ifelse(Source == "MN" & !is.na(Group0_parent), Group0_parent, Group0),
    Group1 = ifelse(Source == "MN" & !is.na(Group1_parent), Group1_parent, Group1),
    Group2 = ifelse(Source == "MN" & !is.na(Group2_parent), Group2_parent, Group2)
  ) %>%
  select(-Group0_parent, -Group1_parent, -Group2_parent)


# Add cluster ID to all_features_df_use
all_features_df_use$ClusterID <- node_to_component[match(all_features_df_use$FeatureID, names(node_to_component))]

# Label all others explicitly as "unknown" if still NA
# all_features_df_use$Group0[is.na(all_features_df_use$Group0)] <- "u"
all_features_df_use$Group0[all_features_df_use$Group0 == ""] <- "unknown"
all_features_df_use$Group1[all_features_df_use$Group1 == ""] <- "unknown"
all_features_df_use$Group2[all_features_df_use$Group2 == ""] <- "unknown"
all_features_df_use$Group1[all_features_df_use$Group0 == "biological"] <- "biological compound"
all_features_df_use$Group2[all_features_df_use$Group0 == "biological"] <- "biological compound"

# Add priority
print(unique(all_features_df_use$Source))
all_features_df_use <- all_features_df_use %>%
  mutate(Priority = case_when(
    Source == "all anthropogenic" ~ 1,
    Source == "MN" ~ 2,
    Source == "nthMN" ~ 3
  )) %>%
  arrange(FeatureID, Priority) %>%
  distinct(FeatureID, .keep_all = TRUE) %>%
  select(-Priority)

# Count the occurrences of each bin using table()
summary <- table(all_features_df_use$Source)
print(summary)

# For Source = "all anthropogenic"
summary_a <- table(all_features_df_use$Group1[all_features_df_use$Source == "all anthropogenic"])
print(summary_a)

# For Source = "MN"
summary_MN <- table(all_features_df_use$Group1[all_features_df_use$Source == "MN"])
print(summary_MN)


# Write results
# Define output file path
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Feature_Annotations_MNneighbors_and_Groups.csv"))

# Write the CSV file
write.csv(all_features_df_use, output_path, row.names = TRUE)

                                                          
##### 3.4 DATA OVERVIEW and RIVER COMPARISON (FIGURE 2) #####                 
                      
long_quant <- blk_avg_fm %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Feature", values_to = "PeakArea")

md_Samples_use$Sample <- rownames(md_Samples_use)

# Extract FeatureID from column names (assuming X123_abcd format) FOR NON-TARGETED
blk_avg_md <- long_quant %>%
  mutate(FeatureID = as.numeric(gsub(".*X([0-9]+).*", "\\1", Feature))) %>%
  left_join(all_features_df_use %>% select(FeatureID, Source, Group0, Group1, Group2), by = "FeatureID")  %>%
  left_join(md_Samples_use, by = "Sample")

blk_avg_md$Group0[is.na(blk_avg_md$Group0)] <- "unknown"
blk_avg_md$Group1[is.na(blk_avg_md$Group1)] <- "unknown"
blk_avg_md$Group2[is.na(blk_avg_md$Group2)] <- "unknown"

# Assuming your dataframe is called df
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggvenn)
library(ggpubr)

# Define color palette for Group1
group1_colors <- c(
  "biological compound" = "#1B9E77",
  "other consumer industrial" = "#7570B3",
  "pesticide"  = "#66A61E",
  "pharmaceutical" = "#377EB8",
  "polymer additive" = "#FF7F00",
  "unknown" = "#999999"
)

plot_pie_venn_box <- function(data, source_filter, label, Directory) {
  
  # Subset dataframe
  df <- data %>%
    filter(!!rlang::enquo(source_filter))  # flexible filtering expression
  
  # --- PIE CHART by Group1 ---
  pie_df <- df %>%
    distinct(FeatureID, Group1) %>%          # unique features per category
    group_by(Group1) %>%
    summarise(FeatureCount = n(), .groups = "drop") %>%
    mutate(Percent = FeatureCount / sum(FeatureCount) * 100,
           Label = paste0(Group1, "\n", FeatureCount, " (", round(Percent, 1), "%)"))
  
  pie_plot <- ggplot(pie_df, aes(x = "", y = FeatureCount, fill = Group1)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.5),
              size = 3.5, color = "black", lineheight = 0.9) +
    scale_fill_manual(values = group1_colors, na.value = "#999999") +
    theme_void() +
    labs(title = paste("Feature distribution by Group1 (", label, ")", sep = "")) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
    )
  
  print(pie_plot)
  
  # Save pie chart
  pie_out <- file.path(Directory, "03_Output_Figures",
                       paste0(Sys.Date(), "_PieChart_", label, "_Group1.pdf"))
  pdf(pie_out, height = 4, width = 4)
  print(pie_plot)
  dev.off()
  
  # --- VENN DIAGRAM ---
  venn_list_df <- df %>%
    filter(PeakArea > 0) %>%
    distinct(FeatureID, ATTRIBUTE_RiverGroup) %>%
    group_by(ATTRIBUTE_RiverGroup) %>%
    summarise(Features = list(unique(FeatureID)), .groups = "drop") %>%
    pivot_wider(names_from = ATTRIBUTE_RiverGroup, values_from = Features)
  
  set1 <- venn_list_df$Spree[[1]]
  set2 <- venn_list_df$Neckar[[1]]
  
  venn_data <- list("Spree" = set1, "Neckar" = set2)
  
  venn_plot <- ggvenn(
    venn_data,
    show_percentage = TRUE,
    fill_color = c('#8795e8', '#60e4c1', "#868686FF"),
    stroke_size = 0.5, set_name_size = 0, stroke_color = "white",
    auto_scale = TRUE, text_size = 7
  )
  
  print(venn_plot)
  
  # Save venn diagram
  venn_out <- file.path(Directory, "03_Output_Figures",
                        paste0(Sys.Date(), "_Venn_", label, ".pdf"))
  
  pdf(file = venn_out, height = 3, width = 3)
  print(venn_plot)
  dev.off()
  
  # --- BOX PLOT ---
  summed_df <- df %>%
    group_by(Sample, ATTRIBUTE_RiverGroup) %>%
    summarise(SummedPeakArea = mean(PeakArea, na.rm = TRUE), .groups = "drop") %>%
    mutate(ATTRIBUTE_RiverGroup = factor(ATTRIBUTE_RiverGroup,
                                         levels = c("Spree", "Neckar")))
  
  # Wilcoxon test
  test_results <- summed_df %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(SummedPeakArea ~ ATTRIBUTE_RiverGroup)$p.value,
        error = function(e) NA
      )
    )
  
  print(test_results)
  
  # Boxplot
  boxplot <- ggplot(summed_df, aes(x = ATTRIBUTE_RiverGroup, y = SummedPeakArea,
                                   fill = ATTRIBUTE_RiverGroup)) +
    geom_boxplot(position = position_dodge(width = 0.5)) +
    theme_bw() +
    scale_fill_manual(values = c("Neckar" = '#60e4c1', "Spree" = '#8795e8')) +
    labs(x = "River", y = "Average Scaled Peak Area", fill = "River") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  print(boxplot)
  
  # Save boxplot
  boxplot_out <- file.path(Directory, "03_Output_Figures",
                           paste0(Sys.Date(), "_Boxplots_", label, "_PeakAreas_byRiver.pdf"))
  
  pdf(file = boxplot_out, height = 3, width = 5)
  print(boxplot)
  dev.off()
  
  return(test_results)
}

# Anthropogenic features
plot_pie_venn_box(
  data = blk_avg_md,
  source_filter = Source %in% c("all anthropogenic", "MN"),
  label = "AnthropogenicFeatures",
  Directory = Directory
)

# Other features
plot_pie_venn_box(
  data = blk_avg_md,
  source_filter = !Source %in% c("all anthropogenic", "MN"),
  label = "OtherFeatures",
  Directory = Directory
)




##### 3.5 HEATMAP of FEATURE GROUPS across SITES (FIGURE 3) #####

# Prepare the data by grouping features 
        # Convert scaled peak areas dataframe to long format
        long_quant <- blk_avg_fm_noErpe %>%
          tibble::rownames_to_column("Sample") %>%
          pivot_longer(-Sample, names_to = "Feature", values_to = "PeakArea")

        # Extract FeatureID from column names (assuming X123_abcd format)
        long_quant <- long_quant %>%
          mutate(FeatureID = as.numeric(gsub(".*X([0-9]+).*", "\\1", Feature)))

        # Join with metadata
        quant_meta <- long_quant %>%
          left_join(all_features_df_use %>% select(FeatureID, Group0, Group1, Group2), by = "FeatureID")

        # Set any remaining NAs to "unknown"
        quant_meta$Group0[is.na(quant_meta$Group0)] <- "unknown"
        quant_meta$Group1[is.na(quant_meta$Group1)] <- "unknown"
        quant_meta$Group2[is.na(quant_meta$Group1)] <- "unknown"

        # Calculate ALL features sums per sample
        groupALL_contributions <- quant_meta %>%
          group_by(Sample) %>%
          summarise(GroupPeak = mean(PeakArea, na.rm = TRUE)) %>%
          mutate(Group = "all") %>%
          relocate(Group, .after = Sample)

        # Calculate OTHER features sums per sample
        groupOTHER_contributions <- quant_meta %>%
          filter(Group0 %in% c("biological", "unknown")) %>%
          group_by(Sample) %>%
          summarise(GroupPeak = mean(PeakArea, na.rm = TRUE)) %>%
          mutate(Group = "other") %>%
          relocate(Group, .after = Sample)

        # Calculate Group0 sums per sample
        group0_contributions <- quant_meta %>%
          group_by(Sample, Group0) %>%
          summarise(GroupPeak = mean(PeakArea, na.rm = TRUE)) %>%
          rename(Group = Group0)

        # Calculate total and group2 sums per sample
        group1_contributions <- quant_meta %>%
          group_by(Sample, Group1) %>%
          summarise(GroupPeak = mean(PeakArea, na.rm = TRUE)) %>%
          rename(Group = Group1)

        # Calculate Group2 sums per sample
        group2_contributions <- quant_meta %>%
          group_by(Sample, Group2) %>%
          summarise(GroupPeak = mean(PeakArea, na.rm = TRUE)) %>%
          rename(Group = Group2)

# Heatmap

        # Combine the grouped feature tables
        corr_df <- rbind(groupALL_contributions, groupOTHER_contributions, group0_contributions, group1_contributions, group2_contributions)

        # Add Sample column from rownames
        md_Samples_use$Sample <- rownames(md_Samples_use)

        # Join grouped feature table values with sample metadata
        corr_df <- corr_df %>%
          left_join(md_Samples_use, by = "Sample")

        # Scale the grouped feature peak areas by the max peak area observed for that group
        corr_df <- corr_df %>%
          group_by(Group) %>%
          mutate(GroupPeak_norm = GroupPeak / max(GroupPeak, na.rm = TRUE)) %>%
          ungroup()  %>%
          relocate(GroupPeak_norm, .after = GroupPeak)

        # Choose names of site characteristics to correlate feature groups with
        colnames(md_Samples_use)
        meta_vars <- c( "WWTP_TWW_m3.Year.Km2", "population_density_km2",
                        "Urban fabric","Industrial or commercial", "Road rail and Trnsport" ,
                       "Agricultural areas", "Forests","Water bodies")

        # Choose groups that you would like to plot
        unique(corr_df$Group)
        groups <- c("all anthropogenic",
                    "pharmaceutical",
                    "organophosphate","HMMM+", "PAG",
                    "other consumer industrial",
                    "fungicide", "herbicide", "insecticide",
                    "biological compound", "other",
                    "all"
        )

        # Set groups levels as the reverse of groups (for plotting orientation)
        groups_levels <- rev(groups)

        # Filter grouped features dataframe to contain only samples from the main stem
        corr_df_filtered <- corr_df %>%
        filter(Group %in% groups) %>%
        # filter(Sample %in% c("N01_both_1.mzML",     "N02_both_1.mzML",     "N03_both_1.mzML",     "N04_both_1.mzML",
        #                      "N05_both_1.mzML",     "N06_both_1.mzML",     "N07_both_1.mzML",     "N08_both_1.mzML",     "N09_both_1.mzML",
        #                      "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
        #                      "S01_both_1.mzML",     "S08_both_1.mzML",     "S11_both_1.mzML",     "S13_both_1.mzML",     "S18_both_1.mzML",
        #                      "S26_both_1.mzML",     "S29_both_1.mzML",     "S51_both_1.mzML",
        #                      "S56_both_1.mzML",     "S57_both_1.mzML",     "S64_both_1.mzML",     "S67_both_1.mzML",
        #                      "S69_both_1.mzML",     "S72_both_1.mzML",     "S73_both_1.mzML",
        #                      "S74_both_1.mzML",     "S75_both_1.mzML",     "S76_both_1.mzML",     "S76a_both_1.mzML"))
        # 
        filter(Sample %in% c("N1_both_1.mzML",     "N2_both_1.mzML",     "N3_both_1.mzML",     "N4_both_1.mzML",
                             "N5_both_1.mzML",     "N6_both_1.mzML",     "N7_both_1.mzML",     "N8_both_1.mzML",     "N9_both_1.mzML",
                             "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
                             "S1_both_1.mzML",     "S2_both_1.mzML",     "S3_both_1.mzML",     "S4_both_1.mzML",     "S5_both_1.mzML",
                             "S6_both_1.mzML",     "S7_both_1.mzML",     "S8_both_1.mzML",
                             "S9_both_1.mzML",     "S10_both_1.mzML",     "S11_both_1.mzML",     "S12_both_1.mzML",
                             "S13_both_1.mzML",     "S14_both_1.mzML",     "S15_both_1.mzML",
                             "S16_both_1.mzML",     "S17_both_1.mzML",     "S18_both_1.mzML",     "S19_both_1.mzML"))

      length(unique(corr_df_filtered$Group))

      # Pivot to wide format for ComplexHeatmap
      heatmap_matrix <- corr_df_filtered %>%
        group_by(Group, Sample) %>%
        summarise(GroupPeak_norm = mean(GroupPeak_norm, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Sample, values_from = GroupPeak_norm) %>%
        column_to_rownames("Group") %>%
        as.matrix()

      # Ensure rownames are a factor in the correct order
      heatmap_matrix_ordered <- heatmap_matrix[groups, , drop = FALSE]

      library(ComplexHeatmap)
      library(circlize)

      set.seed(1235)
      hmap <- Heatmap(
        heatmap_matrix_ordered,  # or use heatmap_matrix if unscaled
        heatmap_legend_param = list(title = "Scaled\nGroupPeak"),
        col = circlize::colorRamp2(c(0, 0.0001,0.5, 1), c("darkgray", "blue", "white", "red")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        width = unit(500, "mm"),
        height = unit(300, "mm"),
        row_km = 1,
        row_km_repeats = 100
      )

      hk <- ComplexHeatmap::draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right")

      # Save plot
      output_path <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Heatmap_AveragedandScaled_GroupedFeatures_PeakAreas_MainRiverSites_noErpe.pdf"))

      pdf(file = output_path, height = 30, width = 30)
      print(hk)
      dev.off()

##### 3.6 CORRELATIONS of FEATURE GROUPS with SITE CHARACTERISTICS (FIGURE 4) #####

# Run correlations between feature groups and site characteristics

      # Check presence and type
      str(corr_df[meta_vars])

      # Loop through metadata variables and calculate correlations
      cor_summary <- map_dfr(meta_vars, function(var) {
        corr_df %>%
          filter(Group %in% groups) %>%
          group_by(Group, ATTRIBUTE_RiverGroup) %>%
          summarise(
            MetadataVar = var,
            Correlation = cor(GroupPeak, .data[[var]], use = "complete.obs", method = "spearman"),
            Pvalue = cor.test(GroupPeak, .data[[var]], method = "spearman")$p.value,
            .groups = "drop"
          )
      })


      # Order MetadataVar factor levels for plotting
      cor_summary <- cor_summary %>%
        mutate(MetadataVar = factor(MetadataVar, levels = meta_vars),
               Group = factor(Group, levels = groups_levels))


      cor_summary <- cor_summary %>%
        mutate(Signif = case_when(
          Pvalue < 0.001 ~ "***",
          Pvalue < 0.01 ~ "**",
          Pvalue < 0.05 ~ "*",
          TRUE ~ ""
        ))

      # Define output file path
      output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Correlation_GroupedFeatsVsSiteChar_SpearmanRho_AllSites_noErpe.csv"))

      # Write the CSV file
      write.csv(cor_summary, output_path, row.names = TRUE)

      # Plot correlations in Neckar
      neckar <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Neckar" & !is.na(cor_summary$Group),], aes(x = MetadataVar, y = Group)) +
        geom_point(aes(size = abs(Correlation), color = Correlation)) +
        geom_text(aes(label = Signif), size = 5, vjust = 0.8) + #vjust = -1.2
        scale_color_gradient2(
          low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
          name = "spearman r"
        ) +
        scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
        theme_bw() +
        labs(
          x = "Metadata Variable",
          y = "Feature Group (Group)"
        ) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "right"
        )

        # Plot correlations in Spree
        spree <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Spree"& !is.na(cor_summary$Group),], aes(x = MetadataVar, y = Group)) +
              geom_point(aes(size = abs(Correlation), color = Correlation)) +
              geom_text(aes(label = Signif), size = 5, color = "black", vjust = 0.8) +
              scale_color_gradient2(
                low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
                name = "spearman r"
              ) +
              scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
              theme_bw() +
              labs(
                # title = "Correlation Between Group PA and Metadata",
                x = "Metadata Variable",
                y = "Feature Group (Group)"
              ) +
              theme(
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
                legend.position = "none"
              )

        # Save correlation results figures together
        output_path <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Correlation_GroupedFeatsVsSiteChar_SpearmanRho_AllSites_noErpe.pdf"))

        pdf(file = output_path, height = 7, width = 13)
        plot_grid( spree,  neckar, ncol = 2, nrow = 1, rel_widths = c(1,1), rel_heights = c(1,1,1,1), label_size = 24, label_y = 0.9,
                   align = "v")
        dev.off()

# Correlate feature groups with each other
        library(dplyr)
        library(tidyr)
        library(purrr)

        # Make sure character
        corr_df_use <- corr_df %>%
          mutate(across(c(Sample, ATTRIBUTE_RiverGroup, Group), as.character))

        # Round to numeric to 4 decimal
        corr_df_use <- corr_df_use %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))

        # Pivot data wider: one column per Group
        corr_wide <- corr_df_use %>%
          filter(Group %in% groups) %>%
          select(Sample, ATTRIBUTE_RiverGroup, Group, GroupPeak) %>%
          pivot_wider(names_from = Group, values_from = GroupPeak)

        # Rename the column in corr_wide
        corr_wide <- corr_wide %>%
          rename(OtherConInd = `other consumer industrial`)

        # Extract the first element from each list in OtherConInd
        corr_wide <- corr_wide %>%
          mutate(OtherConInd = map_dbl(OtherConInd, ~ .x[1]))

        # Ensure group columns are numeric
        corr_wide <- corr_wide %>%
          mutate(across(!c(Sample, ATTRIBUTE_RiverGroup), as.numeric))

        # Function to compute correlations between groups
        get_cor_summary <- function(df, ATTRIBUTE_RiverGroup) {
          groups_cols <- setdiff(names(df), c("Sample", "ATTRIBUTE_RiverGroup"))

          # All unique group pairs
          pairs <- combn(groups_cols, 2, simplify = FALSE)

          map_dfr(pairs, function(pair) {
            x <- df[[pair[1]]]
            y <- df[[pair[2]]]

            ct <- cor.test(x, y, method = "spearman")

            tibble(
              Group1 = pair[1],
              Group2 = pair[2],
              ATTRIBUTE_RiverGroup = ATTRIBUTE_RiverGroup,
              Correlation = unname(ct$estimate),
              Pvalue = ct$p.value
            )
          })
        }

        # Run per ATTRIBUTE_RiverGroup
        cor_summary <- corr_wide %>%
          group_split(ATTRIBUTE_RiverGroup) %>%
          map_dfr(~ get_cor_summary(.x, unique(.x$ATTRIBUTE_RiverGroup)))

        cor_summary <- cor_summary %>%
          mutate(Signif = case_when(
            Pvalue < 0.001 ~ "***",
            Pvalue < 0.01 ~ "**",
            Pvalue < 0.05 ~ "*",
            TRUE ~ ""
          ))

        # Define output file path
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Correlation_BetweenFeatureGroups_SpearmanRho_AllSites_noErpe.csv"))

        # Write the CSV file
        write.csv(cor_summary, output_path, row.names = TRUE)

        # Plot correlations in Neckar
        neckar <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Neckar",], aes(x = Group1, y = Group2)) +
          geom_point(aes(size = abs(Correlation), color = Correlation)) +
          geom_text(aes(label = Signif), size = 5, vjust = 0.8) + #vjust = -1.2
          scale_color_gradient2(
            low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
            name = "spearman r"
          ) +
          scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
          theme_bw() +
          labs(
            # title = "Correlation Between Group PA and Metadata",
            x = "Feature Group",
            y = "Feature Group"
          ) +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "right"
          )

        # Plot correlations in Spree
        spree <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Spree",], aes(x = Group1, y = Group2)) +
          geom_point(aes(size = abs(Correlation), color = Correlation)) +
          geom_text(aes(label = Signif), size = 5, color = "black", vjust = 0.8) +
          scale_color_gradient2(
            low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
            name = "spearman r"
          ) +
          scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
          theme_bw() +
          labs(
            # title = "Correlation Between Group PA and Metadata",
            x = "Feature Group",
            y = "Feature Group"
          ) +
          theme(
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
            legend.position = "none"
          )

        # Save correlation results figures together
        output_path <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Correlation_BetweenFeatureGroups_SpearmanRho_AllSites_noErpe.pdf"))

        pdf(file = output_path, height = 7, width = 13)
        plot_grid(spree,  neckar, ncol = 2, nrow = 1, rel_widths = c(1,1), rel_heights = c(1,1,1,1), label_size = 24, label_y = 0.9,
                   align = "v")
        dev.off()


                                                                                  
                                                                                  
##### 4 PROCESSING For TARGETED COMPOUNDS #####
##### 4.1 BLANK REMOVAL and SUBTRACTION #####
      
    # Data used for targeted data analysis has already been subjected
    # to blank subtraction in a previous, manual process. 
        
##### 4.2 CONCENTRATION SCALING #####
        
        # make copies of data tables
        blk_rem_use <- df_csv_TAorig
        md_Samples_use <- md_Samples
        
        # Format for data processing and analysis
        blk_rem_new <- blk_rem_use %>%
          rename(sample = row.ID) %>%
          mutate(sample_id = paste(sample, "_both_1.mzML", sep = "")) %>%
          column_to_rownames("sample_id") %>% 
          select(-sample)
        
        # Format for data processing and analysis
        blk_rem_new_noErpe <- blk_rem_use %>%
          rename(sample = row.ID) %>%
          mutate(sample_id = paste(sample, "_both_1.mzML", sep = "")) %>%
          filter(!sample_id %in% c("Erpe_both_1.mzML")) %>%
          column_to_rownames("sample_id") %>% 
          select(-sample)
        
        ## Scaling
        # feature max norm function
        normalize <- function(x) {
          return(x / max(x))
        }
        
        # Applying the normalization function to each column
        blk_avg_fm <- as.data.frame(apply(blk_rem_new, 2, normalize))
        
        # Define output file path
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_Concentrations_Scaled.csv'))
        
        # Write the CSV file
        write.csv(blk_avg_fm, output_path, row.names = TRUE)
        
        # Apply normalization to all columns except "Erpe_both_1.mzML" and drop that column
        blk_avg_fm_noErpe <- as.data.frame(apply(blk_rem_new_noErpe, 2, normalize))
        
        # Define output file path
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_Concentrations_Scaled_noErpe.csv'))
        
        # Write the CSV file
        write.csv(blk_avg_fm_noErpe, output_path, row.names = TRUE)
        
##### 4.3 DATA OVERVIEW and RIVER COMPARISON (FIGURE 2) #####
      
      # Set up dataframe                                                             
      long_quant <- blk_avg_fm %>%
        tibble::rownames_to_column("Sample") %>%
        pivot_longer(-Sample, names_to = "Compound", values_to = "Concentration")
      
      # Make sure Sample column is in metadata
      md_Samples_use$Sample <- rownames(md_Samples_use)
      
      # Join metadata and quant data together
      blk_avg_md <- long_quant %>%
        mutate(FeatureID = Compound) %>%
        left_join(md_Samples_use, by = "Sample") 
      
      # Make sure sample names are in a proper column if they are in rownames
      df <- blk_avg_md
      unique(df$Compound)
      
      # Establish prefered order of compounds
      Compound_order <- rev(c("Formyaminoantipyrin", 
                             "o.Desmethylvenlafaxine", "Venlafaxine", 
                             "Carbamazepine",  "Gabapentin.Lactam", "Primidone", 
                             "Irbesartan","Losartan", "Olmesartan","Valsartan", "Valsartansaeure",  "Metoprolol.Tartrate", "Metoprolol.acid",
                             "Sulfamethoxazol",
                             "Caffeine",
                             "Bezafibrat", "Hydrochlorothiazid", "Metformin", "N.Guanylurea",
                             "Triphenylphosphate", "Tris.2.chloroetyl.phosphate", "Tris.2.chloroisopropyl.phosphate", 
                             "Hexa.methoxymethyl.melamine", 
                             "X1.3.Diphenylguanidine", 
                             "Benzothiazol.2.sulfonic.acid", "X2.Hydroxybenzothiazole",  
                             "Benzotriazol","X5.Methyl.benzotriazol",
                             "N.Butylbenzesulfonamide",  "Tributylamine",
                             "DEET", 
                             "Diuron", "Flufenacet", "X5.Metolachlor",  "Terbutryn", 
                             "Imidacloprid" ))
      print(Compound_order)
      
      # Update compound names
      Compound_rename <- c(
        "Tributylamine" = "tributylamine",
        "Diuron" = "diuron",
        "Flufenacet" = "flufenacet",
        "Metoprolol.acid" = "metoprolol acid",
        "Triphenylphosphate" = "triphenylphosphate",
        "Imidacloprid" = "imidacloprid",
        "X5.Metolachlor" = "5-metolachlor",
        "Tris.2.chloroetyl.phosphate" = "TCEP",
        "N.Butylbenzesulfonamide" = "N-butylbenzesulfonamide",
        "Metformin" = "metformin",
        "Terbutryn" = "terbutryn",
        "Losartan" = "losartan",
        "N.Guanylurea" = "guanylurea",
        "Bezafibrat" = "bezafibrate",
        "Tris.2.chloroisopropyl.phosphate" = "TCPP",
        "X2.Hydroxybenzothiazole" = "2-hydroxybenzothiazole",
        "Valsartansaeure" = "valsartan Acid",
        "Sulfamethoxazol" = "sulfamethoxazole",
        "Primidone" = "primidone",
        "X1.3.Diphenylguanidine" = "1,3-diphenylguanidine",
        "Venlafaxine" = "venlafaxine",
        "Metoprolol.Tartrate" = "metoprolol",
        "o.Desmethylvenlafaxine" = "O-desmethylvenlafaxine",
        "DEET" = "DEET",
        "Irbesartan" = "irbesartan",
        "Olmesartan" = "olmesartan",
        "Carbamazepine" = "carbamazepine",
        "Gabapentin.Lactam" = "gabapentin lactam",
        "Caffeine" = "caffeine",
        "X5.Methyl.benzotriazol" = "5-methylbenzotriazole",
        "Hexa.methoxymethyl.melamine" = "HMMM",
        "Hydrochlorothiazid" = "hydrochlorothiazide",
        "Benzothiazol.2.sulfonic.acid" = "BTSA",
        "Valsartan" = "valsartan",
        "Benzotriazol" = "benzotriazole",
        "Formyaminoantipyrin" = "4-formylaminoantipyrine"
      )
      
      
# --- VENN DIAGRAM --- 
      
      # Create df for venn diagram
      venn_list_df <- df %>%
        filter(Concentration > 0) %>%
        group_by(ATTRIBUTE_RiverGroup) %>%
        summarise(Compounds = list(unique(FeatureID)), .groups = "drop") %>%
        pivot_wider(names_from = ATTRIBUTE_RiverGroup, values_from = Compounds)
    
      # Set sets 
      set1 <- venn_list_df$Spree[[1]]
      set2 <- venn_list_df$Neckar[[1]]
      
      # Create the Venn diagram
      venn_data <- list(
        "Spree" = set1,
        "Neckar" = set2
      )
      
      # plot venn
      venn_targeted <- ggvenn(venn_data, show_percentage = TRUE, fill_color = c('#8795e8', '#60e4c1',"#FF8103", "#51a2d5", '#B63400',"#022f8e", '#8795e8', '#60e4c1',"#EFC000FF", "#868686FF", "#CD534CFF"),
                             stroke_size = 0.5, set_name_size = 0, stroke_color = "white",
                             auto_scale = TRUE, text_size = 7
                             
      )
      
      print(venn_targeted)
      
      # Save venn diagram
      venn_out <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Venn_TargetedCompounds.pdf"))
      
      pdf(file = venn_out, height = 3, width = 3)
      print(venn_targeted)
      dev.off()
      
      
# --- BOX PLOT ---
      
      # Create df for box plot
      summed_df <- df %>%
        group_by(Sample, ATTRIBUTE_RiverGroup) %>%
        summarise(SummedConcentration = mean(Concentration, na.rm = TRUE), .groups = "drop") %>%
        mutate(ATTRIBUTE_RiverGroup = factor(ATTRIBUTE_RiverGroup,
                                             levels = c("Spree", "Neckar")))
      
      # Wilcoxon test
      test_results <- summed_df %>%
        summarise(
          p_value = tryCatch(
            wilcox.test(SummedConcentration ~ ATTRIBUTE_RiverGroup)$p.value,
            error = function(e) NA
          )
        )
      
      print(test_results)
      
      # Boxplot
      boxplot <- ggplot(summed_df, aes(x = ATTRIBUTE_RiverGroup, y = SummedConcentration,
                                       fill = ATTRIBUTE_RiverGroup)) +
        geom_boxplot(position = position_dodge(width = 0.5)) +
        theme_bw() +
        scale_fill_manual(values = c("Neckar" = '#60e4c1', "Spree" = '#8795e8')) +
        labs(x = "River", y = "Average Scaled Peak Area", fill = "River") +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
      
      print(boxplot)
      
      # Save boxplot
      boxplot_out <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Boxplots_TargetedConcentrations_byRiver.pdf"))
      
      pdf(file = boxplot_out, height = 3, width = 5)
      print(boxplot)
      dev.off()
      
      
      
      
                                                                                  
##### 4.4 HEATMAP of FEATURE GROUPS across SITES (FIGURE 3) #####                                                                            
      
      # Set up dataframe (use scaled df without Erpe)                                                            
      long_quant <- blk_avg_fm_noErpe %>%
        tibble::rownames_to_column("Sample") %>%
        pivot_longer(-Sample, names_to = "Compound", values_to = "Concentration")
      
      # Make sure Sample column is in metadata
      md_Samples_use$Sample <- rownames(md_Samples_use)
      
      # Join metadata and quant data together
      blk_avg_md <- long_quant %>%
        mutate(FeatureID = Compound) %>%
        left_join(md_Samples_use, by = "Sample") 
      
      # Prepare df for heatmap
      corr_df <- blk_avg_md %>%
      group_by(Sample, Compound, ATTRIBUTE_RiverGroup) %>%
      mutate(Compound = factor(Compound, levels = rev(Compound_order)))
      
      # Rename compounds       
      corr_df <- corr_df %>%
      mutate(Compound = recode(Compound, !!!Compound_rename))
      print(as.vector(unique(corr_df$Compound)))      
      
      # Filter for main stem sites only
      corr_df_filtered <- corr_df %>%
      filter(FeatureID %in% Compound_order) %>%
      # filter(Sample %in% c("N01_both_1.mzML",     "N02_both_1.mzML",     "N03_both_1.mzML",     "N04_both_1.mzML",
      #              "N05_both_1.mzML",     "N06_both_1.mzML",     "N07_both_1.mzML",     "N08_both_1.mzML",     "N09_both_1.mzML",
      #              "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
      #              "S01_both_1.mzML",     "S08_both_1.mzML",     "S11_both_1.mzML",     "S13_both_1.mzML",     "S18_both_1.mzML",
      #              "S26_both_1.mzML",     "S29_both_1.mzML",     "S51_both_1.mzML",
      #              "S56_both_1.mzML",     "S57_both_1.mzML",     "S64_both_1.mzML",     "S67_both_1.mzML",
      #              "S69_both_1.mzML",     "S72_both_1.mzML",     "S73_both_1.mzML",
      #              "S74_both_1.mzML",     "S75_both_1.mzML",     "S76_both_1.mzML",     "S76a_both_1.mzML"))
      # 
      filter(Sample %in% c("N1_both_1.mzML",     "N2_both_1.mzML",     "N3_both_1.mzML",     "N4_both_1.mzML",
                           "N5_both_1.mzML",     "N6_both_1.mzML",     "N7_both_1.mzML",     "N8_both_1.mzML",     "N9_both_1.mzML",
                           "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
                           "S1_both_1.mzML",     "S2_both_1.mzML",     "S3_both_1.mzML",     "S4_both_1.mzML",     "S5_both_1.mzML",
                           "S6_both_1.mzML",     "S7_both_1.mzML",     "S8_both_1.mzML",
                           "S9_both_1.mzML",     "S10_both_1.mzML",     "S11_both_1.mzML",     "S12_both_1.mzML",
                           "S13_both_1.mzML",     "S14_both_1.mzML",     "S15_both_1.mzML",
                           "S16_both_1.mzML",     "S17_both_1.mzML",     "S18_both_1.mzML",     "S19_both_1.mzML"))
      
      # Pivot to wide format for ComplexHeatmap
      heatmap_matrix <- corr_df_filtered %>%
      group_by(Compound, Sample) %>%
      summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Sample, values_from = Concentration) %>%
      column_to_rownames("Compound") %>%
      as.matrix()
      
      library(ComplexHeatmap)
      library(circlize)
      
      # Generate heat map
      set.seed(1235)
      hmap <- Heatmap(
      heatmap_matrix,  # or use heatmap_matrix if unscaled
      heatmap_legend_param = list(title = "Scaled\nGroupPeak"),
      col = circlize::colorRamp2(c(0, 0.0001,0.5, 1), c("darkgray", "blue", "white", "red")),
      show_row_names = TRUE,
      show_column_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_dend_reorder = FALSE,
      column_dend_reorder = FALSE,
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "complete",
      width = unit(500, "mm"),
      height = unit(300, "mm"),
      row_km = 1,
      row_km_repeats = 100
      )
      
      hk <- ComplexHeatmap::draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right")
      
      print(hk)
      
      # Save plot
      output_path <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Heatmap_Scaled_TargetedCompounds_Concentrations_MainRiverSites_noErpe.pdf"))
      
      pdf(file = output_path, height = 30, width = 30)
      print(hk)
      dev.off()
      
##### 4.5 CORRELATIONS of COMPOUNDS with SITE CHARACTERISTICS (FIGURE 4) #####

# List frequently detected compounds to test correlations on      
freq_Compounds <- c("4-formylaminoantipyrine", 
                   "O-desmethylvenlafaxine", 
                   "venlafaxine", 
                   "carbamazepine",  
                   "gabapentin lactam", 
                   "primidone", 
                   "irbesartan", 
                   "olmesartan",
                   "valsartan", 
                   "valsartan Acid", 
                   "metoprolol",
                   "sulfamethoxazole",
                   "metformin", 
                   "hydrochlorothiazide",
                   "HMMM", 
                   "1,3-diphenylguanidine", 
                   "BTSA",
                   "benzotriazole",
                   "5-methylbenzotriazole",
                   "caffeine", 
                   "DEET", 
                   "terbutryn"
                   )


# Loop through metadata variables and calculate correlations
cor_summary <- map_dfr(meta_vars, function(var) {
  corr_df %>%
    filter(Compound %in% freq_Compounds) %>%
    group_by(Compound, ATTRIBUTE_RiverGroup) %>%
    summarise(
      MetadataVar = var,
      Correlation = cor(Concentration, .data[[var]], use = "complete.obs", method = "spearman"),
      Pvalue = cor.test(Concentration, .data[[var]], method = "spearman")$p.value,
      .groups = "drop"
    )
})


# Order MetadataVar factor levels for plotting
cor_summary <- cor_summary %>%
  mutate(MetadataVar = factor(MetadataVar, levels = meta_vars),
         Compound = factor(Compound, levels = rev(freq_Compounds)))


cor_summary <- cor_summary %>%
  mutate(Signif = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))



# Define output file path
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Correlation_TargetedCompounds_SpearmanRho_AllSites_noErpe.csv"))

# Write the CSV file
write.csv(cor_summary, output_path, row.names = TRUE)

# Plot correlations in Neckar
neckar <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Neckar" & !is.na(cor_summary$Compound),], aes(x = MetadataVar, y = Compound)) +
  geom_point(aes(size = abs(Correlation), color = Correlation)) +
  geom_text(aes(label = Signif), size = 5, vjust = 0.8) + #vjust = -1.2
  scale_color_gradient2(
    low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
    name = "spearman r"
  ) +
  scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
  theme_bw() +
  labs(
    # title = "Correlation Between Group PA and Metadata",
    x = "Metadata Variable",
    y = "Compound"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

# Plot correlations in Spree
spree <- ggplot(cor_summary[cor_summary$ATTRIBUTE_RiverGroup == "Spree" & !is.na(cor_summary$Compound),], aes(x = MetadataVar, y = Compound)) +
  geom_point(aes(size = abs(Correlation), color = Correlation)) +
  geom_text(aes(label = Signif), size = 5, color = "black", vjust = 0.8) +
  scale_color_gradient2(
    low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
    name = "spearman r"
  ) +
  scale_size(range = c(0.5, 12), limits = c(0, 1), name = "|p|") +
  theme_bw() +
  labs(
    # title = "Correlation Between Group PA and Metadata",
    x = "Metadata Variable",
    y = "Compound"
  ) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    legend.position = "none"
  )


# Save correlation results figures together
output_path <- file.path(Directory, "03_Output_Figures", paste0(Sys.Date(), "_Correlation_TargetedCompounds_SpearmanRho_AllSites_noErpe.pdf"))

pdf(file = output_path, height = 10.8, width = 12.7)
plot_grid( spree,  neckar, ncol = 2, nrow = 1, rel_widths = c(1,1), rel_heights = c(1,1,1,1), label_size = 24, label_y = 0.9,
           align = "v")
dev.off()



                                                                                
                                                                                  
                                                                                  
                                                                                  
