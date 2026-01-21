### Code to assess and compare organic compounds' targeted concentrations and non-targeted peak areas along and between the Neckar and Spree Rivers in Germany ###
### Code replicates main Hitchhikers script for calling files, blank removal and data set selection for data analysis, which is linked below
### #https://github.com/Functional-Metabolomics-Lab/FBMN-STATS/blob/main/R/.ipynb_checkpoints/Stats_Untargeted_Metabolomics-checkpoint.ipynb

# The main steps are to (1) load the libraries and (2) call the data files to be processed. 
# Subsequently (3) the non-targeted and (4) targeted data area analyzed, which entails producing
# an overview of data in the two rivers, how peak areas or concentrations vary along the rivers,  
# and how peak areas or concentrations correlate with the rivers' site characteristics. 

##### 1. LOAD LIBRARIES #####

library("tidyverse")
library("tidyr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("igraph")
library("ggvenn")
library("patchwork")
library("ComplexHeatmap")
library("rlang")


##### 2. CALL FILES and SET UP #####

#set directory
Directory <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Neckar and Spree Rivers Organic Pollutants Data and Code"
setwd(Directory)

### Call files ###
# Create folders for output if not already there
dir.create(file.path(Directory, "02_Output_DataProducts"))
dir.create(file.path(Directory, "03_Output_Figures"))

# Call file with feature annotation groups
fileName_anthro <- file.path(Directory, "01_Data", "FeatureGroups_Assignment.csv")
df_csv_anthro <- read.csv(fileName_anthro,  head=TRUE, sep=",")
colnames(df_csv_anthro)

# Call file with molecular network results
fileName_edges <- file.path(Directory, "01_Data", "FeatureEdges_MolecularNetwork.csv")
df_csv_edges <- read.csv(fileName_edges,  head=TRUE, sep=",")
dim(df_csv_edges)

# Call file with targeted concentrations
fileName_TA <- file.path(Directory, "01_Data", "Targeted_Concentrations.csv")
df_csv_TAorig <- read.csv(fileName_TA,  head=TRUE, sep=",")
dim(df_csv_TAorig)


### Call files that need reformatting (NTA files) ###
# Load files from directory
file_names <- list.files('.') #list all the files in the working directory (mentioned by 'dot symbol')
file_names

# Identify names of feature quant and metadata tables in directory
input_str <- "NonTargeted_PeakAreas_Raw.csv,Sample_Metadata.txt,FeatureAnnotations_GNPSoutput.tsv"

# Read the files indicated by the indecies 
input <- (strsplit(input_str, ",")[[1]])

# Read feature quant table
first_line <- readLines(file.path(Directory, "01_Data", input[1]), n = 1)
if (length(strsplit(first_line, ';')[[1]]) > 1) {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F, sep = ';') # in case, ';' is the separator
} else {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F)
}

# Read metadata
md <- read.csv(file.path(Directory, "01_Data", input[2]), header = T, check.names = F, sep = '\t') # mention seperator as "/t"(tab-separated) in case of txt or tsv fi

### Annotations ###
# Load GNPS annotations:
an_gnps <- read.csv(file.path(Directory, "01_Data", input[3]), header = T, check.names = F, sep = '\t') 

# Check table dimensions
dim(ft) 
dim(md)
dim(an_gnps) 


### Reformat files ###
# Function to summarize metadata
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

# Check out metadata dimensions & content
ncol(md) #number of columnns of metadata
InsideLevels(md[, 2:ncol(md)]) #excluding 1st filename

# Check compatibility of annotations with feature table ids by checking their class
identical(class(ft$`row ID`),class(an_gnps$`#Scan#`))

# Arrange annotations by Scan number
an_final <- an_gnps %>% arrange(`#Scan#`) #arranging by scan ID

# Reformat GNPS annotations
# Function to compare between gnps_compound_name and its library match (analog)
combine_names <- function(compound_name) {
  return(paste(compound_name))
}

# Consolidate multiple annotations for a single '#Scan#' into one combined name
an_final_single <- an_final %>%
  group_by(`#Scan#`) %>%
  summarise(Combined_Name = combine_names(Compound_Name[1])) %>%
  ungroup() %>%
  as.data.frame()

# Merge annotations with feature table
ft_an <- merge(ft, an_final_single, by.x="row ID", by.y="#Scan#", all.x= TRUE) 
head(ft_an, 2)
dim(ft_an) 


# Arrange feature table and metadata in same order and create duplicate (working) files
new_ft <- ft 
new_md <- md

# Clean the new files
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft)) #removing Peak area extensions from the column names of ft
new_ft <- new_ft[order(new_ft$`row ID`),,drop=F] #arranging the rows of ft file by  by ascending order of row ID
new_ft <- new_ft[,colSums(is.na(new_ft))<nrow(new_ft)] #removing if any NA columns present in the ft file,
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #removing if any NA columns present in the md file,
new_md <- new_md[apply(new_md != "", 1, any), ] # Removing rows that are completely filled with empty strings,
new_md <- new_md[, apply(new_md != "", 2, any)] # Removing columns that are completely filled with empty strings

# Remove the (front & tail) spaces, if any present, from the filenames of md
new_md$filename <- trimws(new_md$filename, which = c("both"))
rownames(new_md) <- new_md$filename
new_md <- new_md[, -which(names(new_md) == "filename")]

# Update row names of feature table
if(exists("ft_an")){identical(ft_an$`row ID`,new_ft$`row ID`)} #should return TRUE if you have annotation file

# Change the row names of the files into the combined name as "XID_mz_RT":
rownames(new_ft) <- paste(paste0("X",new_ft$`row ID`),
                          round(new_ft$`row m/z`,digits = 3),
                          round(new_ft$`row retention time`,digits = 3),
                          if(exists("ft_an")){ft_an$Combined_Name}, 
                          sep = '_') 

# Remove the trailing underscore at rownames
rownames(new_ft) <- sub("_$", "", rownames(new_ft)) 


# In the feature table, identify which columns correspond to samples
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

# Check overlap between the feature table and metadata
new_ft <- new_ft[,order(colnames(new_ft)), drop=F] #ordering the ft by its column names
new_md <- new_md[order(rownames(new_md)),, drop=F] #ordering the md by the 1st column filename

# How many files in the metadata are also present in the feature table
table(rownames(new_md) %in% colnames(new_ft))

# Are the sample names the same
identical(rownames(new_md), colnames(new_ft))

# Which file names in the metadata are not in the feature table?
setdiff(rownames(new_md),colnames(new_ft))
# print(colnames(new_ft)) # uncomment to check the column names of new_ft

# Checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

# Remove rows with only 0 peak area values
new_ft <- new_ft[rowSums(new_ft != 0) > 0,]
dim(new_ft)

# Print summary of the feature table
cat("The number of rows and non-zero columns in our new ft is:",dim(new_ft))

# Transpose and merge feature table and metadata
ft_t <- as.data.frame(t(new_ft)) #transposing the ft
ft_t <- ft_t %>% mutate_all(as.numeric)  #converting all values to numeric
identical(rownames(new_md),rownames(ft_t)) #should return TRUE now


##### 3 PROCESSING For NON-TARGETED FEATURES #####
##### 3.1 BLANK REMOVAL and SUBTRACTION #####

# Get the index levels in your metadata
InsideLevels(new_md)

# Enter indecies
# Sample_attribute <- as.numeric(readline('Enter the index number of the attribute containing sample and blanks information: '))
sample_attribute <- as.numeric("7")
unique_sampletypes <- unique(new_md[, sample_attribute])

# Display the unique sample types along with their index & select desired indecies for blanks and samples
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

# Get the corresponding rows from ft_t
Blank <- ft_t[which(rownames(ft_t) %in% (rownames(md_Blank))), , drop=F]
Samples <- ft_t[which(rownames(ft_t) %in% (rownames(md_Samples))), , drop=F]
InsideLevels(new_md)

# Check out results
dim(Blank) 
dim(Samples)


### Blank Removal Step ###
# Set the Blank removal cutoff value
# When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
# Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
Cutoff <- as.numeric(0.1)

# Get mean for every feature in Blank and Samples in a data frame named 'Avg_ft'
Avg_ft <- data.frame(Avg_blank=colMeans(Blank, na.rm= F)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ft$Avg_samples <- colMeans(Samples, na.rm= F) # adding another column 'Avg_samples' for feature means of samples

# Get the ratio of Blank vs Sample
Avg_ft$Ratio_blank_Sample <- (Avg_ft$Avg_blank+1)/(Avg_ft$Avg_samples+1)

# Create a bin with 1s when the ratio > Cutoff, else put 0s
Avg_ft$Bg_bin <- ifelse(Avg_ft$Ratio_blank_Sample > Cutoff, 1, 0 )

# Calculate the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ft)))
print(paste("No.of Background or noise features:",sum(Avg_ft$`Bg_bin` ==1,na.rm = T)))
print(paste("No.of features after excluding noise:",(ncol(Samples) - sum(Avg_ft$`Bg_bin` ==1,na.rm = T))))

# Create a new, blank removed dataframe
blk_rem_1 <- merge(as.data.frame(t(Samples)), Avg_ft, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) %>% #removing the last 4 columns
  column_to_rownames(var="Row.names") 

# Set blank removed dataframe as data frame and transpose
  blk_rem <- as.data.frame(t(blk_rem_1))
  
  
### Blank Subtraction step ###
# Use max value observed in the blanks 
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
# Define output file path & write the .csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_BlanksRemoved_and_MaxSubtracted_with_cutoff_',Cutoff,'.csv'))
write.csv(t(blk_rem_sub), output_path, row.names = TRUE)
        
# Ensure the final dataframe is formatted correctly
blk_rem_sub <- as.data.frame(blk_rem_sub)

# Review dimensions of feature table without blanks
dim(blk_rem_sub)

# Review dimensions of metadata without the blanks info 
dim(md_Samples)


##### 3.2 PEAK AREA REPLICATE AVERAGING and SCALING #####

# Make copies of data tables
blk_rem_use <- blk_rem_sub
md_Samples_use <- md_Samples


### Average the replicate samples ###
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

# Remove zero-variance columns
blk_avg_pa <- blk_avg_pa[, apply(blk_avg_pa, 2, sd) != 0]
dim(blk_avg_pa)

# Define output file path & write .csv
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged.csv'))
write.csv(t(blk_avg_pa), output_path, row.names = TRUE)


### Scale the averaged feature peak areas - with Erpe sample ###
# Make function to normalize features by the maximum value observed for each across all samples
normalize <- function(x) {
  return(x / max(x))
}

# Apply the normalization function to each column
blk_avg_fm <- as.data.frame(apply(blk_avg_pa, 2, normalize))

# Define output file path & write .csv
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged_Scaled.csv'))
write.csv(t(blk_avg_fm), output_path, row.names = TRUE)


### Scale the averaged feature peak areas - without Erpe sample ###
# Apply normalization to all columns except "Erpe_both_1.mzML" and drop that column
blk_avg_fm_noErpe <- as.data.frame(
  apply(blk_avg_pa[rownames(blk_avg_pa) != c("Erpe_both_1.mzML"),], 2, normalize))

# Define output file path & write .csv
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_PeakAreas_ReplicatesAveraged_Scaled_noErpe.csv'))
write.csv(t(blk_avg_fm_noErpe), output_path, row.names = TRUE)


##### 3.3 SUBSET ANTHROPOGENIC FEATURES USING ANNOTATIONS and MOLECULAR NETWORKS#####

# Create a copy of the original feature annotations group table
features_use <- df_csv_anthro

# Define the group columns
group_cols <- paste0("Group", 0:2)

# Find the rows where cosine is less than 0.8 and the MZerror is greater than 10 - the annotations will not be used from these features
rows_to_blank <- which(!is.na(features_use$Cosine) & features_use$Cosine < 0.8 & features_use$IonMode != "Negative" | !is.na(features_use$Mzerror) & features_use$Mzerror >= 10 & features_use$IonMode != "Negative")
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

# Add these to the dataframe
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

# Define output file path & write .csv
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Feature_Annotations_MNneighbors_and_Groups.csv"))
write.csv(all_features_df_use, output_path, row.names = TRUE)

                                                          
##### 3.4 DATA OVERVIEW and RIVER COMPARISON (FIGURE 2) #####                 

# 
# 1) PREPARE LONG TABLE WITH METADATA  
# 

long_quant <- blk_avg_fm_noErpe %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Feature", values_to = "PeakArea")

md_Samples_use$Sample <- rownames(md_Samples_use)

blk_avg_md <- long_quant %>%
  mutate(FeatureID = as.numeric(gsub(".*X([0-9]+).*", "\\1", Feature))) %>%
  left_join(
    all_features_df_use %>% select(FeatureID, Source, Group0, Group1, Group2),
    by = "FeatureID"
  ) %>%
  left_join(md_Samples_use, by = "Sample")

blk_avg_md$Group0[is.na(blk_avg_md$Group0)] <- "unknown"
blk_avg_md$Group1[is.na(blk_avg_md$Group1)] <- "unknown"
blk_avg_md$Group2[is.na(blk_avg_md$Group2)] <- "unknown"

group1_colors <- c(
  "biological compound" = "#1B9E77",
  "other consumer industrial" = "#7570B3",
  "pesticide"  = "#66A61E",
  "pharmaceutical" = "#377EB8",
  "polymer-related" = "#FF7F00",
  "unknown" = "#999999"
)

# 
# 2) DETECTION TABLE + SETS AT THRESHOLDS
#    (>0 implemented as "gt0" = n_detect >= 1)
# 

build_det_tbl <- function(df) {
  det_tbl <- df %>%
    distinct(Sample, ATTRIBUTE_RiverGroup, FeatureID, PeakArea) %>%
    mutate(Detected = PeakArea > 0) %>%
    group_by(ATTRIBUTE_RiverGroup, FeatureID) %>%
    summarise(n_detect = sum(Detected, na.rm = TRUE), .groups = "drop")
  
  n_sites <- df %>%
    distinct(Sample, ATTRIBUTE_RiverGroup) %>%
    count(ATTRIBUTE_RiverGroup, name = "n_sites")
  
  det_tbl %>%
    left_join(n_sites, by = "ATTRIBUTE_RiverGroup") %>%
    mutate(prop_detect = n_detect / n_sites)
}

make_sets_from_det <- function(det_tbl, threshold) {
  
  get_set <- function(river) {
    x <- det_tbl %>% filter(ATTRIBUTE_RiverGroup == river)
    
    if (is.character(threshold) && threshold == "gt0") {
      x %>% filter(n_detect >= 1) %>% pull(FeatureID) %>% unique()
    } else {
      x %>% filter(prop_detect >= threshold) %>% pull(FeatureID) %>% unique()
    }
  }
  
  list(
    set_spree  = get_set("Spree"),
    set_neckar = get_set("Neckar")
  )
}

make_threshold_sets <- function(df, threshold, type = c("count", "prop")) {
  
  type <- match.arg(type)
  det_tbl <- build_det_tbl(df)
  
  if (type == "count") {
    # threshold is a count (e.g., 10)
    get_set <- function(river) {
      det_tbl %>%
        filter(ATTRIBUTE_RiverGroup == river, n_detect >= threshold) %>%
        pull(FeatureID) %>% unique()
    }
    set_spree  <- get_set("Spree")
    set_neckar <- get_set("Neckar")
    label_suffix <- paste0("minSites_", threshold, "_byRiver")
  } else {
    # threshold is a proportion OR "gt0"
    sets <- make_sets_from_det(det_tbl, threshold = threshold)
    set_spree  <- sets$set_spree
    set_neckar <- sets$set_neckar
    
    label_suffix <- if (is.character(threshold) && threshold == "gt0") {
      "minProp_gt0_byRiver"
    } else {
      paste0("minProp_", as.integer(round(threshold * 100)), "_byRiver")
    }
  }
  
  list(
    set_spree  = set_spree,
    set_neckar = set_neckar,
    union_ids  = union(set_spree, set_neckar),
    label_suffix = label_suffix
  )
}

# 
# 3) PLOTTING FUNCTION (PIE + VENN + BOX + COUNTS)
#    


plot_pie_venn_box <- function(data, source_filter, label, Directory, venn_sets,
                              save_plots = TRUE) {
  
  df <- data %>% filter(!!enquo(source_filter))
  
  # 
  # FEATURE COUNTS PER RIVER
  # 
  feature_counts <- df %>%
    filter(PeakArea > 0) %>%
    distinct(FeatureID, ATTRIBUTE_RiverGroup) %>%
    count(ATTRIBUTE_RiverGroup, name = "n_features") %>%
    complete(
      ATTRIBUTE_RiverGroup = c("Spree", "Neckar"),
      fill = list(n_features = 0)
    )
  
  message("\n==============================")
  message("Plot: ", label)
  message("Detected unique features per river:")
  print(feature_counts)
  
  # 
  # PIE CHART
  # 
  pie_df <- df %>%
    distinct(FeatureID, Group1) %>%
    group_by(Group1) %>%
    summarise(FeatureCount = n(), .groups = "drop") %>%
    mutate(
      Percent = FeatureCount / sum(FeatureCount) * 100,
      Label = paste0(FeatureCount, "\n(", round(Percent, 1), "%)")
    )
  
  pie_plot <- ggplot(pie_df, aes(x = "", y = FeatureCount, fill = Group1)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.7), size = 4) +
    scale_fill_manual(values = group1_colors) +
    theme_void()
  
  if (save_plots) {
    ggsave(
      filename = file.path(Directory, "03_Output_Figures",
                           paste0(Sys.Date(), "_PieChart_", label, "_Group1.pdf")),
      plot = pie_plot, height = 4, width = 4
    )
  }
  
  # 
  # VENN DIAGRAM (RIVER SETS)
  # 
  set_spree  <- venn_sets$set_spree
  set_neckar <- venn_sets$set_neckar
  
  venn_sizes <- data.frame(
    n_spree = length(set_spree),
    n_neckar = length(set_neckar),
    n_overlap = length(intersect(set_spree, set_neckar)),
    n_union = length(union(set_spree, set_neckar))
  )
  
  message("Venn set sizes:")
  print(venn_sizes)
  
  venn_plot <- ggvenn(
    list(Spree = set_spree, Neckar = set_neckar),
    show_percentage = TRUE,
    fill_color = c("#8795e8", "#60e4c1"),
    stroke_color = "white",
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 4
  )
  
  if (save_plots) {
    ggsave(
      filename = file.path(Directory, "03_Output_Figures",
                           paste0(Sys.Date(), "_Venn_", label, ".pdf")),
      plot = venn_plot, height = 4, width = 4
    )
  }
  
  # 
  # BOX PLOT
  # 
  summed_df <- df %>%
    group_by(Sample, ATTRIBUTE_RiverGroup) %>%
    summarise(SummedPeakArea = mean(PeakArea, na.rm = TRUE), .groups = "drop") %>%
    mutate(ATTRIBUTE_RiverGroup = factor(ATTRIBUTE_RiverGroup,
                                         levels = c("Spree", "Neckar")))
  
  box_plot <- ggplot(
    summed_df,
    aes(x = ATTRIBUTE_RiverGroup, y = SummedPeakArea, fill = ATTRIBUTE_RiverGroup)
  ) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(Spree = "#8795e8", Neckar = "#60e4c1")) +
    labs(x = "River", y = "Average scaled peak area")
  
  if (save_plots) {
    ggsave(
      filename = file.path(Directory, "03_Output_Figures",
                           paste0(Sys.Date(), "_Boxplots_", label, "_PeakAreas_byRiver.pdf")),
      plot = box_plot, height = 3, width = 5
    )
  }
  
  invisible(list(
    pie_plot = pie_plot,
    venn_plot = venn_plot,
    box_plot = box_plot,
    feature_counts = feature_counts,
    venn_sizes = venn_sizes
  ))
}

# 
# 4) RUN ALL THRESHOLDS (PIE+VENN+BOX) 
# 

run_plots_for_thresholds <- function(df, Directory,
                                     count_thresholds = c(),
                                     prop_thresholds  = list("gt0", 0.50)) {
  
  plan_specs <- list(
    count = count_thresholds,
    prop  = prop_thresholds
  )
  
  results <- list()
  
  for (type_name in names(plan_specs)) {
    for (th in plan_specs[[type_name]]) {
      
      # Split FIRST (so thresholds computed per subset)
      df_anthro <- df %>% filter(Source %in% c("all anthropogenic", "MN"))
      df_other  <- df %>% filter(!Source %in% c("all anthropogenic", "MN"))
      
      plan_a <- make_threshold_sets(df_anthro, threshold = th, type = type_name)
      plan_o <- make_threshold_sets(df_other,  threshold = th, type = type_name)
      
      df_use_a <- df_anthro %>% filter(FeatureID %in% plan_a$union_ids)
      df_use_o <- df_other  %>% filter(FeatureID %in% plan_o$union_ids)
      
      suffix <- plan_a$label_suffix  # same naming scheme
      
      # Anthro plots
      if (nrow(df_use_a) > 0) {
        plot_pie_venn_box(
          data = df_use_a,
          source_filter = TRUE,
          label = paste0("AnthropogenicFeatures_", suffix),
          Directory = Directory,
          venn_sets = plan_a
        )
      }
      
      # Other plots
      if (nrow(df_use_o) > 0) {
        plot_pie_venn_box(
          data = df_use_o,
          source_filter = TRUE,
          label = paste0("OtherFeatures_", suffix),
          Directory = Directory,
          venn_sets = plan_o
        )
      }
    }
  }
  
  invisible(results)
}

# 
# 5) PURE OVERLAP GRAPH (SCATTER + LINE)
# 

overlap_curve <- function(df_subset,
                          group_label,
                          thresholds = list("gt0", 0.25, 0.50, 0.75, 1.00)) {
  
  det_tbl <- build_det_tbl(df_subset)
  
  bind_rows(lapply(thresholds, function(th) {
    
    sets <- make_sets_from_det(det_tbl, threshold = th)
    
    s <- sets$set_spree
    n <- sets$set_neckar
    inter <- length(intersect(s, n))
    uni   <- length(union(s, n))
    
    x_label <- if (is.character(th) && th == "gt0") ">0" else as.character(th * 100)
    x_value <- if (is.character(th) && th == "gt0") 0 else th * 100
    
    tibble(
      FeatureGroup = group_label,
      Threshold = x_label,
      ThresholdPct = x_value,
      n_spree = length(s),
      n_neckar = length(n),
      n_overlap = inter,
      n_union = uni,
      overlap_pct_jaccard = ifelse(uni == 0, NA_real_, 100 * inter / uni)
    )
  }))
}

# helper: overlap plot object 
make_overlap_plot <- function(overlap_df, feature_group) {
  ggplot(overlap_df[overlap_df$FeatureGroup == feature_group,],
         aes(x = ThresholdPct, y = overlap_pct_jaccard, group = FeatureGroup)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = c(">0", "25", "50", "75", "100")
    ) +
    labs(
      x = "Det. freq. threshold \nin each river (%)",
      y = "Overlap (%)"
    ) +
    theme_bw(base_size = 14)
}

# 
# 6) COMPILED 2x4 PLOTS FOR A SINGLE THRESHOLD
# 
# 1) split once
df_anthro <- blk_avg_md %>% filter(Source %in% c("all anthropogenic", "MN"))
df_other  <- blk_avg_md %>% filter(!Source %in% c("all anthropogenic", "MN"))

# 2) build overlap ONCE (both groups in one table)
thresholds_use <- list("gt0", 0.25, 0.50, 0.75, 1.00)

overlap_df <- bind_rows(
  overlap_curve(df_anthro, "Anthropogenic", thresholds = thresholds_use),
  overlap_curve(df_other,  "Other",          thresholds = thresholds_use)
)

print(overlap_df)

# 3) save overlap table (create folder if needed)
out_dir <- file.path(Directory, "02_Output_DataProducts")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

csv_fp <- file.path(out_dir, paste0(Sys.Date(), "_OverlapCurve_NTAfeatures_Table.csv"))
write.csv(overlap_df, csv_fp, row.names = FALSE)

# 4) overlap plot objects reused in BOTH compiled PDFs
overlap_plot_a <- make_overlap_plot(overlap_df, "Anthropogenic")
overlap_plot_o <- make_overlap_plot(overlap_df, "Other")

compile_detection_rate_plot <- function(df, Directory, threshold,
                                        overlap_df = NULL,
                                        overlap_plot_a = NULL,
                                        overlap_plot_o = NULL) {
  
  if (!"package:patchwork" %in% search()) {
    library(patchwork)
  }
  
  # Split FIRST (so thresholds computed per subset)
  df_anthro <- df %>% filter(Source %in% c("all anthropogenic", "MN"))
  df_other  <- df %>% filter(!Source %in% c("all anthropogenic", "MN"))
  
  # sets for this threshold (per subset)
  plan_a <- make_threshold_sets(df_anthro, threshold = threshold, type = "prop")
  plan_o <- make_threshold_sets(df_other,  threshold = threshold, type = "prop")
  
  df_use_a <- df_anthro %>% filter(FeatureID %in% plan_a$union_ids)
  df_use_o <- df_other  %>% filter(FeatureID %in% plan_o$union_ids)
  
  suffix <- plan_a$label_suffix 
  
  # if overlap not supplied, build it
  if (is.null(overlap_df) || is.null(overlap_plot_a) || is.null(overlap_plot_o)) {
    thresholds_use <- list("gt0", 0.25, 0.50, 0.75, 1.00)
    overlap_df <- bind_rows(
      overlap_curve(df_anthro, "Anthropogenic", thresholds = thresholds_use),
      overlap_curve(df_other,  "Other",          thresholds = thresholds_use)
    )
    overlap_plot_a <- make_overlap_plot(overlap_df, "Anthropogenic")
    overlap_plot_o <- make_overlap_plot(overlap_df, "Other")
  }
  
  # Pie/Venn/Box (do NOT save individually during compilation)
  if (nrow(df_use_a) > 0) {
    out_a <- plot_pie_venn_box(
      data = df_use_a,
      source_filter = TRUE,
      label = paste0("AnthropogenicFeatures_", suffix),
      Directory = Directory,
      venn_sets = plan_a,
      save_plots = FALSE
    )
  } else {
    out_a <- list(
      pie_plot  = ggplot() + theme_void() + labs(title = "No anthropogenic features at this threshold"),
      venn_plot = ggplot() + theme_void(),
      box_plot  = ggplot() + theme_void()
    )
  }
  
  if (nrow(df_use_o) > 0) {
    out_o <- plot_pie_venn_box(
      data = df_use_o,
      source_filter = TRUE,
      label = paste0("OtherFeatures_", suffix),
      Directory = Directory,
      venn_sets = plan_o,
      save_plots = FALSE
    )
  } else {
    out_o <- list(
      pie_plot  = ggplot() + theme_void() + labs(title = "No other features at this threshold"),
      venn_plot = ggplot() + theme_void(),
      box_plot  = ggplot() + theme_void()
    )
  }
  
  # Compiled 2x4 layout: [pie | venn | overlap | box] x 2 rows
  compiled_plot <- (out_a$pie_plot | out_a$venn_plot | overlap_plot_a | out_a$box_plot) /
    (out_o$pie_plot | out_o$venn_plot | overlap_plot_o | out_o$box_plot)
  
  out_fp <- file.path(
    Directory, "03_Output_Figures",
    paste0(Sys.Date(), "_CompiledPlots_", suffix, "_AnthroOther_2x4.pdf")
  )
  
  ggsave(out_fp, compiled_plot, height = 8, width = 16)
  
  invisible(list(
    compiled_plot = compiled_plot,
    overlap_df = overlap_df,
    out_fp = out_fp,
    suffix = suffix
  ))
}

# 5) compile for gt0 + 50% using the SAME overlap plots
compiled_gt0 <- compile_detection_rate_plot(
  df = blk_avg_md,
  Directory = Directory,
  threshold = "gt0",
  overlap_df = overlap_df,
  overlap_plot_a = overlap_plot_a,
  overlap_plot_o = overlap_plot_o
)

compiled_50 <- compile_detection_rate_plot(
  df = blk_avg_md,
  Directory = Directory,
  threshold = 0.50,
  overlap_df = overlap_df,
  overlap_plot_a = overlap_plot_a,
  overlap_plot_o = overlap_plot_o
)





##### 3.5 HEATMAP of FEATURE GROUPS across SITES (FIGURE 3) #####

# Prepare the data by grouping features 
        # Reformat scaled feature quant table from wide to long format - Use all samples except Erpe
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

        # Calculate Group1 sums per sample
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
        filter(Sample %in% c("N1_both_1.mzML",     "N2_both_1.mzML",     "N3_both_1.mzML",     "N4_both_1.mzML",
                             "N5_both_1.mzML",     "N6_both_1.mzML",     "N7_both_1.mzML",     "N8_both_1.mzML",     "N9_both_1.mzML",
                             "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
                             "S1_both_1.mzML",     "S2_both_1.mzML",     "S3_both_1.mzML",     "S4_both_1.mzML",     "S5_both_1.mzML",
                             "S6_both_1.mzML",     "S7_both_1.mzML",     "S8_both_1.mzML",
                             "S9_both_1.mzML",     "S10_both_1.mzML",     "S11_both_1.mzML",     "S12_both_1.mzML",
                             "S13_both_1.mzML",     "S14_both_1.mzML",     "S15_both_1.mzML",
                             "S16_both_1.mzML",     "S17_both_1.mzML",     "S18_both_1.mzML",     "S19_both_1.mzML"))
      # Check number of Groups
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

      # Generate heatmap dataframe
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

      # Plot heatmap
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

      # Add p-value to summary
      cor_summary <- cor_summary %>%
        mutate(Signif = case_when(
          Pvalue < 0.001 ~ "***",
          Pvalue < 0.01 ~ "**",
          Pvalue < 0.05 ~ "*",
          TRUE ~ ""
        ))

      # Define output file path & write .csv
      output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Correlation_GroupedFeatsVsSiteChar_SpearmanRho_AllSites_noErpe.csv"))
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


                                                                                  
                                                                                  
##### 4 PROCESSING For TARGETED COMPOUNDS #####
##### 4.1 BLANK REMOVAL and SUBTRACTION #####
      
    # Data used for targeted data analysis has already been subjected
    # to blank subtraction in a previous, manual process. 
        
##### 4.2 CONCENTRATION SCALING #####
        
        # Make copies of data tables
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
        
        ### Scale the averaged feature peak areas - with Erpe sample ###
        ## Make function to normalize features by the maximum value observed for each across all samples
        normalize <- function(x) {
          return(x / max(x))
        }
        
        # Applying the normalization function to each column
        blk_avg_fm <- as.data.frame(apply(blk_rem_new, 2, normalize))
        
        # Define output file path & write .csv
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_Concentrations_Scaled.csv'))
        write.csv(blk_avg_fm, output_path, row.names = TRUE)
        
        ### Scale the averaged feature peak areas - without Erpe sample ###
        # Apply normalization to all columns except "Erpe_both_1.mzML" and drop that column
        blk_avg_fm_noErpe <- as.data.frame(apply(blk_rem_new_noErpe, 2, normalize))
        
        # Define output file path & write .csv
        output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(),'_Concentrations_Scaled_noErpe.csv'))
        write.csv(blk_avg_fm_noErpe, output_path, row.names = TRUE)
        
##### 4.3 DATA OVERVIEW and RIVER COMPARISON (FIGURE 2) #####
   
      # Establish prefered order of compounds
      Compound_order <- rev(c("Formyaminoantipyrin",
                             "o.Desmethylvenlafaxine", "Venlafaxine",
                             "Carbamazepine",  "Gabapentin.Lactam", "Primidone",
                             "Irbesartan","Losartan", "Olmesartan","Valsartan", "Valsartansaeure",  "Metoprolol.Tartrate", "Metoprolol.acid",
                             "Sulfamethoxazol",
                             "Bezafibrat", "Hydrochlorothiazid", "Metformin", "N.Guanylurea",
                             "Triphenylphosphate", "Tris.2.chloroetyl.phosphate", "Tris.2.chloroisopropyl.phosphate",
                             "Hexa.methoxymethyl.melamine",
                             "X1.3.Diphenylguanidine",
                             "Benzothiazol.2.sulfonic.acid", "X2.Hydroxybenzothiazole",
                             "Benzotriazol","X5.Methyl.benzotriazol",
                             "N.Butylbenzesulfonamide",
                             "Caffeine", "Tributylamine",
                             "DEET", "Flufenacet", "X5.Metolachlor",  "Terbutryn",
                             "Imidacloprid" ))

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
        "Valsartansaeure" = "valsartan acid",
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

      #
      # 1) PREP + CLEAN 
      # 
      
      long_quant_targeted <- blk_avg_fm %>%
        tibble::rownames_to_column("Sample") %>%
        pivot_longer(-Sample, names_to = "Compound", values_to = "Concentration")
      
      md_Samples_use$Sample <- rownames(md_Samples_use)
      
      df_targeted <- long_quant_targeted %>%
        mutate(FeatureID = Compound) %>%            
        left_join(md_Samples_use, by = "Sample")
      
      # Optional: rename compounds (robust; avoids recode/!!! issues)
      # Compound_rename must be a named character vector: names = old, values = new
      df_targeted <- df_targeted %>%
        mutate(
          Compound_clean = dplyr::if_else(
            Compound %in% names(Compound_rename),
            unname(Compound_rename[Compound]),
            Compound
          ),
          FeatureID_clean = Compound_clean
        )
      
      # 
      # 2) DETECTION TABLE + SETS AT THRESHOLDS
      # 
      
      build_det_tbl_targeted <- function(df) {
        det_tbl <- df %>%
          distinct(Sample, ATTRIBUTE_RiverGroup, FeatureID_clean, Concentration) %>%
          mutate(Detected = Concentration > 0) %>%
          group_by(ATTRIBUTE_RiverGroup, FeatureID_clean) %>%
          summarise(n_detect = sum(Detected, na.rm = TRUE), .groups = "drop")
        
        n_sites <- df %>%
          distinct(Sample, ATTRIBUTE_RiverGroup) %>%
          count(ATTRIBUTE_RiverGroup, name = "n_sites")
        
        det_tbl %>%
          left_join(n_sites, by = "ATTRIBUTE_RiverGroup") %>%
          mutate(prop_detect = n_detect / n_sites)
      }
      
      make_sets_from_det_targeted <- function(det_tbl, threshold) {
        get_set <- function(river) {
          x <- det_tbl %>% filter(ATTRIBUTE_RiverGroup == river)
          
          if (is.character(threshold) && threshold == "gt0") {
            x %>% filter(n_detect >= 1) %>% pull(FeatureID_clean) %>% unique()
          } else {
            x %>% filter(prop_detect >= threshold) %>% pull(FeatureID_clean) %>% unique()
          }
        }
        
        list(
          set_spree  = get_set("Spree"),
          set_neckar = get_set("Neckar")
        )
      }
      
      # 
      # 3) OVERLAP CURVE (table) + overlap plot object 
      # 
      
      make_overlap_curve_targeted <- function(det_tbl,
                                              thresholds = list("gt0", 0.25, 0.50, 0.75, 1.00)) {
        
        bind_rows(lapply(thresholds, function(th) {
          sets <- make_sets_from_det_targeted(det_tbl, threshold = th)
          s <- sets$set_spree
          n <- sets$set_neckar
          inter <- length(intersect(s, n))
          uni <- length(union(s, n))
          
          x_label <- if (is.character(th) && th == "gt0") ">0" else as.character(th * 100)
          x_value <- if (is.character(th) && th == "gt0") 0 else th * 100
          
          tibble(
            Threshold = x_label,
            ThresholdPct = x_value,
            n_spree = length(s),
            n_neckar = length(n),
            n_overlap = inter,
            n_union = uni,
            overlap_pct_jaccard = ifelse(uni == 0, NA_real_, 100 * inter / uni)
          )
        }))
      }
      
      make_overlap_plot_targeted <- function(overlap_df) {
        ggplot(overlap_df, aes(x = ThresholdPct, y = overlap_pct_jaccard)) +
          geom_line() +
          geom_point(size = 2) +
          scale_x_continuous(
            breaks = c(0, 25, 50, 75, 100),
            labels = c(">0", "25", "50", "75", "100")
          ) +
          theme_bw(base_size = 14) +
          labs(
            x = "Det. freq. threshold \nin each river (%)",
            y = "Overlap (%)"
          )
      }
      
      # 
      # 4) PLOT OBJECTS PER THRESHOLD (VENN + BOXPLOT)
      # 
      
      make_targeted_venn_box <- function(df, det_tbl, threshold) {
        
        sets <- make_sets_from_det_targeted(det_tbl, threshold = threshold)
        set_spree  <- sets$set_spree
        set_neckar <- sets$set_neckar
        
        th_label <- if (is.character(threshold) && threshold == "gt0") {
          "gt0"
        } else {
          paste0(as.integer(round(threshold * 100)), "pct")
        }
        
        # --- VENN (counts) ---
        venn_plot <- ggvenn(
          list(Spree = set_spree, Neckar = set_neckar),
          show_percentage = FALSE,
          fill_color = c('#8795e8', '#60e4c1'),
          stroke_size = 0.5,
          stroke_color = "white",
          set_name_size = 4,
          text_size = 4,
          auto_scale = TRUE
        )
        
        # BOXPLOT (mean scaled concentration per sample)
        keep_ids <- union(set_spree, set_neckar)
        
        df_use <- df %>%
          filter(FeatureID_clean %in% keep_ids)
        
        summed_df <- df_use %>%
          group_by(Sample, ATTRIBUTE_RiverGroup) %>%
          summarise(SummedConcentration = mean(Concentration, na.rm = TRUE), .groups = "drop") %>%
          mutate(ATTRIBUTE_RiverGroup = factor(ATTRIBUTE_RiverGroup, levels = c("Spree", "Neckar")))
        
        box_plot <- ggplot(
          summed_df,
          aes(x = ATTRIBUTE_RiverGroup, y = SummedConcentration, fill = ATTRIBUTE_RiverGroup)
        ) +
          geom_boxplot() +
          theme_bw(base_size = 14) +
          theme(legend.position = "none") +
          scale_fill_manual(values = c("Spree" = '#8795e8', "Neckar" = '#60e4c1')) +
          labs(x = "River", y = "Average Scaled Concentration")
        
        invisible(list(
          venn_plot = venn_plot,
          box_plot = box_plot,
          th_label = th_label,
          sets = sets
        ))
      }
      
      #
      # 5) COMPILED FIGURE FOR TARGETED (ONE THRESHOLD)
      #    - Layout: venn | overlap | box
      
      compile_detection_rate_plot_targeted <- function(df_targeted, Directory, threshold,
                                                       det_tbl_targeted,
                                                       overlap_plot_targeted) {

        
        out <- make_targeted_venn_box(
          df = df_targeted,
          det_tbl = det_tbl_targeted,
          threshold = threshold
        )
        
        compiled_plot <- (out$venn_plot | overlap_plot_targeted | out$box_plot)
        
        out_fp <- file.path(
          Directory, "03_Output_Figures",
          paste0(Sys.Date(), "_CompiledPlots_TargetedCompounds_", out$th_label, "_VennOverlapBox.pdf")
        )
        
        ggsave(out_fp, compiled_plot, height = 4, width = 12)
        
        invisible(list(
          compiled_plot = compiled_plot,
          out_fp = out_fp,
          th_label = out$th_label
        ))
      }
      
      # 
      # 6) RUN: build overlap once; compile for gt0 + 50%
      #
      
      thresholds_use <- list("gt0", 0.25, 0.50, 0.75, 1.00)
      
      det_tbl_targeted <- build_det_tbl_targeted(df_targeted)
      
      overlap_df_targeted <- make_overlap_curve_targeted(det_tbl_targeted, thresholds = thresholds_use)
      print(overlap_df_targeted)
      
      # save overlap table
      curve_out <- file.path(Directory, "02_Output_DataProducts",
                             paste0(Sys.Date(), "_OverlapCurve_TargetedCompounds_Table.csv"))
      write.csv(overlap_df_targeted, curve_out, row.names = FALSE)
      
      # overlap plot object reused in BOTH compiled PDFs
      overlap_plot_targeted <- make_overlap_plot_targeted(overlap_df_targeted)
      
      compiled_targeted_gt0 <- compile_detection_rate_plot_targeted(
        df_targeted = df_targeted,
        Directory = Directory,
        threshold = "gt0",
        det_tbl_targeted = det_tbl_targeted,
        overlap_plot_targeted = overlap_plot_targeted
      )
      
      compiled_targeted_50 <- compile_detection_rate_plot_targeted(
        df_targeted = df_targeted,
        Directory = Directory,
        threshold = 0.50,
        det_tbl_targeted = det_tbl_targeted,
        overlap_plot_targeted = overlap_plot_targeted
      )
      
##### 4.4 HEATMAP of FEATURE GROUPS across SITES (FIGURE 3) #####                                                                            
      
      # Reformat scaled targeted concentration dataframe from wide to long format - Use all Samples (i.e., with Erpe)                                                                     
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
      filter(Sample %in% c("N1_both_1.mzML",     "N2_both_1.mzML",     "N3_both_1.mzML",     "N4_both_1.mzML",
                           "N5_both_1.mzML",     "N6_both_1.mzML",     "N7_both_1.mzML",     "N8_both_1.mzML",     "N9_both_1.mzML",
                           "N10_both_1.mzML",     "N11_both_1.mzML",     "N12_both_1.mzML",
                           "S1_both_1.mzML",     "S2_both_1.mzML",     "S3_both_1.mzML",     "S4_both_1.mzML",     "S5_both_1.mzML",
                           "S6_both_1.mzML",     "S7_both_1.mzML",     "S8_both_1.mzML",
                           "S9_both_1.mzML",     "S10_both_1.mzML",     "S11_both_1.mzML",     "S12_both_1.mzML",
                           "S13_both_1.mzML",     "S14_both_1.mzML",     "S15_both_1.mzML",
                           "S16_both_1.mzML",     "S17_both_1.mzML",     "S18_both_1.mzML",     "S19_both_1.mzML"))
      
      # set sample order for plot
      sample_order <- c(
        "N1_both_1.mzML", "N2_both_1.mzML", "N3_both_1.mzML", "N4_both_1.mzML",
        "N5_both_1.mzML", "N6_both_1.mzML", "N7_both_1.mzML", "N8_both_1.mzML",
        "N9_both_1.mzML", "N10_both_1.mzML", "N11_both_1.mzML", "N12_both_1.mzML",
        "S1_both_1.mzML", "S2_both_1.mzML", "S3_both_1.mzML", "S4_both_1.mzML",
        "S5_both_1.mzML", "S6_both_1.mzML", "S7_both_1.mzML", "S8_both_1.mzML",
        "S9_both_1.mzML", "S10_both_1.mzML", "S11_both_1.mzML", "S12_both_1.mzML",
        "S13_both_1.mzML", "S14_both_1.mzML", "S15_both_1.mzML",
        "S16_both_1.mzML", "S17_both_1.mzML", "S18_both_1.mzML", "S19_both_1.mzML"
      )
      
      # Pivot to wide format for ComplexHeatmap
      heatmap_matrix <- corr_df_filtered %>%
      group_by(Compound, Sample) %>%
      summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Sample, values_from = Concentration) %>%
      column_to_rownames("Compound") %>%
      as.matrix()
      
      # reorder sites by preferred order
      heatmap_matrix <- heatmap_matrix[, sample_order, drop = FALSE]
    
      # Generate heat map dataframe
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
      
      # Plot heatmap
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
                   "valsartan acid", 
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

# Add p-value signs to correlation summary
cor_summary <- cor_summary %>%
  mutate(Signif = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))



# Define output file path & write .csv
output_path <- file.path(Directory, "02_Output_DataProducts", paste0(Sys.Date(), "_Correlation_TargetedCompounds_SpearmanRho_AllSites_noErpe.csv"))
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



                                                                                
                                                                                  
                                                                                  
                                                                                  
