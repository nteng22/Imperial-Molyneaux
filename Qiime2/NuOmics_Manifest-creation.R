setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT007_PROFOUND/Original-files/")

library(tidyverse)
library(readxl)

# Sequencing list
sequencing <- readxl::read_xlsx("2 - Sample List - PFND, MTN, BRU, ILDCON, PCMS.xlsx", trim_ws = TRUE)
sequencing <- sequencing %>%
  select(1:6)
colnames(sequencing) <- c("PatientID", "Sample-id", "Visit", "SampleType", "Fibrosis", "Diagnosis")
sequencing$`Sample-id` <- gsub("_", "-", sequencing$`Sample-id`)
sequencing$`Sample-id` <- gsub("_(.*)$", "", sequencing$`Sample-id`)

# Manifest created on HPC using command line, to check if sampleIDS match
manifest <- read_tsv("Manifest-Map-final.txt", trim_ws = TRUE) %>%
  select(-`...4`)
manifest$`sequencing-id` <- gsub("^([^_]*_[^_]*)_.*$", "\\1", manifest$`sample-id`)
manifest$`sample-id` <- gsub("_(.*)$", "", manifest$`sample-id`)
manifest_id <- manifest$`sample-id`
missingID <- setdiff(manifest_id, sequencing$`Sample-id`) # controls from sequencing
missingID <- missingID[!str_detect(missingID, "PFND") & !str_detect(missingID, "MTN")]

# Create df with missing IDs 
missingID_visit <- rep(NA, length(missingID)) 
missingID_SampleType <- rep("Sequencing Control", length(missingID))
missingID_Diagnosis <- rep("Control", length(missingID))

controls_df <- data.frame("PatientID" = "Controls",
                          "sample-id" = missingID, 
                          "Visit" = missingID_visit,
                          "SampleType" = missingID_SampleType,
                          "Fibrosis" = "Control",
                          "Diagnosis" = missingID_Diagnosis)
colnames(controls_df) <- colnames(sequencing)

sequencing_df <- rbind(controls_df, sequencing)

duplicated_list <- sequencing_df %>%
  filter(str_detect(`Sample-id`, "duplicate"))

#### Change Diagnosis to Control for samples where it is a control
sequencing_df <- sequencing_df %>%
  mutate(Diagnosis = ifelse(str_detect(SampleType, "Flush"), "Control", Diagnosis),
         Diagnosis = ifelse(str_detect(SampleType, "Control"), "Control", Diagnosis),
         Diagnosis = ifelse(str_detect(`Sample-id`, "POS"), "Positive Control", Diagnosis),
         Diagnosis = ifelse(str_detect(`Sample-id`, "NEG"), "Negative Control", Diagnosis),
         Fibrosis = ifelse(str_detect(SampleType, "Control"), "Control", Fibrosis),
         PatientID = ifelse(str_detect(SampleType, "Control"), "Control", PatientID))

sequencing_df %>%
  group_by(Diagnosis) %>%
  dplyr::summarise(n = n()) 
  
#### Add the picogreen from NU to the list ####
file_path <- "Picogreen/20240613_Imperial-gDNA_PICO.xlsx"
sheet_names <- excel_sheets(file_path)
picogreen <- data.frame()
for (sheet in sheet_names) {
  df <- read_excel(file_path, sheet = sheet)
  df$sheet <- sheet
  picogreen <- rbind(df, picogreen)
  }

picogreen$`Sample ID` <- gsub("_", "-", picogreen$`Sample ID`)
setdiff(picogreen$`Sample ID`, sequencing_df$`Sample-id`)
  # some samples have gaps in between the IDs... 

picogreen$`Sample ID` <- gsub(" ", "", picogreen$`Sample ID`)
setdiff(picogreen$`Sample ID`, sequencing_df$`Sample-id`)

picogreen <- dplyr::rename("Sample-id" = "Sample ID", picogreen)
picogreen <- picogreen %>% 
  select(-Well)

manifest_picogreen <- left_join(sequencing_df, picogreen)
manifest_picogreen <- manifest_picogreen %>%
  distinct()
manifest_picogreen <- dplyr::rename("sample-id" = "Sample-id", manifest_picogreen)

## Remove duplicates from picogreen
#duplicated_picogreen <- manifest_picogreen[duplicated(manifest_picogreen$`sample-id`),]
#duplicated_ID <- duplicated_picogreen$`sample-id`
#duplicated_picogreen <- manifest_picogreen %>% #want to remove duplicates from the picogreen sheet.
#  filter(`sample-id` %in% duplicated_ID)
#manifest_picogreen <- manifest_picogreen %>%
#  filter(!`sample-id` %in% duplicated_ID)

#duplicated_manifest <- manifest %>%
#  filter(`sample-id` %in% duplicated_ID)

# Manually curate the spreadsheet according to which plate and lane it was submitted in...
#write_csv(duplicated_picogreen, "Duplicated_picogreen.csv")
#write_csv(duplicated_manifest, "Duplicated_manifest.csv")

# Manually fixed the files, also changed it on HPC. 
# New manifest includes the ending -a or -b so previous block won't work. 
duplicated_fixed <- read_csv("Duplicated_fixed_sheet.csv") %>% select(1:12)
duplicated_meta <- duplicated_fixed %>% 
  select(`sample-id...1`, PatientID, Visit, SampleType, Fibrosis, Diagnosis, `DNA (ng/ul)`, sheet) %>%
  dplyr::rename("sample-id" = "sample-id...1")

# Join duplicated samples with exsiting picogreen sheet
manifest_picogreen <- rbind(duplicated_meta, manifest_picogreen) %>%
  filter(!str_detect(`sample-id`, "duplicate"))

# Record if it has been sequenced already
manifest_picogreen <- manifest_picogreen %>%
  mutate(Sequenced = ifelse(`sample-id` %in% manifest_id, "yes", "no"))
  ## Samples not sequenced are all duplicates.
write_csv(manifest_picogreen, "Sequencing_metadata_PROFOUND.csv")

## Keep track of samples not sequenced in the end. 
manifest_not_sequenced <- manifest_picogreen %>%
  filter(Sequenced=="no" & !is.na(sheet))

gplots::venn(list(duplicated=duplicated_list$PatientID, manifest_not_sequenced=manifest_not_sequenced$PatientID))

not_sequenced <- readxl::read_xlsx("2 - Sample List - PFND, MTN, BRU, ILDCON, PCMS.xlsx", sheet=3,
                                   trim_ws = TRUE)
not_sequenced <- not_sequenced[,1:4]
colnames(not_sequenced) <- c("PatientID", "sample-id", "Visit", "SampleType")
not_sequenced$`sample-id` <- gsub("_","-", not_sequenced$`sample-id`)

gplots::venn(list(not_sequenced=not_sequenced$`sample-id`, manifest_not_sequenced=manifest_not_sequenced$`sample-id`))

#### Merge barcodes with sequencing ####
final_manifest <- left_join(manifest, manifest_picogreen, by = "sample-id")
final_manifest <- final_manifest[!str_detect(final_manifest$`sample-id`, "^ILDCON(\\d{4})a.BAL"),]

final_manifest <- final_manifest %>%
  mutate(PatientID = ifelse(str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"),
                             str_replace(`sample-id`, "ILDCON(\\d{4})\\.BAL", "ILDCON\\1"), PatientID),
         Visit = ifelse(is.na(Visit), 1, Visit),
         Diagnosis = ifelse(str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"), "Healthy", Diagnosis),
         SampleType = ifelse(str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"), "BAL", SampleType),
         Fibrosis = ifelse(str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"), "No", Fibrosis))

gplots::venn(list(picogreen=manifest_picogreen$`sample-id`, manifest=manifest$`sample-id`))
setdiff(manifest$`sample-id`, manifest_picogreen$`sample-id`) # No controls in manifest picogreen (no picogreen done). 
setdiff(manifest_picogreen$`sample-id`, manifest$`sample-id`) # Duplicated saliva and stool samples

manifest_stool <- final_manifest %>% 
  filter(SampleType == "Stool" | (SampleType == "Sequencing Control" & !str_detect(`sample-id`,"Saliva"))) %>%
  filter(!str_detect(`sample-id`,"SALIVA"))
write.table(manifest_stool, "Manifest-Map_stool.txt", sep="\t", 
            row.names = FALSE, col.names=TRUE, quote = FALSE)

manifest_stool_baseline <- manifest_stool %>%
  filter(Visit == 1 & SampleType=="Stool" & Diagnosis == "IPF" | Diagnosis == "Healthy" & Visit == 1)
write.table(manifest_stool_baseline, "Manifest-baseline.txt", sep = "\t",
            row.names = F, col.names = T, quote=F)

manifest_saliva <- final_manifest %>%
  filter(SampleType == "Saliva" | (SampleType == "Sequencing Control" & !str_detect(`sample-id`, "Stool"))) %>%
  filter(!str_detect(`sample-id`, "STOOL")) 
write.table(manifest_saliva, "Manifest-Map_saliva.txt", sep="\t", 
            row.names = FALSE, col.names=TRUE, quote = FALSE)

### Read in ddPCR file ####
ddPCR <- readxl::read_xlsx("../../NT001_Microbiota-analysis/Original-files/PFND, BRU, historic ILDCON ddPCR.xlsx")
ddPCR$`Sample ID` <- gsub("_", "-", ddPCR$`Sample ID`)
colnames(ddPCR) <- c("sample-id", "DNA (ng/ul)")
ddPCR[duplicated(ddPCR$`sample-id`),]
ddPCR <- ddPCR %>%
  mutate(`DNA (ng/ul)` = ifelse(`sample-id`=="PROCAMS110048-OR", mean(`DNA (ng/ul)`), `DNA (ng/ul)`))
ddPCR <- unique(ddPCR)

manifest_BAL <- final_manifest %>%
  filter(str_detect(SampleType, "BAL") | SampleType=="Sequencing Control" | SampleType=="Reagent Control") %>%
  filter(!str_detect(`sample-id`, "Saliva") & !str_detect(`sample-id`, "SALIVA") & 
           !str_detect(`sample-id`, "Stool") & !str_detect(`sample-id`, "STOOL")) %>%
  select(-`DNA (ng/ul)`)

manifest_BAL <- right_join(ddPCR,manifest_BAL)
manifest_BAL <- manifest_BAL %>%
  mutate(`DNA (ng/ul)` = ifelse(Diagnosis=="Positive Control", 2, `DNA (ng/ul)`),
         `DNA (ng/ul)` = ifelse(Diagnosis=="Negative Control", 0, `DNA (ng/ul)`))
ILDCON_samples <- manifest_BAL %>%
  filter(str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"))
manifest_BAL <- manifest_BAL %>%
  filter(!str_detect(`sample-id`, "ILDCON(\\d{4})\\.BAL"))

POSTCODE_ddPCR <- readxl::read_xlsx("../../NT002_POSTCODE/ddPCR_IDs.xlsx")
POSTCODE_ddPCR <- dplyr::rename("DNA (ng/ul)" = "16S copies/ml", POSTCODE_ddPCR)
POSTCODE_ddPCR <- POSTCODE_ddPCR %>%
  filter(`sample-id` %in% ILDCON_samples$`sample-id`) %>% select(`sample-id`, `DNA (ng/ul)`)

ILDCON_samples$`DNA (ng/ul)` <- POSTCODE_ddPCR$`DNA (ng/ul)`

manifest_BAL <- rbind(manifest_BAL, ILDCON_samples)
manifest_BAL <- unique(manifest_BAL)

write.table(manifest_BAL, "Manifest-Map_BAL.txt", sep="\t",
            row.names=FALSE, col.names=T, quote=F)

manifest_OR <- final_manifest %>%
  filter(str_detect(SampleType, "Oral Rinse") | SampleType=="Sequencing Control" | SampleType=="Reagent Control") %>% 
  filter(!str_detect(`sample-id`, "Saliva") & !str_detect(`sample-id`, "SALIVA") &
           !str_detect(`sample-id`, "Stool") & !str_detect(`sample-id`, "STOOL")) %>%
  select(-`DNA (ng/ul)`)

manifest_OR <- right_join(ddPCR,manifest_OR)
manifest_OR <- manifest_OR %>%
  mutate(`DNA (ng/ul)` = ifelse(Diagnosis=="Positive Control", 2, `DNA (ng/ul)`),
         `DNA (ng/ul)` = ifelse(Diagnosis=="Negative Control", 0, `DNA (ng/ul)`))
write.table(manifest_OR, "Manifest-Map_OR.txt", sep="\t",
            row.names=FALSE, col.names=T, quote=F)
