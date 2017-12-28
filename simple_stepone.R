### revised script to process stepone or quantstudio output
### - takes output from two different qPCR instruments and re-arranges summaries for database
### - uses tidyverse 

require("tidyverse")                                        # requires package for data manipulations

### USER DEFINE ###   ###   ###
fileName <- "my_file_name.txt"                              # provide qPCR file to read
machineType <- "Quantstudio"                                # Use "Quantstudio" or "StepOne"
inhibitionThreshold <- 1                                    # Ct shift for inhibition call
###   ###   ###   ###   ###   ###

if(machineType == "StepOne"){                               # read in file
  data <- read.table(fileName, sep = "\t", skip = 8, header = TRUE, na.strings = "")
} else if (machineType == "Quantstudio") {
  data <- read.table(fileName, sep = "\t", skip = 43, header = TRUE, na.strings = "")
} else {
  stop("Invalid machine type")
}

data_tibble <- as.tibble(data)                              # put data into tibble
if(machineType == "StepOne"){
  data_named <- rename(data_tibble, Sample = Sample.Name,   # if StepOne column rename
                       Target = Target.Name, Ct = CÑ.) 
  } else{
    data_named <- rename(data_tibble, Sample = Sample.Name, # if QuantStudio column rename
           Target = Target.Name, Ct = CT)
  }

### NEED TO FIX no amp to NA for quantity issue! Currently treat no amps as missing data
data_relevant <- data_named %>%
  select(Sample, Target, Ct, Quantity, Task) %>%            # extract only needed columns
  mutate(., Ct = replace(Ct, Ct=="Undetermined", NA)) %>%   # convert "Undetermined" to NA
  mutate(., Ct = as.numeric(as.character(Ct)),              # make Ct numeric
         Quantity = as.numeric(Quantity)) %>% 
  group_by(Sample, Target) %>%                              # group by sample and target
  summarise(AMP = sum(!is.na(Ct)),                          # amplifications
            Mean_CT = mean(Ct, na.rm = T),                  # mean Ct
            Mean_Quant = mean(Quantity, na.rm = T),         # mean quant (note: no amp to NA)
            SD_QUANT = sd(Quantity, na.rm = T),             # SD quant (note: no amp to NA)
            NTC = sum(Task == "NTC") > 0)                   # flag NTC samples

data_taxon <- data_relevant %>%                             # extract taxon assay data
  filter(., Target != "IPC")

data_IPC <- data_relevant %>%                               # extract IPC Ct data
  filter(., Target == "IPC") %>%
  select(., Sample, Mean_CT) %>%
  rename(., IPCCT = Mean_CT)

data_together <- left_join(x = data_taxon, y = data_IPC)    # bind IPC mean Ct to taxon data

IPC_thresholds <- data_together %>%                         # IPC thresholds by target
  filter(., NTC == T) %>%
  group_by(Target) %>%
  select(., IPCCT) %>%
  rename(., IPCTHR = IPCCT)

data_thresh <- data_together %>%                            # data with inhibition indicated
  left_join(., y = IPC_thresholds) %>%                      # if no IPC, simply NA for "INHIB"
  mutate(., INHIB = IPCCT - IPCTHR > inhibitionThreshold) %>%
  filter(., NTC == FALSE) %>%
  select(., -NTC, -IPCCT, -IPCTHR)  %>%
  mutate(., Machine = "stepone") %>%                        # add machine name
  mutate(., Plate = sub("(\\w+)\\.\\w+", "\\1", fileName)) %>%  # add file name
  select(., Machine, Plate, Sample, Target, Mean_CT, AMP, Mean_Quant, SD_QUANT, INHIB) # order columns

write.csv(data_thresh,                                      # output file
          paste("reformat_",sub("(\\w+)\\.\\w+", "\\1",fileName),".csv", sep=""))




