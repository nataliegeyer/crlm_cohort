library(tidyverse)

# Globals
base_dir <-"/home/bibu/Workspace/crlm_cohort/output"
combined_fn <- paste(base_dir, "/combined_annotations.csv", sep="")

combined_slide_fn <- paste(base_dir, "/gp_annotations_by_slide.csv", sep="")
combined_tumor_fn <- paste(base_dir, "/gp_annotations_by_tumor.csv", sep="")
combined_probe_fn <- paste(base_dir, "/gp_annotations_by_probe.csv", sep="")

regression_slide_fn <- paste(base_dir, "/regression_by_slide.csv", sep="")
regression_tumor_fn <- paste(base_dir, "/regression_by_tumor.csv", sep="")
regression_probe_fn <- paste(base_dir, "/regression_by_probe.csv", sep="")

# Read all annotations
df <- read.csv(combined_fn, row.names=NULL)

# Inv front annotations
inv_front <- df %>% select(-percents) %>% filter(annotation_types != "Tumor")

# Inv front annotations, sum and % by GP and slide
sum_inv_front_slide <-  inv_front %>% group_by(ids, tumors, blocks, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_slide <- sum_inv_front_slide %>% group_by(ids, tumors, blocks) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_slide, combined_slide_fn, row.names=FALSE)

# Inv front annotations, sum and % by GP and tumor
sum_inv_front_tumor <- inv_front %>% group_by(ids, tumors, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_tumor <- sum_inv_front_tumor %>% group_by(ids, tumors) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_tumor, combined_tumor_fn, row.names=FALSE)

# Inv front annotations, sum and % by GP and probe - summarizing the GPs in every slide (not by tumors, to diminish bias by tumor size)
sum_inv_front_probe <- inv_front %>% group_by(ids, annotation_types) %>% summarise(length_um = sum(lengths_um))

percent_inv_front_probe <- sum_inv_front_probe %>% group_by(ids) %>% mutate(percent_gp = round(100 / sum(length_um) * length_um, 2))

write.csv(percent_inv_front_probe, combined_probe_fn, row.names=FALSE)


# Tumor regression by slide
regression_slide <- df %>% select(-lengths_um) %>% filter(annotation_types == "Tumor")

write.csv(regression_slide, regression_slide_fn, row.names=FALSE)

# Tumor regression by tumor
regression_tumor <- regression_slide %>% group_by(ids, tumors) %>% summarise(avg_percent = round(mean(percents), 2))

write.csv(regression_tumor, regression_tumor_fn, row.names=FALSE)

# Tumor regression by probe - summarizing the regression in every slide (not by tumors, to diminish bias by tumor size)
regression_probe <- regression_slide %>% group_by(ids) %>% summarise(avg_percent = round(mean(percents), 2))

write.csv(regression_probe, regression_probe_fn, row.names=FALSE)