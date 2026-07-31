args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    "Usage: Rscript audit_psychad_public_metadata.R <metadata.csv> <output_dir>",
    call. = FALSE
  )
}

metadata_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_directory <- args[[2L]]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)

metadata <- utils::read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "DonorID",
  "Cohort",
  "Age",
  "Sex",
  "Ancestry",
  "PMI",
  "Diagnosis",
  "Tier1_crossDis",
  "Tier1_crossDis_dx"
)

missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(metadata$DonorID)) {
  stop("DonorID is not unique in the public metadata.", call. = FALSE)
}

pilot_groups <- c("CTRL", "AD", "DLBD", "SCZ")
pilot <- metadata[
  metadata$Tier1_crossDis == "Y" &
    metadata$Tier1_crossDis_dx %in% pilot_groups,
  ,
  drop = FALSE
]

pilot$group <- factor(
  pilot$Tier1_crossDis_dx,
  levels = pilot_groups
)

pilot$age_top_coded <- grepl("\\+$", pilot$Age)
pilot$age_numeric <- suppressWarnings(
  as.numeric(sub("\\+$", "", pilot$Age))
)
pilot$pmi_numeric <- suppressWarnings(as.numeric(pilot$PMI))

summarize_numeric <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_
    ))
  }

  c(
    median = stats::median(values),
    q1 = unname(stats::quantile(values, 0.25)),
    q3 = unname(stats::quantile(values, 0.75))
  )
}

sample_summary <- do.call(
  rbind,
  lapply(
    pilot_groups,
    function(group_name) {
      group_data <- pilot[pilot$group == group_name, , drop = FALSE]
      age_summary <- summarize_numeric(group_data$age_numeric)
      pmi_summary <- summarize_numeric(group_data$pmi_numeric)

      data.frame(
        group = group_name,
        donors = nrow(group_data),
        female_n = sum(group_data$Sex == "F", na.rm = TRUE),
        female_percent = round(
          100 * mean(group_data$Sex == "F", na.rm = TRUE),
          1
        ),
        age_median = age_summary[["median"]],
        age_q1 = age_summary[["q1"]],
        age_q3 = age_summary[["q3"]],
        age_top_coded_n = sum(group_data$age_top_coded, na.rm = TRUE),
        pmi_median = pmi_summary[["median"]],
        pmi_q1 = pmi_summary[["q1"]],
        pmi_q3 = pmi_summary[["q3"]],
        stringsAsFactors = FALSE
      )
    }
  )
)

count_by <- function(data, variable_name) {
  result <- as.data.frame(
    table(
      group = data$group,
      category = data[[variable_name]],
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  result[result$Freq > 0L, , drop = FALSE]
}

cohort_counts <- count_by(pilot, "Cohort")
ancestry_counts <- count_by(pilot, "Ancestry")
diagnosis_counts <- as.data.frame(
  sort(table(metadata$Diagnosis), decreasing = TRUE),
  stringsAsFactors = FALSE
)
names(diagnosis_counts) <- c("diagnosis", "donors")

audit_metrics <- data.frame(
  metric = c(
    "all_public_donors",
    "tier1_cross_disorder_donors",
    "pilot_donors",
    "pilot_unique_donor_ids",
    "pilot_top_coded_age_donors"
  ),
  value = c(
    nrow(metadata),
    sum(metadata$Tier1_crossDis == "Y"),
    nrow(pilot),
    length(unique(pilot$DonorID)),
    sum(pilot$age_top_coded)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  sample_summary,
  file.path(output_directory, "pilot_sample_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  cohort_counts,
  file.path(output_directory, "pilot_cohort_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  ancestry_counts,
  file.path(output_directory, "pilot_ancestry_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  diagnosis_counts,
  file.path(output_directory, "all_diagnosis_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  audit_metrics,
  file.path(output_directory, "audit_metrics.csv"),
  row.names = FALSE
)

cat("PsychAD public metadata audit passed.\n")
cat("All donors:", nrow(metadata), "\n")
cat("Tier-1 donors:", sum(metadata$Tier1_crossDis == "Y"), "\n")
cat("Pilot donors:", nrow(pilot), "\n")
print(sample_summary, row.names = FALSE)
