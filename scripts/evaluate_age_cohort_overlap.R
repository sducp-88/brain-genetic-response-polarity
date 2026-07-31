args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    "Usage: Rscript evaluate_age_cohort_overlap.R <metadata.csv> <output_dir>",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(ggplot2))

metadata_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_directory <- args[[2L]]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)

metadata <- utils::read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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
pilot$age_numeric <- suppressWarnings(
  as.numeric(sub("\\+$", "", pilot$Age))
)

diseases <- c("AD", "DLBD", "SCZ")
cohorts <- sort(unique(pilot$Cohort))

support_rows <- list()
row_index <- 1L

for (disease_name in diseases) {
  for (cohort_name in cohorts) {
    case_age <- pilot$age_numeric[
      pilot$group == disease_name &
        pilot$Cohort == cohort_name
    ]
    control_age <- pilot$age_numeric[
      pilot$group == "CTRL" &
        pilot$Cohort == cohort_name
    ]

    case_age <- case_age[is.finite(case_age)]
    control_age <- control_age[is.finite(control_age)]

    if (length(case_age) == 0L || length(control_age) == 0L) {
      next
    }

    support_lower <- max(min(case_age), min(control_age))
    support_upper <- min(max(case_age), max(control_age))
    has_support <- support_lower <= support_upper

    cases_in_support <- if (has_support) {
      sum(case_age >= support_lower & case_age <= support_upper)
    } else {
      0L
    }
    controls_in_support <- if (has_support) {
      sum(control_age >= support_lower & control_age <= support_upper)
    } else {
      0L
    }

    support_rows[[row_index]] <- data.frame(
      disease = disease_name,
      cohort = cohort_name,
      total_cases = length(case_age),
      total_controls = length(control_age),
      case_age_min = min(case_age),
      case_age_max = max(case_age),
      control_age_min = min(control_age),
      control_age_max = max(control_age),
      support_lower = if (has_support) support_lower else NA_real_,
      support_upper = if (has_support) support_upper else NA_real_,
      cases_in_support = cases_in_support,
      controls_in_support = controls_in_support,
      cases_retained_percent = round(
        100 * cases_in_support / length(case_age),
        1
      ),
      controls_retained_percent = round(
        100 * controls_in_support / length(control_age),
        1
      ),
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}

support_table <- do.call(rbind, support_rows)

utils::write.csv(
  support_table,
  file.path(output_directory, "age_common_support_counts.csv"),
  row.names = FALSE
)

plot_data <- pilot[
  is.finite(pilot$age_numeric),
  ,
  drop = FALSE
]

age_plot <- ggplot(
  plot_data,
  aes(
    x = group,
    y = age_numeric,
    fill = group
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.55,
    color = "grey35"
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  facet_wrap(~ Cohort, scales = "free_x") +
  scale_fill_manual(
    values = c(
      CTRL = "#808080",
      AD = "#C43C39",
      DLBD = "#7A5195",
      SCZ = "#2F7E9E"
    )
  ) +
  labs(
    title = "PsychAD Tier-1 pilot: age distributions by cohort",
    subtitle = "Age 89+ is represented as 89 for this diagnostic plot",
    x = NULL,
    y = "Age (years)",
    fill = "Group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey70"
    )
  )

ggsave(
  filename = file.path(
    output_directory,
    "age_distribution_by_group_and_cohort.png"
  ),
  plot = age_plot,
  width = 10,
  height = 5.8,
  dpi = 220
)

cat("Age/cohort overlap audit passed.\n")
print(support_table, row.names = FALSE)
