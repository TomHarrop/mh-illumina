#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

############
# FUNCTION #
############

FilenameToAssemblyId <- function(filename){
    my_suffix <- sub(
        "/Volumes/archive/deardenlab/tomharrop/projects/mh-illumina/output/040_meraculous/",
        "",
        filename)
    sub("/meraculous_gap_closure/final.scaffolds.fa",
        "",
        my_suffix)
}

###########
# GLOBALS #
###########

stats_file <- snakemake@input[["stats"]]
plot_file <- snakemake@output[["plot"]]
log_file <- snakemake@log[["log"]]

# dev
# stats_file <- "output/040_meraculous/stats.txt"
# plot_file <- "output/050_assembly-stats/assembly_stats.pdf"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read stats
stats <- fread(stats_file)
stats[, c("read_set", "k", "diplo") :=
          tstrsplit(FilenameToAssemblyId(filename), "_"), by = filename]

# fill in missing values during melt
setkey(stats, read_set, k, diplo)
pd <- melt(stats[CJ(unique(read_set), unique(k), unique(diplo))],
           id.vars = c("read_set", "k", "diplo"),
           measure.vars = c("n_scaffolds", "scaf_bp", "scaf_N50", "scaf_L50"),
           fill = TRUE)

# set up labels
read_set_order <- c( "trim-decon" = "Raw",
                     "norm" = "Normalised")
variable_order <- c("n_scaffolds" = "Scaffolds (K)",
                    "scaf_bp" = "Assembled size (Mb)",
                    "scaf_N50" = "N50 (K)",
                    "scaf_L50" = "L50 (Kb)")

pd[, value := as.double(value)]
pd[variable == "n_scaffolds", value := value / 1000]
pd[variable == "scaf_bp", value := value / 1e6]
pd[variable == "scaf_N50", value := value / 1000]
pd[variable == "scaf_L50", value := value / 1000]


pd[, k := as.numeric(gsub("[^[:digit:]]+", "", k))]
pd[, diplo := as.numeric(gsub("[^[:digit:]]+", "", diplo))]

pd[, read_set := factor(plyr::revalue(read_set, read_set_order),
                        levels = read_set_order)]

pd[, variable := factor(plyr::revalue(variable, variable_order),
                        levels = variable_order)]

# draw the plot
gp <- ggplot(pd, aes(x = as.factor(k), y = value, fill = as.factor(diplo))) +
    theme_minimal() +
    theme(strip.placement = "outside") +
    facet_grid(variable ~ read_set, scales = "free", switch = "y") +
    xlab(expression(italic("k"))) + ylab(NULL) +
    scale_fill_brewer(palette = "Set1",
                      guide = guide_legend(title = "Ploidy")) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.5)

# write output
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
