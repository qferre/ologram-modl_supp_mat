## Collect command line arguments
# NOTE: For now, it is hardcoded. The only argument expected by the script is
# its point of execution, as it will be executed seperately by a Snakefile.
args <- commandArgs(trailingOnly = TRUE)
this_dir <- toString(args[1])

setwd(this_dir)
print(getwd())


# ---------------------------------------------------------------------------- #

# Install packages
allpkgs <- c("UpSetR", "dplyr", "reshape2", "patchwork", "latex2exp")
for (p in allpkgs){
  if (!requireNamespace(p)){ install.packages(c(p), repos="http://cran.us.r-project.org", dependencies = TRUE) }
}

library('UpSetR')
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(grid)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(patchwork)
library(latex2exp)

d <- read.table("all_ologram_results.txt", head=T, sep="\t")

# Compute -log10(pval)
d$summed_bp_overlaps_pvalue[d$summed_bp_overlaps_pvalue == 0] <- 1e-320
d$log10_pval <- -log10(d$summed_bp_overlaps_pvalue)


## Extract ("regexp with capture '()'") the query name
d$query <- str_match(d$query, "-([0-9A-Za-z]+)/00_ologram.*")[,2]

## Extract the combination
tmp <- lapply(strsplit(as.character(d$feature_type), "\\+"), str_match, "_([^_]+)_")

# Quick hack
tmp <- lapply(tmp, "[", i=,j=2)
tmp <- lapply(lapply(lapply(tmp, "[", -1), rev), "[", -1)

all_features <- unique(unlist(tmp))
feature_mat <- matrix(0, nrow=nrow(d), ncol=length(all_features))
colnames(feature_mat) <- all_features
for(i in 1:length(tmp)){
  for(j in tmp[[i]]){
    feature_mat[i, j] <- 1
  }
}

## Prepare a binary matrix with TF as column and combi as row 
feature_mat <- as.data.frame(feature_mat)
feature_mat$degree <- rowSums(feature_mat) + 1
d <- cbind(d, as.matrix(feature_mat))

## Let's plot
d$color <-as.character(factor(d$query, levels=brewer.pal(3,"Paired")))

d$combination <- 1:nrow(d)

dm <- melt(data = d, id.vars = c("degree",
                                 "summed_bp_overlaps_true",
                                 "summed_bp_overlaps_pvalue",
                                 "query",
                                 "combination",
                                 "log10_pval",
                                 "summed_bp_overlaps_log2_fold_change"),
     measure.vars=all_features)
dm$value[dm$value == 0] <- NA
dm$Factors <- dm$variable



## Order combination levels by *summed basepairs in the true overlaps*
dm$combination <- factor(dm$combination, ordered = T, levels = unique(as.character(dm$combination[order(dm$summed_bp_overlaps_true, decreasing = T)])))

dm$Factors <- factor(dm$Factors, levels=c("Ctcf", "Rad21", "Smc1a", "Smc3",
                                            "Irf1", "Irf9", "Stat1", "Stat6",
                                            "Nanog","Pou5f1","Klf4", "Sox2"), order=T)


### Point colors
dm$Dataset <- as.character(dm$variable)
dm$Dataset[dm$Dataset %in% c("Ctcf", "Rad21", "Smc1a", "Smc3")] <- "Insulators"
dm$Dataset[dm$Dataset %in% c("Irf1", "Irf9", "Stat1", "Stat6")] <- "Interferon\nresponse"
dm$Dataset[dm$Dataset %in% c("Nanog","Pou5f1","Klf4", "Sox2")] <- "Stem cells"


### Take the n best combi for each query
n_best <- 20

dm %>% 
  arrange(desc(query)) %>% 
  arrange(desc(summed_bp_overlaps_true)) %>%
  group_by(query) %>% 
  slice(1:(n_best*length(all_features)-1)) %>% 
  ungroup() -> dm_sub

### Order queries
dm_sub$query <- factor(dm_sub$query, levels = c("Nanog", "Ctcf", "Irf1"),
                   ordered=TRUE)

### Plot
p1 <- ggplot(dm_sub, aes(y = Factors, x=combination,  fill=Dataset)) +
    geom_point(aes(size=value), guide='none', shape=22, color="white") +
    theme_bw() +
    scale_size(guide="none") +
    theme(axis.text.x = element_text(angle=0, size=7),
          axis.text.y = element_text(size=7),
          legend.key.height = unit(1,"line"),
          axis.title.y = element_text(size=8)) +
    scale_color_manual(values=c("#D0AB58", "#3F4EF3", "#A63965")) +
    scale_fill_manual(values=c("#D0AB58", "#3F4EF3", "#A63965")) +
    facet_wrap(~query, scale="free_x") +
    scale_x_discrete(name ="Combinations",
                     labels=c(1:n_best)) 

p2 <- ggplot(dm_sub,aes(x=combination, y = log10_pval, fill= factor(degree) )) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(size=8)) +
  scale_fill_manual(name="Combi. order", values=c("#FFC30F", "#1D9A6C", "#104765", "#08062C")) +
  facet_wrap(~query, scale="free_x", strip.position = "bottom")  

p3 <- ggplot(dm_sub,aes(x=combination, y = summed_bp_overlaps_log2_fold_change, fill = factor(degree) )) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size=8)) +
  scale_fill_manual(values=c("#FFC30F", "#1D9A6C", "#104765", "#08062C")) +
  facet_wrap(~query, scale="free_x", strip.position = "bottom") +
  ylab(TeX('$log2(\\sum bp_{obs}/ \\sum bp_{sim})$'))  
  

# Save the combined plot
combined_plot <- p3 + p2 + p1 + plot_layout(ncol = 1, heights = c(1,1,2))
ggsave("murine_fig.png", plot = combined_plot, width = 12, height = 8, dpi = 240)
