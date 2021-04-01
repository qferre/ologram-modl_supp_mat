# NOTE: this will be executed by a Snakefile, so the point of execution is the root of the entire directory
this.dir = "./output/murine_result/"
setwd(this.dir)
print(getwd())




# Install packages
if (!requireNamespace("UpSetR")){
  install.packages(c("UpSetR"))
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

## Extract the combination.
tmp <- lapply(strsplit(as.character(d$feature_type), "\\+"), str_match, "_([^_]+)_")

## Ca c'est très laid mais ça marche...
tmp <- lapply(tmp, "[", i=,j=2)

## Mon dieu :)
tmp <- lapply(lapply(lapply(tmp, "[", -1), rev), "[", -1)

all_features <- unique(unlist(tmp))
feature_mat <- matrix(0, nrow=nrow(d), ncol=length(all_features))
colnames(feature_mat) <- all_features
for(i in 1:length(tmp)){
  for(j in tmp[[i]]){
    feature_mat[i, j] <- 1
  }
}

## Prépare a binary matrix with TF as col and combi as row 
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












### Order combination levels by p-val
dm$combination <- factor(dm$combination, ordered = T, levels = unique(as.character(dm$combination[order(dm$log10_pval, dm$summed_bp_overlaps_log2_fold_change, decreasing = T)])))

















### Order combination levels by SUMMED BP OVERLAPS
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


######################################################
######################################################
######################################################
# TODO 15 should suffice
######################################################
######################################################
######################################################











dm %>% 
  arrange(desc(query)) %>% 
  arrange(desc(log10_pval)) %>%  
  group_by(query) %>% slice(1:(n_best*length(all_features)-1)) %>% ungroup() -> dm_sub

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
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.height = unit(1,"line"),
          axis.title.y = element_text(size=8)) +
    scale_color_manual(values=c("#D0AB58", "#3F4EF3", "#A63965")) +
    scale_fill_manual(values=c("#D0AB58", "#3F4EF3", "#A63965")) +
    facet_wrap(~query,  scale="free_x") +
    scale_x_discrete(name ="Combinations",
                     labels=c(1:n_best)) 

p2 <- ggplot(dm_sub,aes(x=combination, y = log10_pval)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(size=8)) +
  facet_wrap(~query, scale="free_x", strip.position = "bottom")  

p3 <- ggplot(dm_sub,aes(x=combination, y = summed_bp_overlaps_log2_fold_change)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8)) +
  facet_wrap(~query, scale="free_x", strip.position = "bottom") +
  ylab(TeX('$log2(\\sum bp_{obs}/ \\sum bp_{sim})$'))  
  


ggsave("murine_fig.png", plot = p2 / p3 / p1, width = 12, height = 8, dpi = 200)






# TODO : save separately so I can squish more the bars compared to the dots of combination presence