#library####
library(ape)
library(phytools)
library(betapart)
library(dendextend)
library(picante)
library(ggrepel)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(canaper)
library(reshape2)
library(tidyverse)
library(viridis)

# Taxon coverage ####
## Genetically sampled tips iaa vs jetz vs ours ####

#1515, 1018

# df <- data.frame(
#   jetz = c(896, 567),
#   jetzi = c(1380, 957),
#   ours = c(1390, 1006),
#   iaa = c(1795, 1220),
#   area = c("All", "Single-area endemics")
# )

df <- data.frame(
  gen = c(896, 567, 1390, 1006),
  all = c(1380, 957, 1633, 1220),
  label = c("All (2012)", "Single-area endemics (2012)", "All (2023)", "Single-area endemics (2023)")
)
df$diff1 <- round(((df$gen)/df$all), 3)
# df$diff1 <- round(((df$ours - df$jetz)/df$jetz), 1)
# df$diff2 <- round(((df$iaa - df$ours)/df$ours), 1)

# ggplot(df, aes(x = area)) +
#   geom_bar(aes(y=iaa, fill = "iaa"), position = "dodge", stat = "identity") +
#   geom_bar(aes(y=ours, fill = "ours"), position = "dodge", stat = "identity") +
#   geom_bar(aes(y=jetz, fill = "jetz"), position = "dodge", stat = "identity") +
#   geom_text(aes(y=jetz+50, label = paste0(jetz))) +
#   geom_text(aes(y=ours+50, label = paste0(ours)), color = "white") +
#   geom_text(aes(y=iaa+50, label = paste0(iaa))) +
#   labs(y = "No. of genetically sampled tips", x ="") +
#   scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000), expand = c(0,0)) +
#   scale_fill_grey() +
#   theme_minimal() +
#   theme(  # Remove legend
#     axis.line = element_line(color = "black"),  # Add black line to axes
#     axis.text = element_text(size = 12),  # Adjust text size
#     axis.title = element_text(size = 14),
#     axis.title.x = element_text(margin = margin(t = 15)),
#     axis.title.y = element_text(margin = margin(r = 15)))

ggplot(df, aes(x = label)) +
  geom_bar(aes(y=all, fill = "All"),  position = "dodge", stat = "identity") +
  geom_bar(aes(y=gen, fill = "Genetically sampled"),  position = "dodge", stat = "identity") +
  geom_text(aes(y=all+50, label = paste0(all)), size = 5) +
  #geom_text(aes(y=gen+50, label = paste0(gen))) +
  geom_text(aes(y=gen+50, label = paste0(scales::percent(diff1))), size = 5) +
  labs(y = "No. of tips", x ="") +
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000), expand = c(0,0)) +
  scale_fill_manual(values = c("All" = "#FDE725FF", "Genetically sampled" = "#21908CFF"),
                    labels = c("All", "Genetically sampled"), name = "") +
  theme_minimal() +
  theme(  # Remove legend
    legend.position = "top",
    axis.line = element_line(color = "black"),  # Add black line to axes
    axis.text = element_text(size = 12),  # Adjust text size
    axis.title = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    legend.text = element_text(size = 12))

# Genetic Coverage ####
# Define the data
data_text <- "Loci 0 1 2 3 4 5 6 7 8
0 234 2 0 12 0 4 1 1 0
1 120 56 42 36 10 1 3 1 0
2 39 91 74 43 23 28 7 2 1
3 17 37 56 93 55 34 16 3 1
4 16 14 28 16 19 23 29 4 341"

# Convert the text data to a matrix
gen_mx <- as.matrix(read.table(text = data_text, header = TRUE, row.names = 1))
colnames(gen_mx) <-  c(0:8)
#gen_mx[1,1] <- 0
#gen_mx[5,9] <- 0

gen_df <- as.data.frame(as.table(gen_mx))
colnames(gen_df) <- c("MT", "NUC", "Value")

ggplot(data = gen_df, aes(x=MT, y = NUC, fill = Value, label = Value)) +
  geom_tile() +
  geom_text(size =3) +  
  scale_fill_gradient2(high = "lightblue", low = "white", name ="Species") +
  labs(title = "", x="Mitochondrial loci sampled", y="Nuclear loci sampled") +
  theme_minimal() + 
  coord_fixed() +
  theme(axis.line = element_line(color = "black"),  # Add black line to axes
        axis.text = element_text(size = 10),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))

# Trees used ####
## Imputed Tree ####
tr_i <- read.tree("new_all_pass_jetz_sp.NEWICK")

outgroup <- c("Falco_berigora", "Aquila_clanga", "Phodilus_badius",
              "Tadorna_variegata", "Gallus_varius", "Larus_novaehollandiae")

tr_i <- drop.tip(tr_i, outgroup)

tr_i_tips <- tr_i$tip.label

i_sp_mx <- read.csv("jetz_sp_mx.csv")

i_sp_mx <- i_sp_mx[,colnames(i_sp_mx) %in% tr_i_tips]

rownames(i_sp_mx) <- c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA")

#genus and family matrix 
j_g_mx <- read.csv("jetz_g_mx.csv")
j_f_mx <- read.csv("jetz_f_mx.csv")

## Genetic subset of Jetz tree ####
tr_j <- read.tree("new_gen_pass_jetz_sp.NEWICK")

outgroup <- c("Falco_berigora", "Aquila_clanga", "Phodilus_badius",
              "Tadorna_variegata", "Gallus_varius", "Larus_novaehollandiae")

tr_j <- drop.tip(tr_j, outgroup)

tr_j_tips <- tr_j$tip.label

j_sp_mx <- read.csv("jetz_sp_mx.csv")

j_sp_mx <- j_sp_mx[,colnames(j_sp_mx) %in% tr_j_tips]

rownames(j_sp_mx) <- c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA")

## Our tree ####
tr <- read.tree("full_scaled_ut.tree")

tr_tips <- tr$tip.label

geog_data <- read.csv("pass_list_geog.csv")

geog_data <- subset(geog_data, geog_data$Species %in% tr_tips)

colnames(geog_data)

tr <- keep.tip(tr, geog_data$Species)

#write.csv(geog_data, "pass_geog_fl.csv")

sp_mx <- data.frame(t(geog_data[,-13]))
colnames(sp_mx) <- sp_mx[1,]
sp_mx <- sp_mx[c(-1,-2),]
rownames(sp_mx) <- c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA")

sp_mx[] <- lapply(sp_mx, as.numeric)

g_mx <- read.csv("f_mx.csv")
f_mx <- read.csv("g_mx.csv")

# Betadiversity ####
sp_beta <- beta.pair(sp_mx, index.family = "jaccard")
sp_beta

sp_phy_beta <- phylo.beta.pair(sp_mx, tr, index.family = "jaccard")
sp_phy_beta

betadiv.sum <- rbind(phylo.beta.multi(sp_mx, tr, index.family = "jaccard"),
                     beta.multi(sp_mx, index.family = "jaccard"), 
                     beta.multi(g_mx, index.family = "jaccard"), 
                     beta.multi(f_mx, index.family = "jaccard"))
betadiv.sum <- unlist(betadiv.sum)
betadiv.df <- as.data.frame(matrix(betadiv.sum, nrow = 4, ncol = 3, byrow = FALSE))

betadiv.df$level <- c("Phylo-sp","Species", "Genus", "Family")
betadiv.df$level <- factor(betadiv.df$level, levels = c("Phylo-sp","Species", "Genus", "Family"))

colnames(betadiv.df) <- c("beta.JTU", "beta.JNE", "beta.JAC", "level")

j_sp_beta <- beta.pair(j_sp_mx, index.family = "jaccard")
j_sp_beta

j_sp_phy_beta <- phylo.beta.pair(j_sp_mx, tr_j, index.family = "jaccard")
j_sp_phy_beta

j.betadiv.sum <- rbind(phylo.beta.multi(j_sp_mx, tr_j, index.family = "jaccard"),
                       beta.multi(j_sp_mx, index.family = "jaccard"), 
                       beta.multi(j_g_mx, index.family = "jaccard"), 
                       beta.multi(j_f_mx, index.family = "jaccard"))
j.betadiv.sum <- unlist(j.betadiv.sum)
j.betadiv.df <- as.data.frame(matrix(j.betadiv.sum, nrow = 4, ncol = 3, byrow = FALSE))

j.betadiv.df$level <- c("Phylo-sp","Species", "Genus", "Family")
j.betadiv.df$level <- factor(j.betadiv.df$level, levels = c("Phylo-sp","Species", "Genus", "Family"))

colnames(j.betadiv.df) <- c("beta.JTU", "beta.JNE", "beta.JAC", "level")

i_sp_beta <- beta.pair(i_sp_mx, index.family = "jaccard")
i_sp_beta

i_sp_phy_beta <- phylo.beta.pair(i_sp_mx, tr_i, index.family = "jaccard")
i_sp_phy_beta

i.betadiv.sum <- rbind(phylo.beta.multi(i_sp_mx, tr_i, index.family = "jaccard"),
                       beta.multi(i_sp_mx, index.family = "jaccard"), 
                       beta.multi(j_g_mx, index.family = "jaccard"), 
                       beta.multi(j_f_mx, index.family = "jaccard"))
i.betadiv.sum <- unlist(i.betadiv.sum)
i.betadiv.df <- as.data.frame(matrix(i.betadiv.sum, nrow = 4, ncol = 3, byrow = FALSE))

i.betadiv.df$level <- c("Phylo-sp","Species", "Genus", "Family")
i.betadiv.df$level <- factor(i.betadiv.df$level, levels = c("Phylo-sp","Species", "Genus", "Family"))

colnames(i.betadiv.df) <- c("beta.JTU", "beta.JNE", "beta.JAC", "level")

betadiv.c <- 
  ggplot(betadiv.df, aes(x = level, y = beta.JAC, fill = level)) +
  geom_bar(stat = "identity") +
  labs(title = "c) New tree", ylab = expression(beta[jac])) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.1), expand = c(0,0)) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

betadiv.b <- 
  ggplot(j.betadiv.df, aes(x = level, y = beta.JAC, fill = level)) +
  geom_bar(stat = "identity") +
  labs(title = "b) Jetz (genetic)", ylab = expression(beta[jac])) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.1), expand = c(0,0)) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

betadiv.a <- 
  ggplot(i.betadiv.df, aes(x = level, y = beta.JAC, fill = level)) +
  geom_bar(stat = "identity") +
  labs(title = "a) Jetz (imputed)", ylab = expression(beta[jac])) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1.1), expand = c(0,0)) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

grid.arrange(betadiv.a, betadiv.b, betadiv.c, ncol =3)

## Dendrogram plots per area ####
par(mfrow = c(2,3))
par(mar=c(3,5,2,2))
region = data.frame(number = 1:10)
region$region <- c("Sahul","Sahul","Island", "Island", "Island", "Island", "Sunda", "Sunda", "Sunda", "Sunda")
regiontype <- factor(region$region)
n_region_types <- length(unique(regiontype))
cols_6 <- c("#cccccc", "#fde725", "#55cc66")
col_region_type <- cols_6[regiontype]

hc.si <- hclust(i_sp_beta$beta.jac)
dend.si <- as.dendrogram(hc.si)
dend.si <- reorder(dend.si, region$number)
plot(dend.si, ylab = expression(beta[jac]))
abline(a = mean(i_sp_beta$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
title(main = "a) Species")

hc.sj <- hclust(j_sp_beta$beta.jac)
dend.sj <- as.dendrogram(hc.sj)
dend.sj <- reorder(dend.sj, region$number)
plot(dend.sj, ylab = expression(beta[jac]))
abline(a = mean(j_sp_beta$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")

hc.s <- hclust(sp_beta$beta.jac)
dend.s <- as.dendrogram(hc.s)
dend.s <- reorder(dend.s, region$number)
plot(dend.s, ylab = expression(beta[jac]))
abline(a = mean(sp_beta$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")

hc.pi <- hclust(i_sp_phy_beta$phylo.beta.jac)
dend.pi <- as.dendrogram(hc.pi)
dend.pi <- reorder(dend.pi, region$number)
plot(dend.pi, ylab = expression(beta[jac]))
abline(a = mean(i_sp_phy_beta$phylo.beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
title(main = "b) Phylo-species")

hc.pj <- hclust(j_sp_phy_beta$phylo.beta.jac)
dend.pj <- as.dendrogram(hc.pj)
dend.pj <- reorder(dend.pj, region$number)
plot(dend.pj, ylab = expression(beta[jac]))
abline(a = mean(j_sp_phy_beta$phylo.beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")

hc.p <- hclust(sp_phy_beta$phylo.beta.jac)
dend.p <- as.dendrogram(hc.p)
dend.p <- reorder(dend.p, region$number)
plot(dend.p, ylab = expression(beta[jac]))
abline(a = mean(sp_phy_beta$phylo.beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")

## Betadiversity pairwise matrix differences ####
diff.mx1 <- as.matrix(sp_beta$beta.jac-j_sp_beta$beta.jac)
diff.mx1[upper.tri(diff.mx1)] <- 0
diff.mx1 <- melt(diff.mx1)
diff.mx1$value <- diff.mx1$value*100
head(diff.mx1)

diff.mx1.p <- 
  ggplot(data = diff.mx1, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-11,11), name ="% Difference") +
  labs(title = "a) New tree - Jetz (genetic)", x="", y="") +
  theme_minimal() + 
  coord_fixed()

diff.mx2 <- as.matrix(sp_beta$beta.jac-i_sp_beta$beta.jac)
diff.mx2[upper.tri(diff.mx2)] <- 0
diff.mx2 <- melt(diff.mx2)
diff.mx2$value <- diff.mx2$value*100
head(diff.mx2)
diff.mx2.p <- 
  ggplot(data = diff.mx2, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "b) New tree - Jetz (imputed)", x="", y="") +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-11,11), name = "% Difference") +
  theme_minimal() + 
  coord_fixed()

diff.mx3 <- as.matrix(sp_phy_beta$phylo.beta.jac-j_sp_phy_beta$phylo.beta.jac)
diff.mx3[upper.tri(diff.mx3)] <- 0
diff.mx3 <- melt(diff.mx3)
diff.mx3$value <- diff.mx3$value*100
head(diff.mx3)

diff.mx3.p <- 
  ggplot(data = diff.mx3, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-11,11), name ="% Difference") +
  labs(title = "c) New tree - Jetz (genetic)", x="", y="") +
  theme_minimal() + 
  coord_fixed()

diff.mx4 <- as.matrix(sp_phy_beta$phylo.beta.jac-i_sp_phy_beta$phylo.beta.jac)
diff.mx4[upper.tri(diff.mx4)] <- 0
diff.mx4 <- melt(diff.mx4)
diff.mx4$value <- diff.mx4$value*100
head(diff.mx4)
diff.mx4.p <- 
  ggplot(data = diff.mx4, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "d) New tree - Jetz (imputed)", x="", y="") +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-11,11), name = "% Difference") +
  theme_minimal() + 
  coord_fixed()

grid.arrange(arrangeGrob(diff.mx1.p, top = "Taxonomic beta diversity"), 
             arrangeGrob(diff.mx2.p, top = ""), 
             arrangeGrob(diff.mx3.p, top ="Phylogenetic beta diversity"), 
             arrangeGrob(diff.mx4.p, top = ""),
             ncol =2, widths = c(1, 1))

diff.mx5 <- as.matrix(sp_beta$beta.jac-sp_phy_beta$phylo.beta.jac)
diff.mx5[upper.tri(diff.mx5)] <- 0
diff.mx5 <- melt(diff.mx5)
diff.mx5$value <- diff.mx5$value*100
head(diff.mx5)
diff.mx5.p <- 
  ggplot(data = diff.mx5, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "c) New tree TBD - PBD", x="", y="") +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0, 31), name = "% Difference") +
  theme_minimal() + 
  coord_fixed()
diff.mx5.p

beta.sp <- as.matrix(sp_beta$beta.jac)
beta.sp[upper.tri(beta.sp)] <- 0
beta.sp <- melt(beta.sp)
betasp2 <- 
  ggplot(data = beta.sp, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value,2), vjust = 1)) +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white", limit = c(0, 1)) +
  labs(title = "a) New tree TBD", x="", y="") +
  theme_minimal() + 
  coord_fixed()

beta.p.sp <- as.matrix(sp_phy_beta$phylo.beta.jac)
beta.p.sp[upper.tri(beta.p.sp)] <- 0
beta.p.sp <- melt(beta.p.sp)
betasp3 <- ggplot(data = beta.p.sp, aes(x=Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value,2), vjust = 1)) +
  scale_fill_gradient2(low ="blue", high = "red", mid = "white", limit = c(0, 1)) +
  labs(title = "b) New tree PBD", x="", y="") +
  theme_minimal() + 
  coord_fixed()

grid.arrange(betasp2, betasp3, diff.mx5.p, ncol =2)

# PE and PD calculations ####
## Canaper simulations for PE endemism category ####

rand_test_results_i <- cpr_rand_test(i_sp_mx, tr_i, null_model ="curveball")
canape_res_i <- cpr_classify_endem(rand_test_results_i)
# canape_res_i[, "endem_type", drop = FALSE]
# res_sum_i <- canape_res_i[, "endem_type", drop = FALSE]
# 
rand_test_results <- cpr_rand_test(sp_mx, tr, null_model ="curveball")
canape_res <- cpr_classify_endem(rand_test_results)
# canape_res[, "endem_type", drop = FALSE]
# res_sum <- canape_res[, "endem_type", drop = FALSE]
# 
rand_test_results_j <- cpr_rand_test(j_sp_mx, tr_j, null_model ="curveball")
canape_res_j <- cpr_classify_endem(rand_test_results_j)
# canape_res_j[, "endem_type", drop = FALSE]
# res_sum_j <- canape_res_j[, "endem_type", drop = FALSE]
# 
# for (i in 1:99){
#   rand_test_results <- cpr_rand_test(sp_mx, tr, null_model ="curveball")
#   canape_res <- cpr_classify_endem(rand_test_results)
#   res_sum <- cbind(res_sum, canape_res[, "endem_type", drop = FALSE])
# 
#   rand_test_results_j <- cpr_rand_test(j_sp_mx, tr_j, null_model ="curveball")
#   canape_res_j <- cpr_classify_endem(rand_test_results_j)
#   res_sum_j <- cbind(res_sum_j, canape_res_j[, "endem_type", drop = FALSE])
# 
#   rand_test_results_i <- cpr_rand_test(i_sp_mx, tr_i, null_model ="curveball")
#   canape_res_i <- cpr_classify_endem(rand_test_results_i)
#   res_sum_i <- cbind(res_sum_i, canape_res_i[, "endem_type", drop = FALSE])
# }
# 
# endem_c <- c("paleo", "super", "neo", "mixed", "not significant")
# 
# res_mx_i<- matrix(0, nrow = nrow(res_sum_i), ncol = 5)
# 
# for (i in 1:length(endem_c)) {
#   character_to_count <- endem_c[i]
# 
#   # Loop through each row
#   for (j in 1:nrow(res_sum_i)) {
#     # Extract the values for the current row
#     row_values <- res_sum_i[j, ]
# 
#     # Count how many times the specified character appears in the row
#     character_count <- sum(row_values == character_to_count)
# 
#     # Store the count in the result matrix
#     res_mx_i[j, i] <- character_count
#   }
# }
# 
# rownames(res_mx_i) <- rownames(res_sum_i)
# colnames(res_mx_i) <- endem_c
# 
# res_mx<- matrix(0, nrow = nrow(res_sum), ncol = 5)
# 
# for (i in 1:length(endem_c)) {
#   character_to_count <- endem_c[i]
# 
#   # Loop through each row
#   for (j in 1:nrow(res_sum)) {
#     # Extract the values for the current row
#     row_values <- res_sum[j, ]
# 
#     # Count how many times the specified character appears in the row
#     character_count <- sum(row_values == character_to_count)
# 
#     # Store the count in the result matrix
#     res_mx[j, i] <- character_count
#   }
# }
# 
# rownames(res_mx) <- rownames(res_sum)
# colnames(res_mx) <- endem_c
# 
# res_mx_j <- matrix(0, nrow = nrow(res_sum_j), ncol = 5)
# 
# for (i in 1:length(endem_c)) {
#   character_to_count <- endem_c[i]
# 
#   # Loop through each row
#   for (j in 1:nrow(res_sum_j)) {
#     # Extract the values for the current row
#     row_values <- res_sum_j[j, ]
# 
#     # Count how many times the specified character appears in the row
#     character_count <- sum(row_values == character_to_count)
# 
#     # Store the count in the result matrix
#     res_mx_j[j, i] <- character_count
#   }
# }
# 
# rownames(res_mx_j) <- rownames(res_sum_j)
# colnames(res_mx_j) <- endem_c
# 
# res_mx
# res_mx_j
# res_mx_i

## SR calculation (and extra PD) ####
i_area_pd <- pd(i_sp_mx, tr_i)
i_area_pd <- cbind(i_area_pd, region$region)
i_area_pd$area <- rownames(i_area_pd)
rownames(i_area_pd) <-  NULL
colnames(i_area_pd) <- c("PD", "SR", "region", "area")
i_area_pd$region <- factor(i_area_pd$region, levels = c("Sunda", "Island", "Sahul"))
i_area_pd$area <- factor(i_area_pd$area, levels = c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA"))

j_area_pd <- pd(j_sp_mx, tr_j)
j_area_pd <- cbind(j_area_pd, region$region)
j_area_pd$area <- rownames(j_area_pd)
rownames(j_area_pd) <-  NULL
colnames(j_area_pd) <- c("PD", "SR", "region", "area")
j_area_pd$region <- factor(j_area_pd$region, levels = c("Sunda", "Island", "Sahul"))
j_area_pd$area <- factor(j_area_pd$area, levels = c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA"))

area_pd <- pd(sp_mx, tr)
area_pd <- cbind(area_pd, region$region)
area_pd$area <- rownames(area_pd)
rownames(area_pd) <-  NULL
colnames(area_pd) <- c("PD", "SR", "region", "area")
area_pd$region <- factor(area_pd$region, levels = c("Sunda", "Island", "Sahul"))
area_pd$area <- factor(area_pd$area, levels = c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA"))

## raw area PD PE scores ####
canape_sum_obs <- data.frame(rpd_obs = canape_res$rpd_obs, rpe_obs = canape_res$rpe_obs, 
                         rpd_obs_j = canape_res_j$rpd_obs, rpe_obs_j = canape_res_j$rpe_obs,
                         rpd_obs_i = canape_res_i$rpd_obs, rpe_obs_i = canape_res_i$rpe_obs,
                         area = rownames(canape_res), region = region$region,
                         pd_obs = canape_res$pd_obs, pe_obs = canape_res$pe_obs, 
                         pd_obs_j = canape_res_j$pd_obs, pe_obs_j = canape_res_j$pe_obs,
                         pd_obs_i = canape_res_i$pd_obs, pe_obs_i = canape_res_i$pe_obs)

canape_sum_obs$area <- factor(canape_sum_obs$area, 
                                levels = c("SEA", "SUM", "BOR", "JAV",
                                           "LSU", "SUL", "MAU", "EMN", "NG", "AUS"))

canape_long <- gather(canape_sum_obs, key = "variable", value = "value", -area, -region)
canape_long$type <-  ifelse(grepl("^pd_obs", canape_long$variable), "pd_obs",
                            ifelse(grepl("^pe_obs", canape_long$variable), "pe_obs",
                                   ifelse(grepl("^rpd_obs", canape_long$variable), "rpd_obs",
                                          ifelse(grepl("^rpe_obs", canape_long$variable), "rpe_obs", "other"))))
canape_long$tree <- ifelse(grepl("^pd_obs_i$|^pe_obs_i$|^rpd_obs_i$|^rpe_obs_i$", canape_long$variable), "Jetz (imputed)",
                           ifelse(grepl("^pd_obs_j$|^pe_obs_j$|^rpd_obs_j$|^rpe_obs_j$", canape_long$variable), "Jetz (genetic)", "New tree"))
canape_long$area <- factor(canape_long$area, 
                              levels = c("SEA", "SUM", "BOR", "JAV",
                                         "LSU", "SUL", "MAU", "EMN", "NG", "AUS"))

newtree <- 
ggplot(data = subset(canape_long, tree == "New tree"), aes(x = area, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~type, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(pd_obs = "PD", pe_obs = "PE", rpd_obs = "rpd", rpe_obs = "rpe")), strip.position = "left") +
  labs(x = "", y = "" ) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_fill_discrete(labels = c("New Tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),  # Add black line to axes
        strip.text = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

jetzg <- 
  ggplot(data = subset(canape_long, tree == "Jetz (genetic)"), aes(x = area, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~type, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(pd_obs = "PD", pe_obs = "PE", rpd_obs = "rpd", rpe_obs = "rpe")), strip.position = "left") +
  labs(x = "", y = "" ) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_fill_discrete(labels = c("New Tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),  # Add black line to axes
        strip.text = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

jetzi <- 
  ggplot(data = subset(canape_long, tree == "Jetz (imputed)"), aes(x = area, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~type, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(pd_obs = "PD", pe_obs = "PE", rpd_obs = "rpd", rpe_obs = "rpe")), strip.position = "left") +
  labs(x = "", y = "" ) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_fill_discrete(labels = c("New Tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),  # Add black line to axes
        strip.text = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

grid.arrange(jetzi, jetzg, newtree, ncol = 3)

pdp <- ggplot(data = subset(canape_long, type == "pd_obs"), aes(x = area, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "PD") +
  scale_y_continuous(breaks = seq(0, 0.5, by =0.1), limits = c(0, 0.5), expand = c(0,0)) +
  scale_fill_manual(values = c("pd_obs_i" = "#F8766D", "pd_obs_j" = "#00BA38", "pd_obs" = "#619CFF"),
                    labels = c("New tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(legend.position = "none", legend.text = element_text(size = 12),
        axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

pep <- ggplot(data = subset(canape_long, type == "pe_obs"), aes(x = area, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "PE") +
  scale_y_continuous(breaks = seq(0, 0.3, by =0.1), limits = c(0, 0.3), expand = c(0,0)) +
  scale_fill_manual(values = c("pe_obs_i" = "#F8766D", "pe_obs_j" = "#00BA38", "pe_obs" = "#619CFF"),
                    labels = c("New tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(legend.position = "none", axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

rpdp <- ggplot(data = subset(canape_long, type == "rpd_obs"), aes(x = area, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Relative PD") +
  scale_y_continuous(breaks = seq(0, 1.3, by =0.25), limits = c(0, 1.3), expand = c(0,0)) +
  scale_fill_manual(values = c("rpd_obs_i" = "#F8766D", "rpd_obs_j" = "#00BA38", "rpd_obs" = "#619CFF"),
                    labels = c("New tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(legend.position = "none", axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

rpep <- ggplot(data = subset(canape_long, type == "rpe_obs"), aes(x = area, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Relative PE") +
  scale_y_continuous(breaks = seq(0, 1.3, by =0.25), limits = c(0, 1.3), expand = c(0,0)) +
  scale_fill_manual(values = c("rpe_obs_i" = "#F8766D", "rpe_obs_j" = "#00BA38", "rpe_obs" = "#619CFF"),
                    labels = c("New tree", "Jetz (imputed)", "Jetz (genetic)"), name = "") +
  theme_minimal() +
  theme(legend.position = "none", axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

grid.arrange(pdp, pep, rpdp, rpep, ncol = 2)

pde.pd <- 
  ggplot(canape_sum_obs, aes(x = rev(area), fill = rev(region))) +
  geom_bar(aes(y=pd_obs), position = "dodge", stat = "identity") +
  labs(x = "Area", y = "Phylogenetic Diversity") +
  scale_y_continuous(breaks = seq(0, 0.5, by =0.1), limits = c(0, 0.5), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

pde.pe <- 
  ggplot(canape_sum_obs, aes(x = rev(area), fill = rev(region))) +
  geom_bar(aes(y=pe_obs), position = "dodge", stat = "identity") +
  labs(x = "Area", y = "Phylogenetic Endemism") +
  scale_y_continuous(breaks = seq(0, 0.3, by =0.1), limits = c(0, 0.3), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

pde.rpd <- 
  ggplot(canape_sum_obs, aes(x = rev(area), fill = rev(region))) +
  geom_bar(aes(y=rpd_obs), position = "dodge", stat = "identity") +
  labs(x = "Area", y = "Relative Phylogenetic Diversity") +
  scale_y_continuous(breaks = seq(0, 1.2, by =0.25), limits = c(0, 1.1), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

pde.rpe <- 
  ggplot(canape_sum_obs, aes(x = rev(area), fill = rev(region))) +
  geom_bar(aes(y=rpe_obs), position = "dodge", stat = "identity") +
  labs(x = "Area", y = "Relative Phylogenetic Endemism") +
  scale_y_continuous(breaks = seq(0, 1.3, by =0.25), limits = c(0, 1.3), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        axis.line = element_line(color = "black"),  # Add black line to axes
        #panel.grid.major = element_blank(),  # Remove major grid lines
        axis.text = element_text(size = 12),  # Adjust text size
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

grid.arrange(pde.pd, pde.pe, pde.rpd, pde.rpe, ncol =2)

# Area Coverage ####
## Between area SR comparison ####
cf <- data.frame(ours = area_pd$SR, jetz = j_area_pd$SR)
cf$area <- as.factor(area_pd$area)
cf$iaa <- c(346, 382, 198, 158, 140, 180, 265, 197, 268, 465)
cf$diff <- round(cf$jetz/cf$iaa, 3)
cf$diff2 <- round(cf$ours/cf$iaa, 3)
head(cf)

ggplot(cf, aes(x = area)) +
  geom_bar(aes(y=iaa, fill = "Unsampled"), position = "dodge", stat = "identity") +
  geom_bar(aes(y=ours, fill = "New tree"), position = "dodge", stat = "identity") +
  geom_bar(aes(y=jetz, fill = "Jetz (genetic)"), position = "dodge", stat = "identity") +
  geom_text(aes(y=jetz-10, label = paste0(scales::percent(diff)))) +
  geom_text(aes(y=ours-10, label = paste0(" + ", scales::percent((diff2-diff))))) +
  geom_text(aes(y=iaa-8, label = paste0(" + ", scales::percent((1-diff2))))) +
  labs(y = "No. of species", x = "Area") +
 # scale_fill_discrete(labels = c("Unsampled", "Jetz (genetic)", "New tree"), name = "") +
  scale_fill_manual(values = c("Unsampled" = "grey", "Jetz (genetic)" = "#00BA38", "New tree" = "#619CFF"),
                    labels = c("Jetz (genetic)", "New tree", "Unsampled"), name = "") +
  scale_y_continuous(breaks = seq(0, 450, by = 100), limits = c(0, 500), expand = c(0,0)) +
  theme_minimal() +
  theme(  # Remove legend
    axis.line = element_line(color = "black"),  # Add black line to axes
    axis.text = element_text(size = 12),  # Adjust text size
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)), 
    legend.position = "top",
    legend.text = element_text(size = 12))

grid.arrange(tax.coverage, area.coverage, ncol = 2)

# Branch lengths of terminal tips ########
term_tips_i <- tr_i$edge.length[tr_i$edge[,2] <= Ntip(tr_i)]
term_tips_i <- as.data.frame(term_tips_i)
term_tips_a <- 
  ggplot(term_tips_i, aes(x = term_tips_i)) +
  geom_density() +
  labs(title = "a) Jetz (imputed)", x = "Length", y = "Density")+
  xlim(0, 90) +
  ylim(0, 0.15) +
  theme_minimal() +
  theme(  # Remove legend
    axis.line = element_line(color = "black"),  # Add black line to axes
    axis.text = element_text(size = 12),  # Adjust text size
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)))

term_tips_j <- tr_j$edge.length[tr_j$edge[,2] <= Ntip(tr_j)]
term_tips_j <- as.data.frame(term_tips_j)
term_tips_b <- 
  ggplot(term_tips_j, aes(x = term_tips_j)) +
  geom_density() +
  labs(title = "b) Jetz (genetic)", x = "Length", y = "Density")+
  xlim(0, 90) +
  ylim(0, 0.15) +
  theme_minimal() +
  theme(  # Remove legend
    axis.line = element_line(color = "black"),  # Add black line to axes
    axis.text = element_text(size = 12),  # Adjust text size
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)))

term_tips <- tr$edge.length[tr$edge[,2] <= Ntip(tr)]
term_tips <- as.data.frame(term_tips)
term_tips_c <- 
  ggplot(term_tips, aes(x = term_tips)) +
  geom_density() +
  labs(title = "c) New tree", x = "Length", y = "Density")+
  xlim(0, 90) +
  ylim(0, 0.15) +
  theme_minimal() +
  theme(  # Remove legend
    axis.line = element_line(color = "black"),  # Add black line to axes
    axis.text = element_text(size = 12),  # Adjust text size
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)))

grid.arrange(term_tips_a, term_tips_b, term_tips_c, ncol =3)


