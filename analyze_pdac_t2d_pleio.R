# Analyze pdac-t2d pleiotropy data.
# Running this code with the associated data
# (provided in https://odin.mdacc.tmc.edu/~rsun3/pdac_t2d_reproduce/)
# will reproduce results in our paper.

library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(SeqArray)
library(pROC)


# Function to identify different loci from raw hits
find_loci <- function(rejectionDat, threshold=100000) {
  
  # order by chromosome and BP
  rejectionDat <- rejectionDat %>% arrange(Chr, pos38)
  uniqueChrs <- sort(unique(rejectionDat$Chr))
  
  # loop through chromosomes
  res <- c()
  endPosRes <- c()
  lfdrRes <- c()
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- rejectionDat %>% filter(Chr == tempChr)
    
    # hold the results 
    res <- rbind(res, tempDat[1, ]) 
    # start of the BP window
    tempPos <- tempDat$pos38[1]
    # significance
    minLfdr <- tempDat$rawlfdr[1]
    if (nrow(tempDat) == 1) {
      endPosRes <- c(endPosRes, tempDat$pos38[1])
      lfdrRes <- c(lfdrRes, minLfdr)
      next 
    }
    
    # loop from smallest BP to highest BP on the chromosomes
    for (snp_it in 2:nrow(tempDat)) {
      
      # if BP jump is larger than threshold, it's a new locus
      if (tempDat$pos38[snp_it] > tempPos + threshold) {
        
        # add to res 
        res <- rbind(res, tempDat[snp_it, ])
        # update min lfdr
        lfdrRes <- c(lfdrRes, minLfdr)
        minLfdr <- tempDat$rawlfdr[snp_it]
        # update end of BP window
        endPosRes <- c(endPosRes, tempDat$pos38[snp_it - 1])
      }
      
      # update min lfdr 
      if (tempDat$rawlfdr[snp_it] < minLfdr) {minLfdr <- tempDat$rawlfdr[snp_it]}
      # update tempPos
      tempPos <- tempDat$pos38[snp_it]
      
      # end of chr
      if (snp_it == nrow(tempDat)) {
        # update window and lfdr  
        lfdrRes <- c(lfdrRes, minLfdr)
        endPosRes <- c(endPosRes, tempDat$pos38[snp_it]) 
      }
    }
  }
  
  # append lfdr and endpos
  res <- res %>% mutate(minLfdr = lfdrRes, endPos = endPosRes)
  return(res)
}


# load data
mvpHuerta <- fread('/users/rsun3/desktop/pdac_pleio/final_pdac_mvp_diabetes_huerta.txt') %>%
  mutate(lfdr = unlist(fread("/users/rsun3/desktop/pdac_pleio/pdac_mvp_diabetes_huerta_newlfdr.txt"))) %>%
  arrange(lfdr) %>%
  mutate(cumlfdr = cummean(lfdr)) %>%
  filter(abs(Zpdac) > 1.96 | abs(Zdiab) > 1.96)

# smaller data, rejections
huertaSmall <- mvpHuerta %>% filter(abs(Zpdac) > 2.807 | abs(Zdiab) > 2.807)
huertaRej <- mvpHuerta %>% filter(cumlfdr < 0.1)
colnames(huertaRej)[c(2, 3, 21)] <- c("Chr", "pos38", "rawlfdr")
topLoci <- find_loci(huertaRej)

#----------------------------------------------------------------------------------------------------#
# make Table 1

# pick out the genes  around the top loci
# also pick out the top snp at each loci
geneInfo <- fread("/users/rsun3/desktop/pdac_pleio/ensembl_refgene_hg19_20180109_withgrch38.txt")
geneTab <- data.frame(Chr = topLoci$Chr, Start = topLoci$pos38, End = topLoci$endPos,
                      G1=NA, G2=NA)
snpTab <- c()
buffer <- 25000
for (row_it in 1:nrow(geneTab)) {
  
  # find the gene
  tempChr <- geneTab$Chr[row_it]
  tempStart <- geneTab$Start[row_it]
  tempEnd <- geneTab$End[row_it]
  tempDat <- geneInfo %>% filter(Chr == tempChr)
  keepRows1 <- which(tempDat$start38 - buffer <= tempStart & tempDat$end38 >= tempStart)
  keepRows2 <- which(tempDat$start38 <= tempEnd & tempDat$end38 + buffer >= tempEnd)
  keepRows3 <- which(tempDat$start38 >= tempStart & tempDat$end38 <= tempEnd)
  keepRows <- unique(c(keepRows3, keepRows2, keepRows1))
  keepDat <- tempDat[keepRows, ]
  if (length(keepRows) > 0) {
    for (col_it in 1:length(keepRows)) {
      geneTab[row_it, col_it + 3] <- keepDat$HGNC_name[col_it]
    }
  }
  
  # find the top snp
  tempSnps <- huertaRej %>% filter(Chr == tempChr & pos38 >= tempStart & pos38 <= tempEnd) %>%
    arrange(rawlfdr) %>%
    select(Chr, pos38, rsid, effect_allele, other_allele, effect_allele_frequency, 
           pval_pdac, p_value, rawlfdr, cumlfdr) 
  topSnp <- tempSnps %>% mutate(num = nrow(tempSnps)) %>% 
    slice(1)
  snpTab <- rbind(snpTab, topSnp)
}

# make table for paper
tab1 <- cbind(snpTab, geneTab %>% select(-Chr)) %>% mutate(MAF = 1 - effect_allele_frequency) %>%
  select(Chr, Start, End, num, rsid, pos38, MAF, pval_pdac, p_value) %>%
  set_colnames(c("Chr", "Start", "End", "Num", "RS", "Pos", "MAF", "Pval_PC", "Pval_T2D")) 

#------------------------------------------------------------------------------------------------------------#
# Make Figure 1


# make manhattan plot - get positions
chrSize <- fread("/users/rsun3/desktop/pdac_pleio/chr_max_pos38.txt")
cumDF <- data.frame(Chr = 1:22, Size = unlist(chrSize$maxpos), 
                    startPos = c(0, cumsum(as.numeric(unlist(chrSize$maxpos)))[1:21]))

# add total positions
huertaSmall <- huertaSmall %>% mutate(totalPos = cumDF$startPos[chromosome.x] + base_pair_location.x) %>%
  mutate(Even = ifelse(chromosome.x%%2 == 0, 1, 0)) %>%
  mutate(pleio = ifelse(cumlfdr < 0.1, 1, 0)) %>%
  mutate(pval_diab = as.numeric(p_value))


# put it in plot form
huertaPlot <- huertaSmall %>% pivot_longer(cols = c(pval_pdac, pval_diab), 
                                           names_to="Pheno", 
                                           values_to="pvalue") %>%
  mutate(Pheno = ifelse(Pheno == "pval_pdac", "pdac", "diabetes")) %>%
  select(totalPos, Pheno, pvalue, Even, pleio, lfdr) %>%
  mutate(shapeFactor = ifelse(pleio == 1, "Pleiotropic", "Not")) %>% 
  mutate(colFactor = ifelse(Even == 1, "even", "odd")) %>%
  mutate(colFactor = ifelse(pleio == 1, "Pleiotropic", colFactor))

# diabetes manhattan
diabMan <- ggplot(huertaPlot %>% filter(Pheno == "diabetes"), aes(x=totalPos, y = -log10(pvalue))) + 
  geom_point(data = filter(huertaPlot, Pheno == "diabetes" & colFactor != "Pleiotropic"), 
             aes(color = as.factor(Even), shape=shapeFactor)) + 
  geom_point(data = filter(huertaPlot, Pheno == "diabetes" & colFactor == "Pleiotropic"), 
             aes(color = "Pleio", shape=shapeFactor)) + 
  geom_hline(yintercept = -log10(5e-8), color = "black", linetype = "dashed", linewidth = 0.5) +
  ylim(c(0, 320)) +
  scale_shape_manual(values=c(16, 17), labels=c("Non-Pleiotropic", "Pleiotropic"), name="Type") + 
  scale_color_manual(name = "Types", values=
                       c("1" = "grey", "0" = "red", "Pleio" = "blue"), 
                     labels=c("Non-Pleiotropic (Odd Chr)", "Non-Pleiotropic (Even Chr)", "Pleiotropic")) + 
  scale_x_continuous(
    breaks = cumsum(as.numeric(unlist(chrSize$maxpos))), # Explicitly set tick marks for every integer 1 to 22
    labels=c(1:22), # Ensures all ticks are visible
    expand = c(0.01, 0)
  ) +
  labs(title = "", x="Position", y = expression("-log"[10] * "(P-value)")) + 
  guides(color = guide_legend(override.aes = list(
    shape = c(16, 16, 17))),
    shape = "none"
  ) + 
  theme_cowplot()

# pdac manhattan
pdacMan <- ggplot(huertaPlot %>% filter(Pheno == "pdac"), aes(x=totalPos, y = -log10(pvalue))) + 
  geom_point(data = filter(huertaPlot, Pheno == "pdac" & colFactor != "Pleiotropic"), 
             aes(color = as.factor(Even), shape=shapeFactor)) + 
  geom_point(data = filter(huertaPlot, Pheno == "pdac" & colFactor == "Pleiotropic"), 
             aes(color = "Pleio", shape=shapeFactor)) + 
  geom_hline(yintercept = -log10(5e-8), color = "black", linetype = "dashed", linewidth = 0.5) +
  ylim(c(0, 12.5)) +
  scale_shape_manual(values=c(16, 17), labels=c("Non-Pleiotropic", "Pleiotropic"), name="Type") + 
  scale_color_manual(name = "Types", values=
                       c("1" = "grey", "0" = "red", "Pleio" = "blue"), 
                     labels=c("Non-Pleiotropic (Odd Chr)", "Non-Pleiotropic (Even Chr)", "Pleiotropic")) + 
  scale_x_continuous(
    breaks = cumsum(as.numeric(unlist(chrSize$maxpos))), # Explicitly set tick marks for every integer 1 to 22
    labels=c(1:22), # Ensures all ticks are visible
    expand = c(0.01, 0)
  ) +
  labs(title = "", x="Top: Pancreatic Cancer, Bottom: Type 2 Diabetes", y = expression("-log"[10] * "(P-value)")) + 
  guides(color = guide_legend(override.aes = list(
    shape = c(16, 16, 17))),
    shape = "none"
  ) + 
  theme_cowplot()

# reverse diabetes and make other changes
diabMan2 <- diabMan + scale_y_reverse(limits=c(300, 0)) + 
  labs(x="") + 
  xlim(c(0, 2871187872)) + 
  scale_x_continuous(
    breaks = cumsum(as.numeric(unlist(chrSize$maxpos))), # Explicitly set tick marks for every integer 1 to 22
    labels=c(1:18, "", 20, "", 22), # Ensures all ticks are visible
    expand = c(0.01, 0),
    position = "top"
  ) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.justification="center") + 
  theme(plot.margin = unit(c(t=-10, r=5.5, b=5.5, l=2.5), "pt"))

pdacMan2 <- pdacMan + xlim(c(0, 2871187872)) + 
  scale_x_continuous(
    breaks = cumsum(as.numeric(unlist(chrSize$maxpos))), # Explicitly set tick marks for every integer 1 to 22
    labels=c(1:18, "", 20, "", 22), # Ensures all ticks are visible
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = c(3, 6, 9, 12), # Explicitly set tick marks for every integer 1 to 22
    labels=c(3, 6, 9, 12), # Ensures all ticks are visible
    expand = c(0.01, 0),
    limits=c(0, 12.5)
  ) + 
  theme(legend.position ="none")  + 
  theme(plot.margin = unit(c(t=5.5, r=5.5, b=-10, l=9.5), "pt"))


# put together
fullMan <- grid.arrange(arrangeGrob(pdacMan2), 
                        arrangeGrob(diabMan2), padding=0, heights=c(0.5, 0.5))
fullMan


#--------------------------------------------------------------------------------------------#
# Make Figure 2

# downsample the test statistics so that our computer can handle plotting it
mvpHuertaDat <- fread('/users/rsun3/desktop/pdac_pleio/final_pdac_mvp_diabetes_huerta.txt') %>%
  mutate(lfdr = unlist(fread("/users/rsun3/desktop/pdac_pleio/pdac_mvp_diabetes_huerta_newlfdr.txt"))) %>%
  arrange(lfdr) %>%
  mutate(cumlfdr = cummean(lfdr)) 
huertaZdat1 <- mvpHuertaDat %>% filter(abs(Zpdac) > 2.807 | abs(Zdiab) > 2.807)
huertaZdat2 <- mvpHuertaDat %>% filter(abs(Zpdac) < 2.807 & abs(Zdiab) < 2.807)
huertaZdat2 <- huertaZdat2 %>%
  dplyr::slice_sample(., prop=0.01) 
huertaZdat <- rbind(huertaZdat1, huertaZdat2) %>% 
  mutate(pleio = ifelse(cumlfdr < 0.1, 1, 0))

# check how many have opposite signs
huertaZrej <- huertaZdat %>% filter(pleio == 1)
length(which(sign(huertaZrej$Zdiab) == sign(huertaZrej$Zpdac)))

huertaZ <- ggplot(data=huertaZdat, aes(x=Zpdac, y = Zdiab)) + 
  geom_point(data = filter(huertaZdat, pleio == 0), 
             aes(color = as.factor(pleio), shape=as.factor(pleio))) + 
  geom_point(data = filter(huertaZdat, pleio == 1), 
             aes(color = as.factor(pleio), shape=as.factor(pleio))) + 
  scale_shape_manual(values=c(16, 17), labels=c("Non-Pleiotropic", "Pleiotropic"), name="Type") + 
  scale_color_manual(name = "Type", values=
                       c("1" = "blue", "0" = "red"), 
                     labels=c("Non-Pleiotropic", "Pleiotropic")) + 
  scale_x_continuous(
    breaks = seq(from=-8, to=8, by=2)
  ) + 
  scale_y_continuous(
    breaks = seq(from=-40, to=40, by=5)
  ) + 
  labs(x="Z Statistics Type II Diabetes", y="Z Statistics Pancreatic Cancer") +
  theme_cowplot()
huertaZ

#-----------------------------------------------------------------------------------------------#
# Make Figure 3

# read FAVOR annotation data
favorDat <- fread("/users/rsun3/desktop/pdac_pleio/2025-10-06--20-43-31_bothFAVOR_processed.csv")   %>%
  select(VariantVcf.String, Chromosome.String, Position.String, CaddPhred.Float64,
         ApcConservationV2.Float64, ApcConservationV2.Valid, ApcEpigeneticsActive.Float64, 
         ApcEpigeneticsActive.Valid) %>%
  mutate(name = paste0(Chromosome.String, ":", Position.String))
huertaFavorDat <- merge(favorDat, huertaSmall %>% mutate(name = paste0(chromosome.x, ":", base_pair_location.x)) %>%
                          select(name, cumlfdr, lfdr), by="name") %>%
  filter(cumlfdr < 0.1) %>%
  distinct() %>%
  arrange(cumlfdr)
dim(huertaFavorDat)
dim(huertaFAVOR)

# plot the favor epigenetic on chr 16
favor16epi <- ggplot(huertaFavorDat %>% filter(Chromosome.String == 16 & as.numeric(Position.String) > 75000000), 
                     aes(x=as.numeric(Position.String), y = as.numeric(ApcEpigeneticsActive.Float64))) + 
  geom_point(data =huertaFavorDat %>% filter(Chromosome.String == 16 & as.numeric(Position.String) > 75000000), 
             aes(color = -log10(lfdr))) + 
  scale_color_gradient(low="red", high="blue", name=expression("-log"[10] * "(lfdr-value)"))  +
  scale_x_continuous(
    breaks = c(75200000, 75300000, 75400000),
    labels= c(75200000, 75300000, 75400000),
    expand = c(0.01, 0)
  ) +
  labs(title = "", x="Position (Chr 16)", y = "FAVOR Epigenetic") + 
  theme_cowplot() 

# plot the favor conservation on chr 16
favor16cons <- ggplot(huertaFavorDat %>% filter(Chromosome.String == 16 & as.numeric(Position.String) > 75000000), 
                      aes(x=as.numeric(Position.String), y = as.numeric(ApcConservationV2.Float64))) + 
  geom_point(data =huertaFavorDat %>% filter(Chromosome.String == 16 & as.numeric(Position.String) > 75000000), 
             aes(color = -log10(lfdr))) + 
  scale_color_gradient(low="red", high="blue", name=expression("-log"[10] * "(lfdr-value)")) +
  scale_x_continuous(
    breaks = c(75200000, 75300000, 75400000),
    labels= c(75200000, 75300000, 75400000),
    expand = c(0.01, 0)
  ) +
  labs(title = "", x="Position (Chr 16)", y = "FAVOR Conservation") + 
  theme_cowplot() 

# put it together
favor_plot <- plot_grid(favor16epi + theme(legend.position = "none"),
                        favor16cons + theme(legend.position = "none"),
                        labels=c("A", "B"), nrow=2, label_size=22)
favor_legend <- get_legend(favor16epi +  theme(legend.direction="horizontal",
                                               legend.justification="center",
                                               legend.box.just="bottom"))
favor_plot_full <- plot_grid(favor_plot, favor_legend, ncol=1, rel_heights=c(1, 0.1))
favor_plot_full


#--------------------------------------------------------------------------#
# Figures 4 + 5 require individual-level data from the UK Biobank.
# Unfortunately we cannot share this data, and so we cannot demonstrate
# how to reproduce those figures. 
# Qualified researchers can apply to access the data.
# If you have access to the UK Biobank, feel free to email our group,
# and we will share how to reproduce Figures 4 and 5.



