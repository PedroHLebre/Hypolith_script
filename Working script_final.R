#packages to load

library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)


#import tables. The ASVs, taxonomy and tree files used are the output from the QIIME2 pipeline, with taxonomic annotation using the silvav1.32 classifier. The ASV table has been filtered to exclude reads with less than 20 counts prior to this step: 
asv_table=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
taxonomy = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",", quote = "")
metadata = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")

#conver tables into phyloseq objects 
asv_table_p = as.matrix(asv_table)
taxonomy_p = as.matrix (taxonomy)
OTU_tree = compute.brlen(tree, method = "Grafen")
AfSM_tree = read_tree(treefile = "./tree.nwk", errorIfNULL = FALSE)
AfSM_asv = otu_table (asv_table_p, taxa_are_rows = TRUE)
AfSM_tax = tax_table(taxonomy_p)
AfSM_sampledata = sample_data(metadata)

#creation of the phyloseq object
physeq = phyloseq(AfSM_asv, AfSM_sampledata, AfSM_tax, AfSM_tree)

#get rid of misannotations
physeq_filtered = subset_taxa( physeq,tax_table(physeq) != "Eukaryota")

#rarefaction curve to access sample coverage
rarecurve(t(otu_table(physeq_filtered)), step=50, cex=0.5)

#Filter out sample with low coverage 
physeq__trimmed = subset_samples ( physeq_filtered, sample_names(physeq_filtered) != "H17")

#group taxa into phyla
physeq_Phylum= tax_glom (physeq__trimmed, taxrank = "Phylum")

#Pipeline to create relative abundance plot. The first step is to transform the sample counts into fraction of total counts per sample:
physeq_Phylum_rel = transform_sample_counts(physeq_Phylum, function (x) x/sum(x) )

#obtain only the top phyla (representing 99% of reads)
physeq_dominant= filter_taxa(physeq_Phylum_rel, function(x) mean(x) > 0.01, TRUE)


# make the relative abundance barplot
RA_bact = plot_bar(physeq_dominant, fill = "Phylum") + labs(x="Sample", y = "Abundance")  + geom_bar(stat = "identity", position = "fill", colour = "black") 
# change the colours of the different phyla
RA_bact_3 = RA_bact + scale_fill_manual(values = c("Cyanobacteria" = "green", "Acidobacteria" = "orange", "Actinobacteria" = "blue", "Bacteroidetes" = "gold", "Chloroflexi" = "darkgreen", "Deinococcus-Thermus" = "violetred1", "Gemmatimonadetes" = "dark cyan", "Planctomycetes" = "yellow1", "Proteobacteria" = "red", "Verrucomicrobia" = "plum4"))
#Re-order the samples so that they are in numerical order
RA_bact_3$data$Sample = factor(RA_bact_3$data$Sample, levels = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H14", "H15", "H16", "H18", "H19", "H20","H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30"))
#save it as an image
ggsave(filename = "Relative_abundance.png", RA_bact_3,width = 8,height = 6, dpi = 900)

#export the ASV and taxonomy table for phyla taxrank. These two tables will be merged together with the metadata for the grouping from the clustering analysis. This merged table was then transposed, resulting in a table with samples as rows, and metadata, including phyla relative abundances, as columns:  
write.csv(otu_table(physeq_dominant), "ASV_top_phyla.csv")
write.csv(tax_table(physeq_dominant), "tax_top_phyla.csv")

#The merged table is then imported into RStudio: 
merged_phyla_table = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

#checking for the normality of relative abundance distribution. This was performed for each phyla: 
shapiro.test(merged_phyla_table$Actinobacteria)

#If the relative abundance distribution is normal, the impact of the groupings on the relative abundance is tested using ANOVA: 
aov = aov(Actinobacteria ~ Group, data = merged_phyla_table)
summary(aov)
#Alternately, a non-nomal distribution is tested using Kruskal-Wallis: 
krustal.test(Actinobacteria ~ Group, data = merged_phyla_table)

#boxplot of significant 
Phyla_compare = ggplot (merged_phyla_table, aes( x= Group, y = Actinobacteria, fill = Group)) + geom_boxplot() + geom_point (aes(fill = Group), size = 3, position = position_jitterdodge())+theme_classic()
Phyla_compare2 = observed_compare2 + theme(axis.text.x = element_text(size=10, angle=45), axis.text.y = element_text(face="bold",size = 14))


#Pipeline for the alpha and beta-diversity metrics
#In our methodology, we rarefy the ASv table for beta-diversity using the sample with the lowest read count: 
physeq_ra = rarefy_even_depth(physeq__trimmed, rngseed = 1, sample.size = 1*min(sample_sums(physeq__trimmed)), replace = F)
#Beta-diversity metrics assume that the data is normally distributed. We normalize the phyloseq rarefied object with log scale
physeq_transformed = transform_sample_counts(physeq_ra, function (x) log(x+1))

#Clustering analysis to determine the sample groupings if any. The clustering was done using weighted unifrac: 
d = distance(physeq_transformed, method="wunifrac", type="samples")
hell.hclust     = hclust(d, method="ward.D2")
plot(hell.hclust)


#Calculation of beta-diversity and significance of dissimilarity between groups determiined by the clustering analysis
erie_PCoA = phyloseq::distance(physeq_transformed, method = "wunifrac")
sampledf = data.frame(sample_data(physeq_transformed))
adonis(erie_PCoA ~ Group, data = sampledf, permutations = 1000)

#calculation of beta-dispersivity to determine if the variablilty within groups is significant
beta = betadisper(erie_PCoA, sampledf$Group)
permutest(beta)

#plot the weighte unifrac distance on a PCoA plot 
erie_PCoA = ordinate(physeq = physeq_transformed, method = "PCoA", distance = "Unifrac", weighted = TRUE)
p1 = plot_ordination(physeq = physeq_transformed,ordination = erie_PCoA,color = "Group" , axes = 1:2, label = "Description") + theme_bw() 
p2 = p1 +  geom_point(size = 3, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#save the PCoA plot as a image
ggsave(filename = "PCoA_WUnifrac_final.png", p3 , width = 8, height = 6, dpi = 900)


#Zeta diversity pipeline, to calculate the decay in community similarity as a function of distance

library(zetadiv)

#import table with transposed ASVs and location coordinates (for latitude and longitude in metres) in the first two columns 
table_b = read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")
#make a table  just with the distance 
xy.bact = table_b[,1:2]
#specify another table with just the ASvs 
data.spec.bact = table_b[,3:3330]
#calculate zeta decay for a given order, and plot it. The order signifies the number of samples are compared. SO order 2 means the zeta diversity decays in comparisons between any two samples over distance.
zeta.ddecay.bact = Zeta.ddecay(xy.bact, data.spec.bact, sam = 200, order = 2 ,method.glm = "glm.fit2", confint.level = 0.95, plot = FALSE, normalize = "Sorensen")
zeta_plot = Plot.zeta.ddecay(zeta.ddecay.bact)

#export the plot as an image
ggsave(filename = "Zeta_cuttoff20reads.png", zeta_plot , width = 8, height = 6, dpi = 900)


#DESeq2 was used to find significant differences in abundance between taxa at the genus level accross the groups; \

library(DESeq2)

#cluster the ASV table in the phyloseq object into genus taxrank
physeq_genus = tax_glom (physeq__trimmed, taxrank = "Genus")
#Conversion of the phyloseq into a DEseq object using the groups as the experimental design
diagdds = phyloseq_to_deseq2(physeq_genus, ~ Group )
#Transform deseq object. This is required because the ASV table contains many zeros, which will affect the analysis
dds = estimateSizeFactors(diagdds, type= c("poscount")) # The poscount transformation is adapted specifically for metagenomic data. 
#Perform Deseq analYsis
diagdds_deseq = DESeq(dds, test="Wald", fitType="local")
# Create a table containing pair-wise comparison between 2 groups specified in the script. These step was done for all group combinations
res= results(diagdds_deseq, contrast=c("Group","A","B"))
#filter the results to include only significant p-adj values 
alpha = 0.01 # this is your p-value treshold
res_filtered = res[which(res$padj < alpha), ]
# add taxonomy information to the results
sigtab = cbind(as(res_filtered, "data.frame"), as(tax_table(physeq_genus)[rownames(res_filtered), ], "matrix"))

# export the table as csv 
write.csv(sigtab, "Deseq_results_genus_A_B.csv") 
#At this point, the tables for each pair-wise combination are concatenated in excel in order to identify the taxa that are over-represented in on group vs all others. Excel was also used to draw the barplot in figure 4. 

