library(MoonlightR)

# Load data 
DEGs <- read.csv("../results/biclass_trueDEGs_data.csv")
cluster1_DEGs <- read.csv("../results/cola_cluster1.csv")
cluster2_DEGs <- read.csv("../results/cola_cluster2.csv")
BRCA_DEGs <- get(load("../../Omics_integration/data/BRCA_dataDEGs.rda"))
rm(dataDEGs)


# Extract DEGs
DEG_expr <- BRCA_DEGs[DEGs$Gene.Symbol,]
cluster1_DEG_expr <- BRCA_DEGs[cluster1_DEGs$x,]
cluster2_DEG_expr <- BRCA_DEGs[cluster2_DEGs$x,]


# Functional Enrichment Analysis
dataFEA <- FEA(DEGsmatrix = DEG_expr)
save(dataFEA, file = "../results/trueDEGs_FEA.rda")

cluster1_FEA <- FEA(DEGsmatrix = cluster1_DEG_expr)
save(cluster1_FEA, file = "../results/cluster1_FEA.rda")

cluster2_FEA <- FEA(DEGsmatrix = cluster2_DEG_expr)
save(cluster2_FEA, file = "../results/cluster2_FEA.rda")

# Plotting
# Find significant BPs
index <- which(abs(dataFEA$Moonlight.Z.score) >= 1 & dataFEA$FDR <= 0.01)
index1 <- which(abs(cluster1_FEA$Moonlight.Z.score) >= 1 & cluster1_FEA$FDR <= 0.01)
index2 <- which(abs(cluster2_FEA$Moonlight.Z.score) >= 1 & cluster2_FEA$FDR <= 0.01)

# Plot
plotFEA(dataFEA = dataFEA[index,],
        topBP = nrow(dataFEA[index,]), 
        additionalFilename = "../figures/trueDEGs_FEAplot.pdf", 
        # height = 10,
        width = 11, 
        offsetValue = 5, angle = 90, 
        xleg = 0, yleg = 0)

plotFEA(dataFEA = cluster1_FEA[index1,],
        topBP = nrow(cluster1_FEA[index1,]), 
        additionalFilename = "../figures/cola_cluster1_FEAplot.pdf", 
        # height = 10,
        width = 11, 
        offsetValue = 5, angle = 90, 
        xleg = 0, yleg = 0)

plotFEA(dataFEA = cluster2_FEA[index2,],
        topBP = nrow(cluster2_FEA[index2,]), 
        additionalFilename = "../figures/cola_cluster2_FEAplot.pdf", 
        # height = 10,
        width = 11, 
        offsetValue = 5, angle = 90, 
        xleg = 0, yleg = 0)
