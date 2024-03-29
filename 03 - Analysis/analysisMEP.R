library(cglasso)
library(VennDiagram)

source("../01 - RCode/datajcggm.R")
source("../01 - RCode/jcglasso.R")
source("../01 - RCode/jcggm.R")
source("../01 - RCode/gof.R")
source("../01 - RCode/to_graph.R")

dati <- t(read.table("data/GSE79331_non-normalized.txt", header = TRUE, row.names = 1))
dim(dati)
head(dati)

apply(dati, 2, function(x) table(factor(1*(x == 40), levels = c(0, 1))))

# Non-informative assays with standard deviation of zero were excluded (n = 3).
idZeroVariance <- which(apply(dati, 2, sd) == 0)
dati2 <- dati[, -idZeroVariance]

# Cells were excluded if more than 70 assays did not result in amplification in that cell. 
id70 <- which(apply(dati2, 1, function(x) sum(x == 40)) >= 70)
dati3 <- dati2[-id70, ]

# Cells that displayed low levels of B2M or GAPDH were excluded, using cutoffs of 13 and 15 Ct cycles, respectively.
par(mfrow=c(1,2))
hist(dati2[, "B2M"], breaks = 100); abline(v = 13, lty = 3, col = 2)
hist(dati2[, "GAPDH"], breaks = 100); abline(v = 15, lty = 3, col = 2)
idB2MandGAPDH <- which(dati3[, "B2M"] >= 13 | dati3[, "GAPDH"] >= 15)
dati4 <- dati3[-idB2MandGAPDH, ]

# Finally, the mean Ct for each cell was calculated including Ct values for all assays that yielded detectable amplification, 
# and cells that displayed mean Ct value greater than 20 were removed after visual inspection of the data to exclude outliers.
par(mfrow=c(1,1))
hist(apply(dati4, 1, function(x) mean(x[x != 40])))
sum(apply(dati4, 1, function(x) mean(x[x != 40])) >= 20)

dim(dati4)
# [1] 681  87

# Normalization was performed to the mean of B2M and GAPDH (as per the single-cell analysis)
offset <- rowMeans(dati4[, c("B2M", "GAPDH")])
hist(offset)

# Ct values were normalized to the mean of B2M and GAPDH expression. 
# Normalized Ct values and raw data are listed in Additional file 2: Table S2 and S3, respectively.
dati5 <- t(sapply(1:nrow(dati4), function(i) {
  dati4[i, dati4[i, ] != 40] <- dati4[i, dati4[i, ] != 40] - offset[i]
  dati4[i, ]
}))
rownames(dati5) <- rownames(dati4)

set.seed(7)
PCAs2 <- prcomp(dati5, center = TRUE, scale. = TRUE)
objKmeans2 <- kmeans(PCAs2$x[,1:2], centers = 3, iter.max = 100, nstart = 50)
table(objKmeans2$cluster)
varExpl2 <- round(PCAs2$sdev^2 / sum(PCAs2$sdev^2) * 100, 2)
par(mfrow = c(1,2))
plot(PCAs2$x[,1:2], xlim = c(-6, 6), ylim = c(-10, 5), 
     xlab = "", ylab = "", col = objKmeans2$cluster, pch = 16)
mtext(text = paste0("Principal Component 1\n (", varExpl2[1], "% variance)"), side = 1, line = 3)
mtext(text = paste0("Principal Component 2\n (", varExpl2[2], "% variance)"), side = 2, line = 2)
plot(PCAs2$rotation[order(rowMeans(abs(PCAs2$rotation[,1:2])), decreasing = TRUE), 1:2][1:18,], cex = .1)
text(PCAs2$rotation[order(rowMeans(abs(PCAs2$rotation[,1:2])), decreasing = TRUE), 1:2][1:18, 1], 
     PCAs2$rotation[order(rowMeans(abs(PCAs2$rotation[,1:2])), decreasing = TRUE), 1:2][1:18, 2], 
     labels = rownames(PCAs2$rotation[order(rowMeans(abs(PCAs2$rotation[,1:2])), decreasing = TRUE), 1:2][1:18,]), 
     cex = 1)

heatmap(t(dati4[order(objKmeans2$cluster), 
                rownames(PCAs2$rotation[order(rowMeans(abs(PCAs2$rotation[,1:2])), decreasing = TRUE), 1:2][1:18,])]),
        scale = "none", Colv = NA, col = hcl.colors(25, "blues3", rev = FALSE))

dati6 <- t(sapply(1:nrow(dati5), function(i) {
  dati5[i, dati5[i, ] != 40] <- dati5[i, dati5[i, ] != 40] + offset[i]
  dati5[i, ]
}))
rownames(dati6) <- rownames(dati5)
# c("Pre-MEP", "E-MEP", "MK-MEP")

dati7 <- list("Pre-MEP" = dati6[objKmeans2$cluster == 1, ], 
              "E-MEP" = dati6[objKmeans2$cluster == 2, ], 
              "MK-MEP" = dati6[objKmeans2$cluster == 3, ])
offset2 <- split(offset, objKmeans2$cluster)
sapply(dati7, function(x) colMeans(x == 40))

# IDs of non-detected transcript in at least one of the MEP population + 2 housekeeping
id <- c(which(apply(sapply(dati7, 
                           function(y) apply(y, 2, 
                                             function(x) length(x) - sum(x == 40))), 1, 
                    function(z) any(z < 2))),
        match(c("B2M", "GAPDH"), colnames(dati7[[1]])))
# id <- which(apply(sapply(dati9, function(x) colMeans(is.na(x))) > .95, 1, any))
length(id)
nmsGenes2 <- colnames(dati7[[1]])[-id]
length(nmsGenes2)

# data description and checking censoring assumption
dati8 <- dati7
for(i in 1:length(dati7)) {
  dati8[[i]] <- dati8[[i]][, -id]
  colnames(dati8[[i]]) <- gsub("/", "", colnames(dati8[[i]]))
}

Z0 <- lapply(dati8, datacggm, up = 40)

up <- sapply(Z0, function(x) upper(x))
ymu <- sapply(Z0, function(x) ColMeans(x)$Y)
ysd <- sqrt(sapply(Z0, function(x) ColVars(x)$Y))
frc <- sapply(Z0, function(x) colMeans(event(x) == +1L))
prc <- pnorm(q = (up - ymu) / ysd, lower.tail = FALSE)

# pdf("ctvalues.pdf", width = .65*14, height = .65*6)
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
plot(ymu[, 1L], frc[, 1L], col = 1, cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex = 1.5,
     xlab = "Average cycle-threshold", ylab = "Proportion of missing data", main = "Pre-MEP")
points(ymu[, 1L], prc[, 1L], col = 2, pch = 4, cex = 1.5)
plot(ymu[, 2L], frc[, 2L], col = 1, cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex = 1.5,
     xlab = "Average cycle-threshold", ylab = "Proportion of missing data", main = "E-MEP")
points(ymu[, 2L], prc[, 2L], col = 2, pch = 4, cex = 1.5)
plot(ymu[, 3L], frc[, 3L], col = 1, cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex = 1.5,
     xlab = "Average cycle-threshold", ylab = "Proportion of missing data", main = "MK-MEP")
points(ymu[, 3L], prc[, 3L], col = 2, pch = 4, cex = 1.5)
# dev.off()

par(mfrow=c(1,1))
# pdf("histPRE.pdf", width = .8*14/3, height = .8*6)
hist(Z0[[1]], which = 55, max.hist = 1L, polygon.col = NA, legend = FALSE, xlab = "Ct values", 
     cex.axis = 1.6, cex.lab = 1.6, points.cex = 2, cex.main = 1.6, main = "Pre-MEP: MYB")
legend("topright", pch = c(4L, 1L), bty = "n", 
       pt.cex = 1.8, c("Proportion of\ncensored values\n", 
                       "Probability of a\ncensored value"))
# dev.off()
# pdf("histE.pdf", width = .8*14/3, height = .8*6)
hist(Z0[[2]], which = 55, max.hist = 1L, polygon.col = NA, legend = FALSE, xlab = "Ct values", 
     cex.axis = 1.6, cex.lab = 1.6, points.cex = 2, cex.main = 1.6, main = "E-MEP: MYB")
# dev.off()
# pdf("histMK.pdf", width = .8*14/3, height = .8*6)
hist(Z0[[3]], which = 55, max.hist = 1L, polygon.col = NA, legend = FALSE, xlab = "Ct values", 
     cex.axis = 1.6, cex.lab = 1.6, points.cex = 2, cex.main = 1.6, main = "MK-MEP: MYB")
# dev.off()

# normalizing data by using averages by tissues of the two housekeeping
dati9 <- dati8
for(i in 1:length(dati8)) {
  dati9[[i]][dati9[[i]] == 40] <- NA
  dati9[[i]] <- sweep(dati9[[i]], MARGIN = 1, STATS = offset2[[i]], FUN = "-")
}

loc <- read.csv2("../localizzazione.csv")
codifica <- loc$compartimento3
names(codifica) <- loc$X

# datajcggm data used for the selection of the optimal lambda and rho values
Z <- datajcggm(Y = lapply(dati9, \(x) x[, codifica == 1]), 
               X = lapply(dati9, \(x) x[, codifica == 0]))
summary(Z, which = 1L)
summary(Z, which = 2L)
summary(Z, which = 3L)

sapply(event2(Z), \(x) colMeans(x == 9))
sapply(ColMeans2(Z), unlist)

# datajcggm data used for the selection of the optimal nu value
ZX <- datajcggm(Y = lapply(dati9, \(x) x[, codifica == 0]))
summary(ZX, which = 1L)
summary(ZX, which = 2L)
summary(ZX, which = 3L)

names(Z) <- names(ZX) <- names(Z0)

# IPA and conversione gene names
conversion <- read.csv2("auxiliary/conversion_genes.csv", header = TRUE)
conversion$IDs2 <- gsub("/", "", conversion[, 1])
relationships <- read.csv2("auxiliary/Relationships.csv", header = TRUE)

data.frame.true <- data.frame(from = relationships$From.Molecule.s., to = relationships$To.Molecule.s.)
data.frame.true <- data.frame.true[data.frame.true[,1] != data.frame.true[,2], ]
for(i in seq_len(nrow(conversion))) {
  data.frame.true[, 1] <- gsub(conversion$IPA.names[i], conversion$IDs2[i], data.frame.true[, 1])
  data.frame.true[, 2] <- gsub(conversion$IPA.names[i], conversion$IDs2[i], data.frame.true[, 2])
}

# validated connections
g_true <- graph_from_data_frame(data.frame.true, directed = FALSE)
V(g_true)$color <- NA
V(g_true)$frame.color <- NA
V(g_true)$size <- 12
V(g_true)$label.cex <- 1
V(g_true)$label.font <- 2
V(g_true)$label.dist <- 0
E(g_true)$type <- "undirect"
E(g_true)$color <- "gray50"
E(g_true)$width <- 1
E(g_true)$arrow.mode <- 0

par(mfrow=c(1,1))
plot(g_true)

gtrue_link <- c(apply(data.frame.true, 1, function(x) paste0(x[1], ":", x[2])),
                apply(data.frame.true, 1, function(x) paste0(x[2], ":", x[1])))

sapply(c(.25, .5, .75), \(a) 
       jcglasso(data = ZX, nrho = 1L, alpha1 = a, trace = 1L, 
                thr.bcd = 1E-4, thr.em = 1)$rho)

nu.seq <- exp(seq(log(0.823198), log(0.823198 * 1E-6), l = 50L))

temp <- jcglasso(data = ZX, rho = nu.seq, alpha1 = .25, trace = 1L, 
                 thr.bcd = 1E-4, thr.em = 1E-3)
tempQfun <- QFun2(temp, mle = TRUE, trace = 1L)

par(mfrow=c(2,2))
plot(AIC(temp))
plot(AIC(temp, Qfun = tempQfun))

plot(BIC(temp))
plot(BIC(temp, Qfun = tempQfun))

temp.2 <- jcglasso(data = ZX, rho = nu.seq, alpha1 = .5, trace = 1L, 
                 thr.bcd = 1E-4, thr.em = 1E-3)
tempQfun.2 <- QFun2(temp.2, mle = TRUE, trace = 1L)

par(mfrow=c(2,2))
plot(AIC(temp.2))
plot(AIC(temp.2, Qfun = tempQfun.2))

plot(BIC(temp.2))
plot(BIC(temp.2, Qfun = tempQfun.2))

temp.3 <- jcglasso(data = ZX, rho = nu.seq, alpha1 = .75, trace = 1L, 
                   thr.bcd = 1E-4, thr.em = 1E-3)
tempQfun.3 <- QFun2(temp.3, mle = TRUE, trace = 1L)

par(mfrow=c(2,2))
plot(AIC(temp.3))
plot(AIC(temp.3, Qfun = tempQfun.3))

plot(BIC(temp.3))
plot(BIC(temp.3, Qfun = tempQfun.3))

summary(temp, GoF = BIC(temp, Qfun = tempQfun))
summary(temp.2, GoF = BIC(temp.2, Qfun = tempQfun.2))
summary(temp.3, GoF = BIC(temp.3, Qfun = tempQfun.3))

temp2 <- select.jcglasso(temp, BIC(temp, Qfun = tempQfun))
plot(temp2, type = "Gyy")
par(mfrow=c(1,3))
plot(to_graph2(temp2), type = "Gyy", highlight.connections = 3L)

temp2.2 <- select.jcglasso(temp.2, BIC(temp.2, Qfun = tempQfun.2))
plot(temp2.2, type = "Gyy")
par(mfrow=c(1,3))
plot(to_graph2(temp2.2), type = "Gyy", highlight.connections = 3L)

temp2.3 <- select.jcglasso(temp.3, BIC(temp.3, Qfun = tempQfun.3))
plot(temp2.3, type = "Gyy")
par(mfrow=c(1,3))
plot(to_graph2(temp2.3), type = "Gyy", highlight.connections = 3L)

# pdf(file = "bic1_200923.pdf", height = .7*10, width = .7*10)
par(mfrow=c(1,1))
plot(BIC(temp.3, Qfun = tempQfun.3), cex.axis = 1.2, cex.lab = 1.2,
     arg.points = list(cex = 1.5, pch = 17), main = "", xlab = expression(log(nu)))
# dev.off()

# optimal nu value = 0.06508239
nu.opt <- temp2.3$rho
nu.opt
alpha3.opt <- temp2.3$alpha1
alpha3.opt

rho.seq <- exp(seq(log(.25), log(.01), l = 10L))
# lambda.seq <- exp(seq(log(.25), log(.01), l = 10L))
lambda.seq <- exp(seq(log(.20), log(.01), l = 10L))
# lambda.seq <- exp(seq(log(.15), log(.01), l = 10L))

temp3 <- jcglasso(data = Z, 
                  rho = rho.seq, lambda = lambda.seq, nu = nu.opt, 
                  alpha1 = .25, alpha2 = .25, alpha3 = alpha3.opt,  
                  trace = 1L, thr.em = .01, thr.bcd = .001)
tempQfun3 <- QFun2(temp3, mle = TRUE, trace = 1L)

par(mfrow=c(2, 3))
plot(AIC(temp3))
plot(BIC(temp3))
plot(BIC(temp3, g = .5, type = "FD"))

plot(AIC(temp3, Qfun = tempQfun3))
plot(BIC(temp3, Qfun = tempQfun3))
plot(BIC(temp3, Qfun = tempQfun3, g = .5, type = "FD"))

temp3.2 <- jcglasso(data = Z, 
                    rho = rho.seq, lambda = lambda.seq, nu = nu.opt, 
                    alpha1 = .5, alpha2 = .5, alpha3 = alpha3.opt,  
                    trace = 1L, thr.em = .01, thr.bcd = .001)
tempQfun3.2 <- QFun2(temp3.2, mle = TRUE, trace = 1L)

par(mfrow=c(2, 3))
plot(AIC(temp3.2))
plot(BIC(temp3.2))
plot(BIC(temp3.2, g = .5, type = "FD"))

plot(AIC(temp3.2, Qfun = tempQfun3.2))
plot(BIC(temp3.2, Qfun = tempQfun3.2))
plot(BIC(temp3.2, Qfun = tempQfun3.2, g = .5, type = "FD"))


temp3.3 <- jcglasso(data = Z, 
                    rho = rho.seq, lambda = lambda.seq, nu = nu.opt, 
                    alpha1 = .75, alpha2 = .75, alpha3 = alpha3.opt,  
                    trace = 1L, thr.em = .01, thr.bcd = .001)
tempQfun3.3 <- QFun2(temp3.3, mle = TRUE, trace = 1L)

par(mfrow=c(2, 3))
plot(AIC(temp3.3))
plot(BIC(temp3.3))
plot(BIC(temp3.3, g = .5, type = "FD"))

plot(AIC(temp3.3, Qfun = tempQfun3.3))
plot(BIC(temp3.3, Qfun = tempQfun3.3))
plot(BIC(temp3.3, Qfun = tempQfun3.3, g = .5, type = "FD"))

summary(temp3, GoF = AIC(temp3, Qfun = tempQfun3), print.all = FALSE)
summary(temp3.2, GoF = AIC(temp3.2, Qfun = tempQfun3.2), print.all = FALSE)
summary(temp3.3, GoF = AIC(temp3.3, Qfun = tempQfun3.3), print.all = FALSE)

summary(temp3, GoF = BIC(temp3, Qfun = tempQfun3), print.all = FALSE)
summary(temp3.2, GoF = BIC(temp3.2, Qfun = tempQfun3.2), print.all = FALSE)
summary(temp3.3, GoF = BIC(temp3.3, Qfun = tempQfun3.3), print.all = FALSE)

# pdf(file = "bic2_200923.pdf", height = .7*10, width = .7*10)
par(mfrow=c(1,1))
plot(BIC(temp3.3, Qfun = tempQfun3.3), nlevels = 10, cex.axis = 1.2, cex.lab = 1.2,
     labcex = 1.1, arg.points = list(cex = 1.5, pch = 17), main = "", ylab = expression(log(rho)))
# dev.off()

temp4 <- select.jcglasso(temp3.3, GoF = BIC(temp3.3, Qfun = tempQfun3.3))
summary(temp4)

plot(temp4, type = "Gyy")
plot(temp4, type = "Gxy")
plot(temp4, type = "Gxx")

# MLE of the optimal model
temp_obj <- jcggm(temp4, trace = 1L)
summary(temp_obj)

par(mfrow=c(1,3))
plot(to_graph2(temp_obj), type = "bipartite", which = 1L)
plot(to_graph2(temp_obj), type = "bipartite", which = 2L)
plot(to_graph2(temp_obj), type = "bipartite", which = 3L)


# estimated connections
g_est <- getGraph2(to_graph2(temp_obj, weighted = TRUE), "bipartite")
names(g_est) <- c("Pre-MEP", "E-MEP", "MK-MEP")#names(Z)

gpre_link <- apply(as_data_frame(g_est$`Pre-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
ge_link <- apply(as_data_frame(g_est$`E-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
gmk_link <- apply(as_data_frame(g_est$`MK-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))

edge_venn <- list("Pre-MEP" = apply(as_data_frame(g_est$`Pre-MEP`)[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")), 
                  "E-MEP" = apply(as_data_frame(g_est$`E-MEP`)[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")),
                  "MK-MEP" = apply(as_data_frame(g_est$`MK-MEP`)[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")))
sapply(edge_venn, length)
round(sapply(edge_venn, length) / (nresp(Z[[1]]) * (nresp(Z[[1]])-1) / 2 + npred(Z[[1]]) * (npred(Z[[1]])-1) / 2 + npred(Z[[1]]) * nresp(Z[[1]])) * 100, 2L)

a <- venn.diagram(
  x = edge_venn,
  main = "All links",
  main.cex = 1.25,
  filename = NULL,
  output = TRUE ,
  imagetype = "png" ,
  height = .7*960 , 
  width = .7*960 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  fill = c("red", "blue", "green"),
  alpha = .35,
  cex = 1.25,
  fontfamily = "sans",
  cat.cex = .9,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

edge_venn2 <- list("Pre-MEP" = sort(gpre_link[gpre_link %in% gtrue_link]), 
                  "E-MEP" = sort(ge_link[ge_link %in% gtrue_link]),
                  "MK-MEP" = sort(gmk_link[gmk_link %in% gtrue_link]))
sapply(edge_venn2, length)
round(sapply(edge_venn2, length) / (nresp(Z[[1]]) * (nresp(Z[[1]])-1) / 2 + npred(Z[[1]]) * (npred(Z[[1]])-1) / 2 + npred(Z[[1]]) * nresp(Z[[1]])) * 100, 2L)
round(sapply(edge_venn2, length) / sapply(edge_venn, length) * 100, 2L)

b <- venn.diagram(
  x = edge_venn2,
  main = "Validated links",
  main.cex = 1.25,
  filename = NULL,
  output = TRUE ,
  imagetype = "png" ,
  height = .7*960 , 
  width = .7*960 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  fill = c("red", "blue", "green"),
  alpha = .35,
  cex = 1.25,
  fontfamily = "sans",
  cat.cex = .9,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

par(mfrow=c(1,1))
# pdf("venn_diagram_links_mep_200923.pdf", width = .7*12, height = .7*8)
plot.new()
pushViewport(plotViewport(layout = grid.layout(1, 2)))
pushViewport(plotViewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(a)
popViewport()
pushViewport(plotViewport(layout.pos.col = 2, layout.pos.row = 1))
grid.draw(b)
# dev.off()

gpre_link <- edge_venn2$`Pre-MEP`
ge_link <- edge_venn2$`E-MEP`
gmk_link <- edge_venn2$`MK-MEP`

gpre <- intersection(g_est$`Pre-MEP`, 
                     graph_from_data_frame(
                       as.data.frame(matrix(unlist(strsplit(gpre_link, ":")), 
                                            ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to"))))))
gpre <- delete.vertices(gpre, degree(gpre) == 0)
ge <- intersection(g_est$`E-MEP`, 
                   graph_from_data_frame(
                     as.data.frame(matrix(unlist(strsplit(ge_link, ":")), 
                                          ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to"))))))
ge <- delete.vertices(ge, degree(ge) == 0)
gmk <- intersection(g_est$`MK-MEP`, 
                    graph_from_data_frame(
                      as.data.frame(matrix(unlist(strsplit(gmk_link, ":")), 
                                           ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to"))))))
gmk <- delete.vertices(gmk, degree(gmk) == 0)

g_pre_e_mk <- union(gpre, ge, gmk)
g_pre_e_mk <- delete.vertices(g_pre_e_mk, degree(g_pre_e_mk) == 0)
g_pre_e_mk_int <- intersection(gpre, ge, gmk)
g_pre_e_mk_int <- delete.vertices(g_pre_e_mk_int, degree(g_pre_e_mk_int) == 0)
common_link <- apply(as_data_frame(g_pre_e_mk_int, what = "both")$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":"))

set.seed(5)
coords <- layout_(g_pre_e_mk, with_dh())
rownames(coords) <- V(g_pre_e_mk)$name
# plot(g_pre_e_mk, layout = coords)
plot(coords); text(coords[,1], coords[,2], labels = rownames(coords))

V(gmk)$size <- 20
V(ge)$size <- 20
V(gpre)$size <- 20
V(gmk)$label.cex <- 1.4
V(ge)$label.cex <- 1.4
V(gpre)$label.cex <- 1.4

gpre2 <- as_data_frame(gpre, what = "both")
gpre2$edges[!(apply(gpre2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- ifelse(gpre2$edges[!(apply(gpre2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] == "gray85", "green", "red")
gpre2$edges[(apply(gpre2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- "gray70"
gpre2$vertices <- rbind.data.frame(gpre2$vertices, 
                                   data.frame(type = "response", 
                                              color = NA, frame.color = NA, size = 20, 
                                              label.cex = 1.4, label.font = 2, label.dist = 0, 
                                              label.color = "gray65", 
                                              name = setdiff(V(g_pre_e_mk)$name, gpre2$vertices$name),
                                              row.names = setdiff(V(g_pre_e_mk)$name, gpre2$vertices$name)))

ge2 <- as_data_frame(ge, what = "both")
ge2$edges[!(apply(ge2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- ifelse(ge2$edges[!(apply(ge2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] == "gray85", "green", "red")
ge2$edges[(apply(ge2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- "gray70"
ge2$vertices <- rbind.data.frame(ge2$vertices, 
                                 data.frame(type = "response", 
                                            color = NA, frame.color = NA, size = 20, 
                                            label.cex = 1.4, label.font = 2, label.dist = 0, 
                                            label.color = "gray65", 
                                            name = setdiff(V(g_pre_e_mk)$name, ge2$vertices$name),
                                            row.names = setdiff(V(g_pre_e_mk)$name, ge2$vertices$name)))

gmk2 <- as_data_frame(gmk, what = "both")
gmk2$edges[!(apply(gmk2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- ifelse(gmk2$edges[!(apply(gmk2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] == "gray85", "green", "red")
gmk2$edges[(apply(gmk2$edges[,1:2], 1, \(x) paste(x[1], x[2], sep = ":")) %in% common_link), "color"] <- "gray70"
gmk2$vertices <- rbind.data.frame(gmk2$vertices, 
                                  data.frame(type = "response", 
                                             color = NA, frame.color = NA, size = 20, 
                                             label.cex = 1.4, label.font = 2, label.dist = 0, 
                                             label.color = "gray65", 
                                             name = setdiff(V(g_pre_e_mk)$name, gmk2$vertices$name),
                                             row.names = setdiff(V(g_pre_e_mk)$name, gmk2$vertices$name)))

for(i in seq_len(nrow(conversion))){
  V(g_pre_e_mk)$name <- gsub(conversion$IDs2[i], conversion$IDs[i], V(g_pre_e_mk)$name)
  
  gpre2$vertices$name <- gsub(conversion$IDs2[i], conversion$IDs[i], gpre2$vertices$name)
  rownames(gpre2$vertices) <- gsub(conversion$IDs2[i], conversion$IDs[i], rownames(gpre2$vertices))
  gpre2$edges$from <- gsub(conversion$IDs2[i], conversion$IDs[i], gpre2$edges$from)
  gpre2$edges$to <- gsub(conversion$IDs2[i], conversion$IDs[i], gpre2$edges$to)
  
  ge2$vertices$name <- gsub(conversion$IDs2[i], conversion$IDs[i], ge2$vertices$name)
  rownames(ge2$vertices) <- gsub(conversion$IDs2[i], conversion$IDs[i], rownames(ge2$vertices))
  ge2$edges$from <- gsub(conversion$IDs2[i], conversion$IDs[i], ge2$edges$from)
  ge2$edges$to <- gsub(conversion$IDs2[i], conversion$IDs[i], ge2$edges$to)
  
  gmk2$vertices$name <- gsub(conversion$IDs2[i], conversion$IDs[i], gmk2$vertices$name)
  rownames(gmk2$vertices) <- gsub(conversion$IDs2[i], conversion$IDs[i], rownames(gmk2$vertices))
  gmk2$edges$from <- gsub(conversion$IDs2[i], conversion$IDs[i], gmk2$edges$from)
  gmk2$edges$to <- gsub(conversion$IDs2[i], conversion$IDs[i], gmk2$edges$to)
}

gpre2$vertices[gpre2$vertices$name %in% names(codifica[codifica == 1]), "label.color"] <- "blue"
gpre2$vertices[gpre2$vertices$name %in% names(codifica[codifica == 0]), "label.color"] <- "black"
gpre2$vertices[names(which(degree(graph_from_data_frame(gpre2$edges, vertices = gpre2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"] <- adjustcolor(gpre2$vertices[names(which(degree(graph_from_data_frame(gpre2$edges, vertices = gpre2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"], .5)

ge2$vertices[ge2$vertices$name %in% names(codifica[codifica == 1]), "label.color"] <- "blue"
ge2$vertices[ge2$vertices$name %in% names(codifica[codifica == 0]), "label.color"] <- "black"
ge2$vertices[names(which(degree(graph_from_data_frame(ge2$edges, vertices = ge2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"] <- adjustcolor(ge2$vertices[names(which(degree(graph_from_data_frame(ge2$edges, vertices = ge2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"], .5)

gmk2$vertices[gmk2$vertices$name %in% names(codifica[codifica == 1]), "label.color"] <- "blue"
gmk2$vertices[gmk2$vertices$name %in% names(codifica[codifica == 0]), "label.color"] <- "black"
gmk2$vertices[names(which(degree(graph_from_data_frame(gmk2$edges, vertices = gmk2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"] <- adjustcolor(gmk2$vertices[names(which(degree(graph_from_data_frame(gmk2$edges, vertices = gmk2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])) == 0)), "label.color"], .5)


gpre2 <- graph_from_data_frame(gpre2$edges, vertices = gpre2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])
ge2 <- graph_from_data_frame(ge2$edges, vertices = ge2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])
gmk2 <- graph_from_data_frame(gmk2$edges, vertices = gmk2$vertices[V(g_pre_e_mk)$name, c(9, 1:8)])

m <- matrix(0, 4, 4); m[1:2, 1:2] <- 1; m[1:2, 3:4] <- 2; m[3:4, 2:3] <- 3
# pdf("links_mep_only_200923.pdf", width = .8*12, height = .8*12)
layout(m)
par(mar = c(2.1, 2.7, 2.1, 2.7))
plot(gpre2, layout = coords, xlim = c(-.85,.85), ylim = c(-.85,.85))
title("Pre-MEP", cex.main = 1.5)
plot(ge2, layout = coords, xlim = c(-.85,.85), ylim = c(-.85,.85))
title("E-MEP", cex.main = 1.5)
plot(gmk2, layout = coords, xlim = c(-.85,.85), ylim = c(-.85,.85))
title("MK-MEP", cex.main = 1.5)
# dev.off()


### Table 2
sapply(c("Theta", "B", "Omega"), \(type) {
  temp <- sapply(1:length(Z), \(i) {
    x <- coef(temp_obj, type = type, rho.id = 1L, lambda.id = 1L, class.id = i)
    if(type == "B") sum(x[-1,] != 0) else sum(x[upper.tri(x, diag = FALSE)] != 0)
  })
  names(temp) <- names(Z)
  rbind(temp, round(temp / switch(type, 
                                  "Theta" = nresp(Z) * (nresp(Z) - 1) / 2,
                                  "Omega" = npred(Z) * (npred(Z) - 1) / 2,
                                  "B" = nresp(Z) * npred(Z)) * 100, 2L), deparse.level = 0)
}, simplify = FALSE)

sapply(c("Gyy", "Gxy", "Gxx"), \(type) {
  g_est2 <- getGraph2(to_graph2(temp_obj, weighted = TRUE), type = type)
  names(g_est2) <- names(Z)
  
  gpre_link2 <- apply(as_data_frame(g_est2$`Pre-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  ge_link2 <- apply(as_data_frame(g_est2$`E-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  gmk_link2 <- apply(as_data_frame(g_est2$`MK-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  
  edge_venn3 <- list("Pre-MEP" = sort(gpre_link2[gpre_link2 %in% gtrue_link]), 
                     "E-MEP" = sort(ge_link2[ge_link2 %in% gtrue_link]),
                     "MK-MEP" = sort(gmk_link2[gmk_link2 %in% gtrue_link]))
  rbind(sapply(edge_venn3, length),
        round(sapply(edge_venn3, length) / switch(type, 
                                                  "Gyy" = nresp(Z) * (nresp(Z) - 1) / 2,
                                                  "Gxx" = npred(Z) * (npred(Z) - 1) / 2,
                                                  "Gxy" = nresp(Z) * npred(Z)) * 100, 2L))
}, simplify = FALSE)

sapply(edge_venn, length)
round(sapply(edge_venn, length) / (nresp(Z) * (nresp(Z)-1) / 2 + npred(Z) * (npred(Z)-1) / 2 + npred(Z) * nresp(Z)) * 100, 2)

sapply(edge_venn2, length)
round(sapply(edge_venn2, length) / (nresp(Z) * (nresp(Z)-1) / 2 + npred(Z) * (npred(Z)-1) / 2 + npred(Z) * nresp(Z)) * 100, 2)

df_pre_MEP <- sapply(c("Gyy", "Gxy", "Gxx"), \(type) {
  g_est2 <- getGraph2(to_graph2(temp_obj, weighted = TRUE), type = type)
  names(g_est2) <- names(Z)
  
  gpre_link2 <- apply(as_data_frame(g_est2$`Pre-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  
  as_data_frame(g_est2$`Pre-MEP`)[gpre_link2 %in% gtrue_link, c(1,2,8)]
}, simplify = FALSE)

df_e_MEP <- sapply(c("Gyy", "Gxy", "Gxx"), \(type) {
  g_est2 <- getGraph2(to_graph2(temp_obj, weighted = TRUE), type = type)
  names(g_est2) <- names(Z)
  
  ge_link2 <- apply(as_data_frame(g_est2$`E-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  
  as_data_frame(g_est2$`E-MEP`)[ge_link2 %in% gtrue_link, c(1,2,8)]
}, simplify = FALSE)

df_mk_MEP <- sapply(c("Gyy", "Gxy", "Gxx"), \(type) {
  g_est2 <- getGraph2(to_graph2(temp_obj, weighted = TRUE), type = type)
  names(g_est2) <- names(Z)
  
  gmk_link2 <- apply(as_data_frame(g_est2$`MK-MEP`)[, 1:2], 1, \(x) paste0(x[1], ":", x[2]))
  
  as_data_frame(g_est2$`MK-MEP`)[gmk_link2 %in% gtrue_link, c(1,2,8)]
}, simplify = FALSE)

do.call(rbind, df_pre_MEP)

# clip <- pipe("pbcopy", "w")
# write.table(do.call(rbind, lapply(df_pre_MEP, \(x){
#   x <- x[order(abs(x$weight), decreasing = TRUE), ]
#   x$weight <- round(x$weight, digits = 3L)
#   x
# })), file = clip, sep="\t", row.names = TRUE)
# close(clip)

do.call(rbind, df_e_MEP)

# clip <- pipe("pbcopy", "w")
# write.table(do.call(rbind, lapply(df_e_MEP, \(x){
#   x <- x[order(abs(x$weight), decreasing = TRUE), ]
#   x$weight <- round(x$weight, digits = 3L)
#   x
# })), file = clip, sep="\t", row.names = TRUE)
# close(clip)

do.call(rbind, df_mk_MEP)

# clip <- pipe("pbcopy", "w")
# write.table(do.call(rbind, lapply(df_mk_MEP, \(x){
#   x <- x[order(abs(x$weight), decreasing = TRUE), ]
#   x$weight <- round(x$weight, digits = 3L)
#   x
# })), file = clip, sep="\t", row.names = TRUE)
# close(clip)
