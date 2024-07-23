## scRNA-seq GEX save
write.csv(as.character(pbmc@active.ident),"./filtered_gene_bc_matrices/meta.csv")
write.csv(as.matrix(pbmc@assays$RNA$scale.data),"./filtered_gene_bc_matrices/exp.csv")

## BUILD deconvolution dataset

cell <- read.csv("./filtered_gene_bc_matrices/meta.csv")
exp <- read.csv("./filtered_gene_bc_matrices/exp.csv")
rownames(exp) <- exp$X
exp <- exp[,-1]
ucell <- unique(cell$x)
ref <- {}
for(n in ucell){
  l <- which(cell$x == n)
  out <- apply(exp[,l], 1, median, na.rm=TRUE)
  ref <- cbind(ref, out)
}
colnames(ref) <- ucell
write.csv(ref,"./filtered_gene_bc_matrices/ref.csv")

## Deconvolution method
ref <- read.csv("./filtered_gene_bc_matrices/ref.csv")
rownames(ref) <- ref$X
ref <- ref[,-1]
pre <- read.csv("./filtered_gene_bc_matrices/tumor.csv")
pre <- pre[!duplicated(pre$X), ]
gene <- pre$X
pre <- as.matrix(pre[,-1])
rownames(pre) <- gene
ana <- merge(ref,pre, by = 'row.names')
gene <- ana[,1]
ana <- ana[,-1]
relation <- lm(V1 ~ .,ana)
process_proportion <- relation$coefficients[-1]
process_proportion[process_proportion<0] <- 0
final_proportion <- process_proportion/sum(process_proportion)
final_proportion
