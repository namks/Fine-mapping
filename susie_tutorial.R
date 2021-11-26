###############################
## Fine-mapping using susieR ##
## FOR ONLY ONE PHENO/REGION ##
###############################

# Import packages
library(data.table)
library(seqminer)
library(susieR)

# Define the window size 
ws <- 100000    # (Total window size = ws * 2)
phenocode <- "1"    # 1 is "systolic blood pressure"

# Read phenotype and GWAS summary statistics
fname <- paste0("/media/leelabsg-storage0/KBN_WORK/plinkfile/pheno/KBN_SAIGE_", phenocode, "_CA2.txt")
f <- fread(fname)
ssname <- paste0("/media/leelabsg-storage0/pheweb/data/SAIGE.pheno", phenocode, ".ca2.header.txt.gz")
ss <- fread(ssname)

# Calculate z-scores
z_scores <- ss$beta / ss$SE

# Pick the top SNP and region
top_idx <- which(abs(ss$Tstat) == max(abs(ss$Tstat)))
chr <- ss[top_idx,]$CHROM
pos <- ss[top_idx,]$POS
top_reg_idx <- which((ss$CHROM == chr) & (abs(ss$POS - pos) < ws))
top_reg <- ss[top_reg_idx,]
z <- z_scores[top_reg_idx]

# Read VCF files using seqminer
start_pos <- min(top_reg$POS)
end_pos <- max(top_reg$POS)
vcfname <- paste0("/media/leelabsg-storage0/KBN/유전정보/KCHIP_72298/CHR", chr, "_annoINFO_filINFO0.8_72K.vcf.gz")
pos_string <- paste0(chr, ":", start_pos, "-", end_pos)
vcf <- readVCFToMatrixByRange(vcfname, pos_string)
geno <- t(vcf[[1]])
geno_sort <- geno[sort(rownames(geno)),]

# Remove NA
no_na_idx <- which(is.na(f$y_quant) == F)
Y <- f[no_na_idx,]$y_quant
X <- geno_sort[no_na_idx,]

# Calculate pairwise correlation matrix
R <- cor(X)
fitted_rss <- susie_rss(z, R, L = 10)

# Save PIP plot as a png file
png(filename="susie_plot.png", width = 1000, height = 500)
susie_plot(fitted_rss, y="PIP")
dev.off()
