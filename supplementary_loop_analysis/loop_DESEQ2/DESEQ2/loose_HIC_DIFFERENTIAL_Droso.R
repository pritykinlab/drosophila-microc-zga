#### FOR DROSOPHILA

library(dplyr)
library(glue)
library(DESeq2)


col.names=c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'fc')

final_table <- read.csv('LOOSE_hic_values_close.csv')
print("USING HIC VALUES – WE CAN ALSO USE HIC VALUES CLOSE!!")

conditions <- colnames(final_table)[1:(length(colnames(final_table))-2)]
countData <- final_table[1:length(conditions)]


shifted <- final_table$shifted
places <- final_table$places

unparsed_conditions <- strsplit(conditions, '_')
parsed_conditions <- c()
for(i in unparsed_conditions){
	v <- i[1]
	if(v == 'WT' | v =='yw'){
		v <- 'yw'
	}
	parsed_conditions <- c(parsed_conditions, v)
}

colData <- cbind(colnames(countData), parsed_conditions)
colnames(colData) <- c('name', 'condition')
rownames(countData) <- paste(shifted, places, sep="_")

indices <- rownames(countData)
split <- strsplit(indices, '_')

starts <- as.numeric(matrix(unlist(split), ncol = length(split), byrow=FALSE)[3, ])
ends <- as.numeric(matrix(unlist(split), ncol = length(split), byrow=FALSE)[6, ])
dists <- (ends-starts)/1000
quartiles <- ntile(dists, 4)

# normalizationFactors(dds) <- normFactors
# need to get size factors too!

process_conditions <- function(c1, c2, dds, prefix, output='output'){
	res = results(dds, contrast=c("condition",  c1, c2))
	sum(res[!is.na(res$padj), ]$padj < .05)
	path <- glue('plots/{prefix}_dispersion_estimation_{c1}_vs_{c2}.png')
	png(path)
	plotDispEsts(dds)
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_unshifted.png')
	png(path)
	plotMA(res[shifted==0,], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_shifted.png')
	png(path)
	plotMA(res[shifted==1,], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_shifted_quart1.png')
	png(path)
	plotMA(res[(shifted==1)&(quartiles==1),], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_shifted_quart2.png')
	png(path)
	plotMA(res[(shifted==1)&(quartiles==2),], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_shifted_quart3.png')
	png(path)
	plotMA(res[(shifted==1)&(quartiles==3),], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}_shifted_quart4.png')
	png(path)
	plotMA(res[(shifted==1)&(quartiles==4),], ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()

	path <- glue('plots/{prefix}_ma_plot_{c1}_vs_{c2}.png')
	png(path)
	plotMA(res, ylim=c(-2,2), xlim=c(.1, 40000))
	dev.off()
	path <- glue("{output}/{prefix}_{c1}_vs_{c2}.csv")
	write.table(res, path)
	print(glue("writing to {path}"))
	return(res)
}

run_deseq <- function(normFactors, countData, prefix, output='output'){
	n <- length(colnames(normFactors))
	combined_scaling <- 1/normFactors 

	not_na <- (rowSums(is.na(combined_scaling))==0)
	# combined_scaling <- combined_scaling[not_na, ] 
	# sub_combined_scaling <- combined_scaling
	# sub_combined_scaling['nc14_1'] <- sub_combined_scaling['nc14_1_1'] + sub_combined_scaling['nc14_2_1'] 
	# sub_combined_scaling['nc14_2'] <- sub_combined_scaling['nc14_1_2'] + sub_combined_scaling['nc14_2_2'] 
	# sub_combined_scaling['nc1.8_1'] <- sub_combined_scaling['nc1.8_1_1'] + sub_combined_scaling['nc1.8_2_1'] 
	# sub_combined_scaling['nc1.8_2'] <- sub_combined_scaling['nc1.8_1_2'] + sub_combined_scaling['nc1.8_2_2']
	# sub_combined_scaling['s10.12_1'] <- sub_combined_scaling['s10.12_1_1'] + sub_combined_scaling['s10.12_2_1'] 
	# sub_combined_scaling['s10.12_2'] <- sub_combined_scaling['s10.12_1_2'] + sub_combined_scaling['s10.12_2_2']
	# sub_combined_scaling <- subset(sub_combined_scaling, select = c(nc14_1, nc14_2, nc1.8_1, nc1.8_2, s10.12_1, s10.12_2))

	sub_countData <- countData[not_na, ]
	# sub_countData['nc14_1'] <- sub_countData['nc14_1_1'] + sub_countData['nc14_2_1'] 
	# sub_countData['nc14_2'] <- sub_countData['nc14_1_2'] + sub_countData['nc14_2_2'] 
	# sub_countData['nc1.8_1'] <- sub_countData['nc1.8_1_1'] + sub_countData['nc1.8_2_1'] 
	# sub_countData['nc1.8_2'] <- sub_countData['nc1.8_1_2'] + sub_countData['nc1.8_2_2']
	# sub_countData['s10.12_1'] <- sub_countData['s10.12_1_1'] + sub_countData['s10.12_2_1'] 
	# sub_countData['s10.12_2'] <- sub_countData['s10.12_1_2'] + sub_countData['s10.12_2_2']
	# sub_countData <- subset(sub_countData, select = c(nc14_1, nc14_2, nc1.8_1, 
	# 												nc1.8_2, s10.12_1, s10.12_2
													# nc14.FED.rerun_1_1, nc14.FED.rerun_1_2, 
													# nc14.FED_1_1, nc14.FED_1_2,
													# nc12.mitotic.FED_1_1, nc12.mitotic.FED_1_2 
													# nc12.mitotic.FED.rerun_1_1, nc12.mitotic.FED.rerun_1_2 
													# ))

	shifted <- shifted[not_na]
	quartiles <- quartiles[not_na]
	nonmerge_sub_countData <- countData[not_na, ]

	sub_unparsed_conditions <- c()
	sub_parsed_conditions <- c()
	for(i in colnames(sub_countData)){
		v <- strsplit(i, '_')[[1]]
		n = length(v)
		v <- paste(v[1:n-1], collapse='_')

		sub_parsed_conditions <- c(sub_parsed_conditions, v)
		sub_unparsed_conditions <- c(sub_unparsed_conditions, i)
	}
	condition <- sub_parsed_conditions
	sub_colData <- cbind(colnames(sub_countData), condition)
	dds = DESeqDataSetFromMatrix(countData = nonmerge_sub_countData, colData = colData, design = ~condition)
	normalizationFactors(dds) <- as.matrix(combined_scaling)
	dds <- DESeq(dds)
	prefix <- 'loose'
	resultsNames(dds)
	all_conditions <- unique(colData[, 'condition'])

	 for(condition1 in all_conditions){
	 	for(condition2 in all_conditions){
	 		if(condition1 != condition2){
	 			res <- process_conditions(condition1, condition2, dds, prefix)
	 		}
	 	}
	 }
	return(dds)
}



# normFactors <- read.csv('balanced_normalization_df.csv')
# normFactors <- normFactors[c(-13, -14, -17, -18)]
normFactors <- read.csv('LOOSE_raw_normalization_df.csv')
dds <- run_deseq(normFactors, countData, 'LOOSE_normalization_factors', output='replicate_results')
# sizefacs <- sizeFactors(dds)
# write.csv(sizefacs, './sizefacs_test.csv')

save.image(file='rdata_for_evan.RData') 

# sizefacs <- sizeFactors(dds)
# write.csv(sizefacs, './sizefacs_test.csv')

