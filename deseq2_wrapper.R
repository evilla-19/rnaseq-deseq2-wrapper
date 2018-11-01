
#setwd(#Set your wd here)

#############################################
#####	#	Count file 				#########
#############################################

txt_files = list.files(pattern = 'Sample*.*geneCounts*')

data_list = lapply(txt_files, read.delim, sep = "\t", header = FALSE)

data_frame = data_list[[1]]

for (i in 2:length(data_list))
{
	data_frame = cbind(data_frame,data_list[[i]][,2] )
	}
	

colnames(data_frame)[2:length(data_frame)] = txt_files

data_frame = data_frame[grep('ENS*', data_frame[,1]),]

colnames(data_frame)[1] = ''


# file = 'WT_jq1_geneCounts.txt'
write.table(data_frame, file  = file, sep = '\t', row.names = FALSE, quote = FALSE)


#Libraries

library(DESeq)
library(DESeq2)
library(annotate)
library(AnnotationFuncs)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)


#############################################
#####	#	Wrapper DESeq 			#########
#############################################


# genelist_txt = 'WT_jq1_geneCounts.txt'
# control = 'CA1wtVEH'
# treated = 'CA1jq1'
# colnames(countData)


col = colorRampPalette(brewer.pal(10, 'RdBu'))(99)


wrapper_deseq = function(genelist_txt, control, treated, type = '', outlier = 'none')
{
	
#####################
#### Analysis #######
#####################

countData = read.table(genelist_txt, row.names = 1, header = TRUE)
countData = countData[,c(grep(as.character(control), colnames(countData)), grep(as.character(treated), colnames(countData)))]


#get rid of outliers if they are given
if (as.character(outlier) == 'none')
{
countData = countData
} else {
countData = countData[,-c(grep(outlier[1], colnames(countData)), grep(outlier[2], colnames(countData)))]
	}

colData = data.frame(condition = c(rep(as.character(control), length(grep(as.character(control), colnames(countData)))),rep(as.character(treated), length(grep(as.character(treated), colnames(countData))))), row.names = colnames(countData))
dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
colData(dds)$condition = factor(colData(dds)$condition,levels=c(as.character(control), as.character(treated)))

dds = DESeq(dds)
res = results(dds, cooksCutoff = FALSE)
res = as.data.frame(res)
genesymbols = unlist(AnnotationFuncs::translate(rownames(res), from = org.Mm.egENSEMBL2EG, to = org.Mm.egSYMBOL, reduce = 'first'))
genesymbols = data.frame(id = names(genesymbols), genesymbol = genesymbols)
res = data.frame(res, id = rownames(res))

res_annot = merge(res, genesymbols, by = 'id', all.x = TRUE)
res_annot = res_annot[order(res_annot$padj),]


write.table(file = paste(paste(paste(as.character(control), as.character(treated), sep = 'vs'),as.character(type), sep = '_'),'_deseq2.xls', sep=''),res_annot, sep = '\t',row.names = FALSE)

#####################
#### Plotting #######
#####################


rld = rlog(dds)
vsd = varianceStabilizingTransformation(dds)

distsRL = dist(t(assay(rld)))

mat = as.matrix(distsRL)

class(mat)

rownames(mat) = colnames(mat) 


pdf(file = paste(paste(paste(as.character(control), as.character(treated), sep = 'vs'), as.character(type), sep = '_'),'distace_plot.pdf', sep = '_'))
heatmap.2(mat, trace="none", margin=c(13, 13), col = col)
dev.off()


pdf(file = paste(paste(paste(as.character(control), as.character(treated), sep = 'vs'), as.character(type), sep = '_'),'pca.pdf', sep = '_'))
print(plotPCA(rld, intgroup=c("condition")))
dev.off()


objects = list(dds, res_annot)
save(dds, res_annot, file = paste(paste(paste(as.character(control), as.character(treated), sep = 'vs'), as.character(type), sep = '_'),'objects_deseq2', sep=''))


}


