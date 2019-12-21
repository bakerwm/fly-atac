library(mclust)
args = commandArgs(TRUE)

#Prefix of current experiment
currexpt = args[1]
#cellfloor = No. reads to be considered a cell
#Either specify a number or "mclust" to calculate automatically
cellfloor = args[2]

# report2 = read.table(paste0(currexpt,".report.txt"), header=T)
report2 = read.table(currexpt, header = TRUE)
print(currexpt)
report2$Tag = report2$Experiment
report2$Total = report2$ReadCount
bkgd.ind = grep("bkgd", report2$Tag)
if(length(bkgd.ind) > 0){
  nobkgdmat = report2[-bkgd.ind,]
} else {
  nobkgdmat = report2
}
cellcall = Mclust(data.frame(log10(nobkgdmat$Total)), G=2)
if(cellfloor == "auto"){
  cellfloor = min(nobkgdmat[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05), "Total"])
}else{
  cellfloor = as.numeric(cellfloor)
}
#cellfloor = 1000

subsetmat = nobkgdmat[which(nobkgdmat$Total >= cellfloor),]
write.table(cbind(as.character(rownames(subsetmat)),as.character(subsetmat[,1])),paste0(currexpt,".readdepth.cells.indextable.txt"),row.names=F,col.names=F,sep="\t",quote=F)
subsamples = levels(nobkgdmat$Tag)

# pdf(paste0(currexpt,".results.hists.pdf"),height=12,width=12)
pdfFile <- gsub(".report.txt", ".results.hists.pdf", currexpt)
pdf(pdfFile, width = 12, height = 12)
par(mfrow=c(2,2))
for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  currind = grep(subsamples[i], nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  hist(log10(currsub$Total),breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2)
  legend("topright",c(paste0("Total Reads: ", sum(currsub$Total)),
                      paste0("\n Total Reads (cells only): ", sum(currsub[currsubcells, "Total"])),
                      paste0("\n Total Barcodes: ",length(currsub$Total)),
                      paste0("\n Number of Cells: ",length(subsetmat$Total[currsubcells])),
                      paste0("\n Median Reads/Cell: ",median(currsub[currsubcells, "Total"])),
                      paste0("\n Range of Reads/Cell: ",min(currsub[currsubcells, "Total"])," - ",max(currsub[currsubcells, "Total"]))),bty="n")
}

hist(log10(nobkgdmat$Total),breaks=60,col="mediumseagreen",main="Overall",
     xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
abline(v=log10(cellfloor),lty="dashed",  lwd=2)
legend("topright",c(paste0("Total Reads: ",sum(nobkgdmat$Total)),
                    paste0("Total Reads (cells only): ",sum(subsetmat$Total)),
                    paste0("\n Total Barcodes: ",length(nobkgdmat$Total)),
                    paste0("\n Number of Cells: ",length(subsetmat$Total)),
                    paste0("\n Median Reads/Cell: ",median(subsetmat[, "Total"])),
                    paste0("\n Range of Reads/Cell: ",min(subsetmat[, "Total"])," - ",max(subsetmat[, "Total"]))),bty="n")
dev.off()
