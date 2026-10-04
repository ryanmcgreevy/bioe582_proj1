args <- commandArgs(trailingOnly = TRUE)

library(affy)
library(affydata)
library(hgu133a.db)
input_txt = args[1] 
output_path = args[2]
dir_name = dirname(input_txt)


samples = read.table(input_txt,header=T)

sample_files = paste(samples$sampleId,".CEL",sep="")
rawdata = ReadAffy(filenames=sample_files,celfile.path=dir_name,phenoData=samples,sampleNames=samples$sampleId)

eset = rma(rawdata,normalize=T,background=T)

pdf(paste(output_path,'/boxplot_raw_vs_normalized.pdf',sep=''))
par(mfrow=c(1,2))
boxplot(rawdata, col=c("red","red","blue","blue"), names=sampleNames(rawdata), main="Pre-Normalization Boxplot")
boxplot(exprs(eset), col=c("red","red","blue","blue"), names=sampleNames(eset), main="Post-Normalization Boxplot")
dev.off()

mas5 = mas5calls(rawdata,tau=0.015,alpha1=0.04,alpha2=0.06)
pma = exprs(mas5)

all_names = row.names(pma)
#get probe sets that are P in at least 1 sample of each phenotype
mas5filtered = pma[(pma[,1] == "P" | pma[,2]=="P") & (pma[,3] == "P" | pma[,4]=="P"),]
#names of probes present
present_names = row.names(mas5filtered)
#names of probes which were filtered out
filtered_out = pma[!(all_names %in% present_names),]
out_file = file(paste(output_path,"/mas5_filtered_out.txt",sep=""),"w")
#for each probe set filtered out, sum the numbers of "absent" for each phenotype
for(probe in row.names(filtered_out))
	write(sprintf("%s\t%i\t%i",probe,sum(filtered_out[probe,1:2] == "A"),sum(filtered_out[probe,3:4]=="A")),file=out_file)
close(out_file)

#find control probe indices
matches = grep("AFFY*",present_names)
#filter out all control probes
no_control_and_mas5 = present_names[!(1:length(present_names) %in% matches)]

#get entrez and gene symbol
entrez = as.list(hgu133aENTREZID[no_control_and_mas5])
symbol = as.list(hgu133aSYMBOL[no_control_and_mas5])

#calculate the average expression for each probe set
es = exprs(eset)
avg1 = apply(es[no_control_and_mas5,1:2],1,mean)
avg0 = apply(es[no_control_and_mas5,3:4],1,mean)

#make dataframe for easier sorting
avg_and_names = data.frame(as.numeric(entrez),as.character(symbol),as.character(no_control_and_mas5),as.numeric(avg1),as.numeric(avg0))
#remove NA's (from EntrezID)
complete = avg_and_names[complete.cases(avg_and_names),]
#sort in ascending based on entrezid,gene symbol and probe id
sorted = complete[with(complete,order(complete[,1],as.character(complete[,2]),as.character(complete[,3]))),]
names(sorted)= c("entrez","symbol","probe_id","avg1","avg0")

#write out average expression
out_file2 = file(paste(output_path,"/annotated_probesets_avg_expr.txt",sep=""),"w")
for(i in 1:dim(sorted)[1]){
	#check if there is an "absent" sample, if so replace avg with '--'
	if(pma[as.character(sorted[i,3]),1] == "A" | pma[as.character(sorted[i,3]),2] == "A") {
		avg1i = "--"
	}else{avg1i = sorted[i,4]}
	if(pma[as.character(sorted[i,3]),3] == "A" | pma[as.character(sorted[i,3]),4] == "A") {
		avg0i = "--"
	}else{avg0i = sorted[i,5]}
	write(sprintf("%i\t%s\t%s\t%s\t%s",sorted[i,1],as.character(sorted[i,2]),as.character(sorted[i,3]),as.character(avg1i),as.character(avg0i)),file=out_file2)
}

close(out_file2)
