inputfilename = "GSE25066"
clinicaldatafilename = "GSE25066_Clinical_Data"

lrmc_score = function(x){
	sc = 0
	mx = max(x,na.rm=TRUE)
	pmx = which(x==mx)
	pmx = pmx[length(pmx)]
	
	if(pmx>1 && pmx<length(x)){
		for(i in 2:pmx){
			if(x[i]>=x[i-1]){
				sc = sc + 1
			}else{
				sc = sc - 1
			}
		}
		for(i in (pmx+1):length(x)){
			if(x[i]<=x[i-1]){
				sc = sc + 1
			}else{
				sc = sc - 1
			}
		}
	}else if(pmx==1){
		for(i in 2:length(x)){
			if(x[i]<=x[i-1]){
				sc = sc + 1
			}else{
				sc = sc - 1
			}
		}		
	}else{
		for(i in 2:length(x)){
			if(x[i]>=x[i-1]){
				sc = sc + 1
			}else{
				sc = sc - 1
			}
		}
	}
	
	return(sc/(length(x)-1))
}


lrmc = function(ex1,gene,osd,evos){
	lrps = c()
	thrs = c()
	cpcs = c()
	
	sortedex1 = unique(sort(ex1))
	
	for(j in 3:(length(ex1)-3)){
		ex6 = ex1
		
		ex6[which(ex1<=sortedex1[j])] = 0
		ex6[which(ex1>sortedex1[j])] = 1
		
		if(length(unique(ex6))>1){
			lr3 = survdiff(Surv(osd,evos)~ex6)
			lrp3 = 1 - pchisq(lr3$chisq, length(lr3$n) - 1)
		}else{
			lrp3 = 1
		}
		
		lrps = c(lrps,lrp3)
		thrs = c(thrs,sortedex1[j])
		
		cp3 = coxph(Surv(osd,evos)~ex1)
		scp3 = summary(cp3)
		cpc3 = scp3$coef[2]
		cpcs = c(cpcs,cpc3)
	}
	
	mycol = rep("gray",length(lrps))
	mycol[which(lrps<0.05 & cpcs>1)] = "red"
	mycol[which(lrps<0.05 & cpcs<1)] = "blue"
	mycol[which(lrps==0)] = "green"
	
	lrps[which(lrps==0)] = 10^(-16)
	
	plot(thrs,-log10(lrps),pch=21,bg=mycol,type="b",xlab="Expression Threshold",ylab="-log10(p-value",main=paste0("Gene: ",gene),ylim=c(0,max(-log10(lrps),na.rm=TRUE)),cex=1.6)
	abline(h=-log10(0.05),col="gray",v=quantile(ex1)[2:4])
	return(lrmc_score(lrps))
}


usat = function(x,gene,osd,evos){
	ex1 = as.numeric(x)
	
	cp = coxph(Surv(osd,evos)~ex1)
	scp = summary(cp)
	cpp = scp$coef[5]
	cpc = scp$coef[2]
	cpl = scp$conf.int[3]
	cpu = scp$conf.int[4]
	
	ms = maxstat.test(Surv(osd,evos)~ex1, data=data.frame(osd=c(),evos=c(),ex1=c()), smethod="LogRank", pmethod="Lau94")
	msp = ms$p.value
	msc = ms$estimate
	
	ex2 = ex1
	ex2[which(ex1<=msc)] = 0
	ex2[which(ex1>msc)] = 1
	
	if(length(unique(ex2))>1){
		lr = survdiff(Surv(osd,evos)~ex2)
		lrp = 1 - pchisq(lr$chisq, length(lr$n) - 1)
	}else{
		lrp = 1
	}
	
	km = survfit(Surv(osd,evos)~ex2)
	
	sig = 0
	lrmcs = 0
	if(cpp<0.05 && msp<0.05 && lrp<0.05){
		par(mfrow=c(1,2))
		plot(km,col=c("blue","red"),mark.time=TRUE,lwd=2,xlab="Overall Survival After Diagnosis (months)",ylab="Survival",main=paste0("Training Set - Gene: ",gene,"\nLog-Rank p-value: ",lrp,"\nThrehold: ",msc))
		legend("bottomleft",lwd=2,col=c("blue","red"),legend=c("Expression Below Threshold","Expression Above Threshold"))
		sig = 1
		
		lrmcs = lrmc(ex1,gene,osd,evos)
	}
	return(c(gene,cpc,cpl,cpu,cpp,msc,msp,lrp,sig,lrmcs))
}

setwd("F:/Murat RNA seq Poland/Heatmap")
d = read.delim(file=paste0(inputfilename,".txt"),sep="\t",header=TRUE,row.names=1)
cd = read.delim(file=paste0(clinicaldatafilename,".txt"),sep="\t",header=TRUE)

cn0 = colnames(d)
cn0 = cn0[-length(cn0)]
posbc0 = match(cn0,cd[,1])
posbc = posbc0[which(!is.na(posbc0))]

d1 = d[,c(posbc,ncol(d))]

os = cd$Survival
evos = cd$Status

library(survival)
library(maxstat)

mat1 = c()

pdf(paste0(inputfilename," KM and LRMC Graphs.pdf"),width=16,height=8)

for(i in 1:nrow(d1)){
	cat(paste0(i,"/",nrow(d1),"\r"))
	flush.console()
	gene = paste0(rownames(d1)[i],"-",d1[i,ncol(d1)])
	x = as.numeric(d1[i,-ncol(d1)])
	
	mat1 = rbind(mat1,usat(x,gene,os,evos))
}
dev.off()

colnames(mat1) = c("Gene","CoxPH Coef","Coef Lower 95% CI","Coef Upper 95% CI","CoxPH P-value","Maxstat Cut-off","Maxstat P-value","Log-rank P-value","Candidate","LRMC Score")

write.table(file=paste0(inputfilename," Statistical Test Results.txt"),mat1,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

