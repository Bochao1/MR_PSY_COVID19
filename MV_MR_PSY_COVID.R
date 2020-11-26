#----------------------------multivariable MR
gsmr_data_plot <- read_gsmr_data('Uni_MR_PGC_COVID.eff_plot.gz')
gsmr_res=as.data.frame(gsmr_data_plot$snp_effect)
head(gsmr_res)
dim(gsmr_res)
pheno<- gsmr_data_plot$pheno[3:length(gsmr_data_plot$pheno)]
betas<- paste('beta',pheno,sep = '_'); ses<-paste('se',pheno,sep = '_')
exp.id<- as.character(gsmr_res_c$Exposure[!duplicated(gsmr_res_c$Exposure)])
out.id<- as.character(gsmr_res_c$Outcome[!duplicated(gsmr_res_c$Outcome)])
colnames(gsmr_res)[1:(length(gsmr_data_plot$pheno)+2)]<- c('SNP',pheno,'A1','A2','FREQ')
colnames(gsmr_res)[seq((length(gsmr_data_plot$pheno)+3), ncol(gsmr_res),by=2)]<- c(betas)
colnames(gsmr_res)[seq((length(gsmr_data_plot$pheno)+4), ncol(gsmr_res),by=2)]<- c(ses)
colnames(gsmr_res)

###forward MVMR 

exp.id.sub<-c("SCZ","BIP","AD") 
betas.exp<- paste('beta',exp.id.sub,sep = '_'); ses.exp<-paste('se',exp.id.sub,sep = '_')
betas.out<- paste('beta',out.id,sep = '_'); ses.out<-paste('se',out.id,sep = '_')
dim(gsmr_res)

gsmr_res_mv<-subset(gsmr_res, gsmr_res$AD!= "0000000" | gsmr_res$BIP!= "0000000"  | gsmr_res$MDD_23me!= "0000000"  |gsmr_res$SCZ!= "0000000")

gsmr_res_mv_prune<- clump_data(gsmr_res)

gsmr_res_mv_prune <-subset(gsmr_res_mv_prune, select = c(exp.id.sub, betas.exp,ses.exp, betas.out,ses.out ) )
head(gsmr_res_mv)

for (var in 1: length(out.id)){
  by= as.numeric(paste(gsmr_res[,betas.out[var]]))  ;byse= as.numeric(paste(gsmr_res[,ses.out[var]]))
  bx=as.matrix(bx); bxse=as.matrix(bxse)
  
  MR_input=MendelianRandomization::mr_mvinput(bx= bx, bxse=bxse,
                                              by=by ,byse=  byse)
  
  resultsivw=MendelianRandomization::mr_mvivw(MR_input,distribution="t-dist",model="fixed")
  resultsivw
  
  resultsegger=MendelianRandomization::mr_mvegger(MR_input, orientate = 1,distribution="t-dist")
  resultsegger
  }

