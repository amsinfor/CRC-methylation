#- for the study 'Robust performance of the methylated NTMT1 and MAP3K14-AS1 dual-target test for colorectal cancer detection in plasma by using sense-antisense and dual-MGB probe technique'
library(parallel)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(pROC);
library(glmnet)
library(ggpubr);
library(cowplot)
library(reshape2)
library(ggsci)
library(ggplotify);
library(superheat)

#-----------functions:
map_HumanMethylation450_annotats<-function(expd,db){
	as.character(expd$DATA)->probes
	as.character(db$UCSC_RefGene_Name)->gnames;
	db$Name->names(gnames);
	gnames[probes]->expd_probes_symbol;
	#--------
	which(expd_probes_symbol=="")->na_index;
	probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=probes,"gName"=expd_probes_symbol,expd[,-1],stringsAsFactors=F)->expd_t;
	which(is.na(expd_t$gName))->na_index;
	expd_t$Probe[na_index]->expd_t$gName[na_index];
	#-------
	grep(";",expd_t$gName)->gName_index;
	lapply(gName_index,function(x){
		unlist(strsplit(as.character(expd_t$gName[x]),split=";"))->x_split;
		unique(x_split)->x_split;
		rep(x,length(x_split));
	})->gName_index_res;
	unlist(gName_index_res)->gName_index_res;
	#-----
	lapply(gName_index,function(x){
		unlist(strsplit(as.character(expd_t$gName[x]),split=";"))->x_split;
		unique(x_split);
	})->gName_split_res;
	unlist(gName_split_res)->gName_split_res;
	#-----
	expd_t[-gName_index,]->expd_t1;
	expd_t[gName_index_res,]->expd_t2;
	gName_split_res->expd_t2$gName;
	rbind(expd_t1,expd_t2,stringsAsFactors=F)->expd_t;
	#----------#------------na
	return(expd_t);
}
#---------minfi: 450k of idat
process_TCGA_pancancer_data<-function(expd,clinical_d){
	colnames(expd)[-c(1,2)]->expd_samples;
	unlist(lapply(expd_samples,function(x){
		unlist(strsplit(x,split="-"))->x_split;
		paste(x_split[1:3],collapse="-")->x_sample;
		x_split[4]->x_sample_code;
		c(x_sample,x_sample_code);
	}))->expd_samples;
	data.frame("A0_Samples"=expd_samples[seq(1,length(expd_samples),2)],"SampleTypeCode"=expd_samples[seq(2,length(expd_samples),2)],stringsAsFactors=F)->expd_sample_df;
	#-----------------diff samples 
	setdiff(expd_sample_df$A0_Samples,clinical_d$A0_Samples)->diff_samples;
	#----------------------shared samples:
	clinical_d$A0_Samples->rownames(clinical_d);
	clinical_d[as.character(expd_sample_df$A0_Samples),]->new_clinical_d;
	expd_sample_df$A0_Samples->new_clinical_d$SampleID;
	expd_sample_df$SampleTypeCode->new_clinical_d$SampleTypeCode;
	paste(new_clinical_d$A0_Samples,new_clinical_d$SampleTypeCode,sep="-")->new_clinical_d$A0_Samples;
	unique(new_clinical_d)->new_clinical_d;
	#------add diff samples :
	expd_sample_df[which(expd_sample_df$A0_Samples%in%diff_samples),]->diff_samples;
	for(xi in 1:nrow(diff_samples)){
		diff_samples$A0_Samples[xi]->ds;
		diff_samples$SampleTypeCode[xi]->ds_split_code;
		paste(ds,ds_split_code,sep="-")->ds;
		c(ds,rep(NA,ncol(new_clinical_d)-2),ds_split_code)->ds_res;
		rbind(new_clinical_d,ds_res)->new_clinical_d;
	}
	##---------------
	new("SampleObj",DATA=expd,CliniInfo=new_clinical_d)->new_clinical_d_objRef;
	return(new_clinical_d_objRef);
}
change_values<-function(myd,column,values){
	which(is.na(myd[,column]))->na_index;
	if(length(na_index)>0){
		myd[na_index,]->myd_na;
		as.character(myd_na[,column])->myd_na[,column];
		"Un"->myd_na[is.na(myd_na[,column]),column];
	}else{
		NULL->myd_na;
	}
	myd[!is.na(myd[,column]),]->myd
	colnames(myd)[column]->column_name;
	table(myd[,column])->N_stage.table;
	as.factor(names(N_stage.table))->N_stage.names;
	data.frame(column_name=N_stage.names,"Value"=values)->N_stage.df;
	c()->tmp.value;
	for(i in 1:nrow(myd)){
		for(j in 1:nrow(N_stage.df)){
			if(myd[i,column]==N_stage.df[j,1]){
				c(tmp.value,as.character(N_stage.df[j,2]))->tmp.value;
			}
		}
	}
	tmp.value->myd[,column];
	if(!is.null(myd_na)){
		rbind(myd,myd_na)->myd;
	}
	return(myd);
}
prepare_series_matrix_data<-function(infile){
	readLines(gzfile(infile,'r'))->infile_lines;
	grep("series_matrix_table_begin",infile_lines)->skip_lines;
	read.table(infile,skip=skip_lines,sep="\t",header=T,stringsAsFactors=F,nrows=length(infile_lines)-skip_lines-2)->myd;
	apply(myd[,-1],2,as.numeric)->myd_v;
	data.frame("DATA"=myd[,1],myd_v)->myd;
	return(myd);
}
process_series_matrix<-function(series_file,infors){
	readLines(gzfile(series_file,open='r'))->myd_info;
	c()->res;
	c()->res_colnames;
	for(s in infors){
		myd_info[grep(s,myd_info)]->s_lines;
		for(sl in s_lines){
			unlist(strsplit(sl,split="\t"))->s_lines.split;
			gsub("\"","",s_lines.split)->s_lines.split;
			if(grepl(": ",sl)){	
				unlist(strsplit(s_lines.split[2],split=": "))[1]->s;
				unlist(lapply(s_lines.split,function(x){
					unlist(strsplit(x,split=": "))[2]
				}))->x_res;
				c(res,x_res[-1])->res;
			}else if(grepl("ftp",sl)){
				unlist(lapply(s_lines.split,function(x){
					gsub("\\.gz","",basename(x));
				}))->x_basenames;
				c(res,x_basenames[-1])->res;
			}else{
				c(res,s_lines.split[-1])->res;
			}
			c(res_colnames,s)->res_colnames;
		}
	}
	matrix(res,ncol=length(res_colnames),byrow=F)->res.matrix;
	res_colnames->colnames(res.matrix);
	data.frame(res.matrix,stringsAsFactors=F)->res.matrix
	return(res.matrix);
}
do_ttest_diff<-function(expd,groups,type,group_order=NULL){
	#method: t test for two groups;
	if(class(expd)!="matrix"){
		print("input data is not matrix!");
		return(NULL);
	}
	names(table(groups$Condition))->condition_table;
	if(!is.null(group_order)){
		intersect(group_order,condition_table)->condition_table;
	}
	detectCores()->no_cores;
	makeCluster(no_cores-1)->c1;
	c()->test_results;
	#-----
	as.character(groups[groups$Condition==condition_table[1],1])->group1_samples;
	as.character(groups[groups$Condition==condition_table[2],1])->group2_samples;
	clusterExport(c1,c("expd","group1_samples","group2_samples"),envir=environment());#
	if(type!="pair"){
		parSapply(c1,rownames(expd),function(g){	
			
				expd[g,group1_samples]->group1_values;
				expd[g,group2_samples]->group2_values;
				mean(group1_values,na.rm=T)->group1_mean;
				mean(group2_values,na.rm=T)->group2_mean;
				#---initialize
				1->group_p.value;
				1->group_statistic;
				tryCatch({
					t.test(group1_values,group2_values,alternative="two.sided", var.equal=FALSE)->group1_group2_test;
					group1_group2_test$p.value->group_p.value;
					group1_group2_test$statistic->group_statistic;
				},error=function(e){
					"error";
				})
				if(group1_mean==0){
					0.001->group1_mean;
				}
				#-------
				c(group2_mean/group1_mean,group2_mean,group_statistic,group_p.value,group1_mean)->tmp_res;
				tmp_res;
		})->test_results;
	}else if(type == "pair"){
		parSapply(c1,rownames(expd),function(g){
				expd[g,group1_samples]->group1_values;
				expd[g,group2_samples]->group2_values;
				mean(group1_values,na.rm=T)->group1_mean;
				mean(group2_values,na.rm=T)->group2_mean;
				#---initialize 
				1->group_p.value;
				1->group_statistic;
				tryCatch({
					t.test(group1_values,group2_values,alternative="two.sided", paired=T)->group1_group2_test;
					group1_group2_test$p.value->group_p.value;
					group1_group2_test$statistic->group_statistic;
				},error=function(e){
					"error";
				})
				#----
				if(group1_mean==0){
					0.001->group1_mean;
				}
				#-------
				c(group2_mean/group1_mean,group2_mean,group_statistic,group_p.value,group1_mean)->tmp_res;
				tmp_res;
			
		})->test_results;
		
	}
	unlist(test_results)->test_results;
	matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
	c("logFC","AveExpr","t","P.value","B")->colnames(test_results_matrix);
	#------------
	stopCluster(c1);
	data.frame("TransID"=rownames(expd),test_results_matrix,stringsAsFactors=F)->res;
	unlist(lapply(res$TransID,function(xi){unlist(strsplit(xi,split="\\|"))[1]}))->res$gName;
	unlist(lapply(res$TransID,function(xi){unlist(strsplit(xi,split="\\|"))[2]}))->res$Probe;
	res[!is.na(res$P.value),]->res;
	res[order(res$P.value),]->res;
	p.adjust(res$P.value)->res$FDR;
	return(res);

}
change_GEObjRef_to_expd<-function(objRef,select_f="gName",genes=NULL){
	objRef$DATA->objRef_dat;
	objRef$CliniInfo->objRef_info;
	if(!is.null(genes)){
		subset(objRef_dat,gName%in%genes)->objRef_dat;
	}
	#---
	t(objRef_dat[,-c(1,2)])->dat_;
	if(select_f=="gName"){
		objRef_dat$gName->colnames(dat_);
	}else{
		objRef_dat$Probe->colnames(dat_);
	}
	data.frame("SampleID"=rownames(dat_),dat_)->dat_;
	merge(objRef_info,dat_,by.x="A0_Samples",by.y="SampleID")->res;
	return(res);
}
#--#------------------------ROC curve analysis
add_Y_labels<-function(expd,f,f_labels){
	if(length(f_labels)>1){
		0:(length(f_labels)-1)->f_values;
		f_labels->names(f_values);
		#---
		unlist(lapply(f_labels,function(x){
			which(expd[,f]==x)->x_index;
			x_index;
		}))->keep_rows;
		#-----
		expd[keep_rows,]->expd;
		f_values[as.character(expd[,f])]->expd$Y_label;
	}else if(length(f_labels)==1){
		which(expd[,f]==f_labels[1])->f_index;
		0->expd$Y_label;
		1->expd$Y_label[f_index];
	}else{
		NULL->expd$Y_label;
	}
	return(expd);
}
draw_multiClass_roc<-function(dat,Y_label,variables,myd_colors,gtitle){
	variables[1]->v1;
	as.formula(paste(Y_label,v1,sep="~"))->v1_formula;
	glm(v1_formula, data = dat, family = "binomial")->model_glm
	predict(model_glm, newdata = dat, type = "response")->v1.test_prob
	pROC::roc(dat[,Y_label]~v1.test_prob,plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->v1.test_roc#
	plot(v1.test_roc,main=gtitle, col = myd_colors[1], print.thres = "best", print.auc = T,ci.type="shape");
	#--------
	2->i;
	for(vx in variables[-1]){
		as.formula(paste(Y_label,vx,sep="~"))->vx_formula;
		glm(vx_formula, data = dat, family = "binomial")->model_glm
		predict(model_glm, newdata = dat, type = "response")->vx.test_prob
		pROC::roc(dat[,Y_label]~vx.test_prob, plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->vx.test_roc
		plot(vx.test_roc,main=gtitle, col = myd_colors[i],print.thres = "best", print.auc = T,ci.type="shape",add=T);
		i+1->i;
	}
	legend("bottomright",legend=variables,lty=1,lwd=1,col=myd_colors[1:length(variables)]);
}
draw_multiClass_roc_v2<-function(dat,Y_label,variables,myd_colors,gtitle){
	c()->v_best
	#------------------
	variables[1]->v1;
	as.formula(paste(Y_label,v1,sep="~"))->v1_formula;
	glm(v1_formula, data = dat, family = "binomial")->model_glm
	predict(model_glm, newdata = dat, type = "response")->v1.test_prob
	pROC::roc(dat[,Y_label]~v1.test_prob,plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->v1.test_roc#
	plot(v1.test_roc,main=gtitle, col = myd_colors[1], print.thres = "best", print.auc = T,ci.type="shape");
	#----summary 
	which.max(v1.test_roc$sensitivities+v1.test_roc$specificities)[1]->v_max;
	c(v_best,round(v1.test_roc$sensitivities[v_max],2),round(v1.test_roc$specificities[v_max],2),round(v1.test_roc$auc,2),round(v1.test_roc$ci[c(1,3)],2),paste(round(v1.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
	#--------
	2->i;
	for(vx in variables[-1]){
		as.formula(paste(Y_label,vx,sep="~"))->vx_formula;
		glm(vx_formula, data = dat, family = "binomial")->model_glm
		predict(model_glm, newdata = dat, type = "response")->vx.test_prob
		pROC::roc(dat[,Y_label]~vx.test_prob, plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->vx.test_roc
		plot(vx.test_roc,main=gtitle, col = myd_colors[i],print.thres = "best", print.auc = T,ci.type="shape",add=T);
		i+1->i;
		#-----summary 
		which.max(vx.test_roc$sensitivities+vx.test_roc$specificities)[1]->vx_max;
		c(v_best,round(vx.test_roc$sensitivities[vx_max],2),round(vx.test_roc$specificities[vx_max],2),round(vx.test_roc$auc,2),round(vx.test_roc$ci[c(1,3)],2),paste(round(vx.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
	}
	legend("bottomright",legend=variables,lty=1,lwd=1,col=myd_colors[1:length(variables)]);
	#---
	matrix(v_best,nrow=length(variables),byrow=T)->v_best;
	c("best_sensi","best_spec","auc","auc95cil","auc95ciu","auc_ci")->colnames(v_best);
	data.frame("Feature"=variables,v_best,stringsAsFactors=F)->v_best;
	return(v_best);
}
#---------sensitivity/specificity+95%CI：
draw_multiClass_roc_v3<-function(dat,Y_label,variables,myd_colors,gtitle,ci_value=0.95){
	#---coords:
	c("sensitivity","specificity","ppv","npv","accuracy")->coords_ret;
	#-------------------
	c()->v_best
	c()->v_ci;
	#------------------
	variables[1]->v1;
	as.formula(paste(Y_label,v1,sep="~"))->v1_formula;
	glm(v1_formula, data = dat, family = "binomial")->model_glm
	predict(model_glm, newdata = dat, type = "response")->v1.test_prob
	pROC::roc(dat[,Y_label]~v1.test_prob,plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->v1.test_roc#
	plot(v1.test_roc,main=gtitle, col = myd_colors[1], print.thres = "best", print.auc = T,ci.type="shape");
	#----summary: youden index
	which.max(v1.test_roc$sensitivities+v1.test_roc$specificities)[1]->v_max;
	#----coords:
	coords(v1.test_roc,ret=coords_ret)->v1_coords;
	c(v_best,round(as.numeric(v1_coords[v_max,]),2),round(v1.test_roc$auc,2),round(v1.test_roc$ci[c(1,3)],2),paste(round(v1.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
	#---for sensi/speci ci 95%:
	v1_coords$sensitivity[v_max]/100->v_sensi_max;
	v1_coords$specificity[v_max]/100->v_speci_max;
	length(which(dat$Y_label==1))->v_sensi_size;
	length(which(dat$Y_label==0))->v_speci_size;
	c(v_sensi_max-qnorm((1+ci_value)/2)*sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size),v_sensi_max+qnorm((1+ci_value)/2)*sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size))->v1_sensi_ci;
	c(v_speci_max-qnorm((1+ci_value)/2)*sqrt(v_speci_max*(1-v_speci_max)/v_speci_size),v_speci_max+qnorm((1+ci_value)/2)*sqrt(v_speci_max*(1-v_speci_max)/v_speci_size))->v1_speci_ci;
	#---ppv, npv
	v1_coords$ppv[v_max]/100->v_ppv_max;
	v1_coords$npv[v_max]/100->v_npv_max;
	length(which(dat[,v1]==1))->v_ppv_size;
	length(which(dat[,v1]==0))->v_npv_size;
	c(v_ppv_max-qnorm((1+ci_value)/2)*sqrt(v_ppv_max*(1-v_ppv_max)/v_ppv_size),v_ppv_max+qnorm((1+ci_value)/2)*sqrt(v_ppv_max*(1-v_ppv_max)/v_ppv_size))->v1_ppv_ci;
	c(v_npv_max-qnorm((1+ci_value)/2)*sqrt(v_npv_max*(1-v_npv_max)/v_npv_size),v_npv_max+qnorm((1+ci_value)/2)*sqrt(v_npv_max*(1-v_npv_max)/v_npv_size))->v1_npv_ci;
	#--accuracy
	v1_coords$accuracy[v_max]/100->v_accu;
	v_ppv_size+v_npv_size->v_size;
	c(v_accu-qnorm((1+ci_value)/2)*sqrt(v_accu*(1-v_accu)/v_size),v_accu+qnorm((1+ci_value)/2)*sqrt(v_accu*(1-v_accu)/v_size))->v1_accu_ci;
	c(v_ci,v1_sensi_ci,v1_speci_ci,v1_ppv_ci,v1_npv_ci,v1_accu_ci)->v_ci;
	#print(sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size));flush.console();
	#--------
	2->i;
	for(vx in variables[-1]){
		as.formula(paste(Y_label,vx,sep="~"))->vx_formula;
		glm(vx_formula, data = dat, family = "binomial")->model_glm
		predict(model_glm, newdata = dat, type = "response")->vx.test_prob
		pROC::roc(dat[,Y_label]~vx.test_prob, plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->vx.test_roc
		plot(vx.test_roc,main=gtitle, col = myd_colors[i],print.thres = "best", print.auc = T,ci.type="shape",add=T);
		i+1->i;
		#-----summary 
		which.max(vx.test_roc$sensitivities+vx.test_roc$specificities)[1]->vx_max;
		#----coords:
		coords(vx.test_roc,ret=coords_ret)->vx_coords;
		c(v_best,round(as.numeric(vx_coords[vx_max,]),2),round(vx.test_roc$auc,2),round(vx.test_roc$ci[c(1,3)],2),paste(round(vx.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
		#---sens, speci: ci95%
		vx_coords$sensitivity[vx_max]/100->vx_sensi_max;
		vx_coords$specificity[vx_max]/100->vx_speci_max;
		c(vx_sensi_max-qnorm((1+ci_value)/2)*sqrt(vx_sensi_max*(1-vx_sensi_max)/v_sensi_size),vx_sensi_max+qnorm((1+ci_value)/2)*sqrt(vx_sensi_max*(1-vx_sensi_max)/v_sensi_size))->vx_sensi_ci;
		c(vx_speci_max-qnorm((1+ci_value)/2)*sqrt(vx_speci_max*(1-vx_speci_max)/v_speci_size),vx_speci_max+qnorm((1+ci_value)/2)*sqrt(vx_speci_max*(1-vx_speci_max)/v_speci_size))->vx_speci_ci;
		#---ppv, npv
		vx_coords$ppv[v_max]/100->vx_ppv_max;
		vx_coords$npv[v_max]/100->vx_npv_max;
		length(which(dat[,vx]==1))->vx_ppv_size;
		length(which(dat[,vx]==0))->vx_npv_size;
		c(vx_ppv_max-qnorm((1+ci_value)/2)*sqrt(vx_ppv_max*(1-vx_ppv_max)/vx_ppv_size),vx_ppv_max+qnorm((1+ci_value)/2)*sqrt(vx_ppv_max*(1-vx_ppv_max)/vx_ppv_size))->vx_ppv_ci;
		c(vx_npv_max-qnorm((1+ci_value)/2)*sqrt(vx_npv_max*(1-vx_npv_max)/vx_npv_size),vx_npv_max+qnorm((1+ci_value)/2)*sqrt(vx_npv_max*(1-vx_npv_max)/vx_npv_size))->vx_npv_ci;
		#--accuracy
		vx_coords$accuracy[v_max]/100->vx_accu;
		vx_ppv_size+vx_npv_size->vx_size;
		c(vx_accu-qnorm((1+ci_value)/2)*sqrt(vx_accu*(1-vx_accu)/vx_size),vx_accu+qnorm((1+ci_value)/2)*sqrt(vx_accu*(1-vx_accu)/vx_size))->vx_accu_ci;
		c(v_ci,vx_sensi_ci,vx_speci_ci,vx_ppv_ci,vx_npv_ci,vx_accu_ci)->v_ci;
	}
	legend("bottomright",legend=variables,lty=1,lwd=1,col=myd_colors[1:length(variables)]);
	#---
	matrix(v_best,nrow=length(variables),byrow=T)->v_best;
	c("best_sensi","best_spec","ppv","npv","accuracy","auc","auc95cil","auc95ciu","auc_ci")->colnames(v_best);
	matrix(v_ci,nrow=length(variables),byrow=T)->v_ci;
	c("se.low","se.high","sp.low","sp.high","ppv.low","ppv.high","npv.low","npv.high","accuracy.low","accuracy.high")->colnames(v_ci);
	data.frame("Feature"=variables,v_best,v_ci,stringsAsFactors=F)->v_best;
	return(v_best);
}
do_logistic_fit<-function(dat,Y,X,group=NULL){
	as.formula(paste(Y,X,sep="~"))->tmp_formula;
	glm(tmp_formula, data = dat, family = "binomial")->model_glm;#,nlambda=100
	#----
	predict(model_glm,type = "response")->model_glm_pred#for glm type=c("link", "response", "terms")
	model_glm_pred->dat$PredictedProb
	roc(response=dat[,Y], predictor=dat$PredictedProb,ci=T,ci.method="boot",boot.n=100,quiet=T)->test_roc;#,smooth=T
	if(is.null(group)){
		paste(Y,X,sep="~")->group;
	}
	plot(test_roc,main=group, col = "blue", print.thres = "best", print.auc = T,ci=TRUE, ci.type="shape");
	#---add counts :
	table(dat[,Y])->Y_numbers;
	paste(Y_numbers,collapse="/")->Y_numbers;
	text(x=0.4,y=0.4,labels=paste("N/T",Y_numbers,sep=":"));
	#---
	all.values <- c("threshold","specificity","sensitivity","accuracy","tn","tp","fn","fp","npv","ppv","1-specificity","1-sensitivity","1-npv","1-ppv")
	t(coords(test_roc, "all", ret = all.values))->res;
	return(list("model"=model_glm,"roc"=test_roc,"Score"=model_glm_pred,"Coords"=res));
}
do_logistic_predict<-function(dat,train_model){
	as.character(train_model$formula)->x_y_formula;
	x_y_formula[2]->Y;
	x_y_formula[3]->X;
	#----
	predict(train_model,type = "response",newdata=dat)->model_glm_pred#for glm type=c("link", "response", "terms")
	model_glm_pred->dat$PredictedProb
	#---
	if(length(intersect(Y,colnames(dat)))>0){
		roc(response=dat[,Y], predictor=dat$PredictedProb,ci=T,ci.method="boot",boot.n=100,quiet=T)->test_roc;#,smooth=T
		#---
		paste(Y,X,sep="~")->group;
		plot(test_roc,main=group, col = "blue", print.thres = "best", print.auc = T,ci.type="shape");
		#---add counts :
		table(dat[,Y])->Y_numbers;
		paste(Y_numbers,collapse="/")->Y_numbers;
		text(x=0.4,y=0.4,labels=paste("N/T",Y_numbers,sep=":"));
		#---
		all.values <- c("threshold","specificity","sensitivity","accuracy","tn","tp","fn","fp","npv","ppv","1-specificity","1-sensitivity","1-npv","1-ppv")
		return(list("roc"=test_roc,"coords"=coords(test_roc, "all", ret = all.values),"Score"=model_glm_pred));
	}else{
		return(list("roc"=NULL,"coords"=NULL,"Score"=model_glm_pred));
	}
}
sigmoid<-function(x,coefs){
	1/(1+exp(-x*coefs[2]-coefs[1]))->res;
	return(res);
}
get_x_cutoff<-function(coefs,sigmoid_cutoff){
	-(log(1+1/sigmoid_cutoff)+coefs[1])/coefs[2]->res;
	return(res);
}
#--
calculate_multiP_speci_sensi_auc<-function(expd,features,Y_label){
	intersect(features,colnames(expd))->features;
	which(expd[,Y_label]=="0")->n_index;
	which(expd[,Y_label]=="1")->t_index;
	#------------for normal samples
	if(length(features)==1){
		rep(features,2)->features;
	}else if(length(features)==0){
		return(NULL);
	}
	apply(expd[n_index,features],1,function(x){
		length(which(x==1))->x_positive_count;
		if(x_positive_count>0){
			1;
		}else{
			0;
		}
	})->n_labels;
		#------------for tumor samples 
	apply(expd[t_index,features],1,function(x){
		length(which(x==1))->x_positive_count;
		if(x_positive_count>0){
			1;
		}else{
		0;
			}
	})->t_labels;
	#---------------
	#NA->expd$MethyPredict;
	#n_labels->expd$MethyPredict[n_index];
	#t_labels->expd$MethyPredict[t_index];
	#-------table : for 0
	c(length(which(n_labels==0)),length(which(n_labels==1)),length(which(t_labels==0)),length(which(t_labels==1)))->expd_values;
	matrix(expd_values,nrow=2,byrow=T)->expd_table;
	paste(expd_table[1,],collapse="/")->n_table;#N/P
	paste(expd_table[2,],collapse="/")->t_table;#N/P
	paste(table(expd[,Y_label]),collapse="/")->normal_tumor_count;
	#---
	round(expd_table[1,1]/sum(expd_table[1,]),4)->expd_speci;
	round(expd_table[2,2]/sum(expd_table[2,]),4)->expd_sensi;
	#----------
	c(normal_tumor_count,expd_speci,expd_sensi,n_table,t_table)->res;
	return(res);
}
batch_calculate_multiP<-function(expd,features,Y_label){
	calculate_multiP_speci_sensi_auc(expd,features,Y_label)->expd_res;
	for(fs in features){
		calculate_multiP_speci_sensi_auc(expd,fs,Y_label)->fs_res;
		c(expd_res,fs_res)->expd_res;
	}
	#---------
	matrix(expd_res,ncol=5,byrow=T)->expd_res;
	c("Samples","Specificity","Sensitivity","Normal","Tumor")->colnames(expd_res);
	data.frame("Feature"=c("Combined",features),expd_res,stringsAsFactors=F)->expd_res;
	(as.numeric(expd_res$Specificity)+as.numeric(expd_res$Sensitivity)-1)->expd_res$Youden_index;
	return(expd_res);
}

#---objects---
setRefClass("SampleObj",
	fields=list(DATA="data.frame",CliniInfo="data.frame"),
	methods=list(
		initialize=function(DATA,CliniInfo){
			intersect(colnames(DATA),CliniInfo$A0_Samples)->shared_samples;
			CliniInfo[which(CliniInfo$A0_Samples%in%shared_samples),]->>CliniInfo;
			DATA[,c("gName","Probe",shared_samples)]->>DATA;
		},
		getSYMBOL=function(probes){
			unlist(lapply(probes,function(i){
				which(DATA$Probe==i)->i_index
			}))->probes_index;
			return(as.character(DATA$gName[probes_index]))
		},
		getConditionGroup=function(f,groups){
			c()->groups_index;
			for(g in groups){
				which(CliniInfo[,f]==g)->g_index;
				c(groups_index,g_index)->groups_index;
			}
			CliniInfo[groups_index,c("A0_Samples",f)]->condition_groups;
			c("SampleID","Condition")->colnames(condition_groups);
			return(condition_groups);
		},
		getGroupMatrix=function(f,groups){
			getConditionGroup(f,groups)->condition_groups;
			#----
			as.matrix(DATA[,as.character(condition_groups$SampleID)])->dat_matrix;
			paste(DATA$gName,DATA$Probe,sep="|")->rownames(dat_matrix);
			return(dat_matrix);
		},
		getDataMatrix=function(genes,select_f="gene"){
			as.matrix(DATA[,-c(1,2)])->dat_matrix;
			paste(DATA$gName,DATA$Probe,sep="|")->rownames(dat_matrix);
			if(select_f=="gene"){
				which(DATA$gName%in%genes)->g_index;
			}else if(select_f=="probe"){
				which(DATA$Probe%in%genes)->g_index;
			}
			dat_matrix[g_index,]->dat_matrix;
			return(dat_matrix);
		},
		getSubSet=function(f,groups){
			c()->groups_index;
			for(g in groups){
				which(CliniInfo[,f]==g)->g_index;
				c(groups_index,g_index)->groups_index;
			}
			CliniInfo[groups_index,]->sub_CliniInfo;
			as.character(sub_CliniInfo[,f])->sub_CliniInfo[,f];
			#---
			colnames(DATA)[1:2]->two_cols;
			DATA[,c(two_cols,as.character(sub_CliniInfo$A0_Samples))]->sub_DATA;
			new("SampleObj",DATA=sub_DATA,CliniInfo=sub_CliniInfo)->sub_SampleObj;
			return(sub_SampleObj);
		},
		addFeatures=function(feature_df){
			merge(CliniInfo,feature_df,by.x="A0_Samples",by.y="A0_Samples",all.x=T)->>CliniInfo;
		},
		getDataSubSet=function(genes,select_f="gene"){
			if(select_f=="gene"){
				which(DATA$gName%in%genes)->g_index;
			}else if(select_f=="probe"){
				which(DATA$Probe%in%genes)->g_index;
			}
			DATA[g_index,]->sub_DATA;
			new("SampleObj",DATA=sub_DATA,CliniInfo=CliniInfo)->sub_SampleObj;
			return(sub_SampleObj);
		}
))->SampleObj;
#--
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
get_palette("npg",10)->myd_colors_npg;
c("METTL11A","C9orf50","LOC100133991")->target_genes;#: METTL11A+C9orf50 => NTMT1（M55）和 LOC100133991 => MAP3K14-AS1（M58）
#-
HumanMethylation450_annotats[which(HumanMethylation450_annotats$CHR==9),]->chr9_annots;
chr9_annots[intersect(which(chr9_annots$MAPINFO>=132371163),which(chr9_annots$MAPINFO<=132398209)),]#NTMT1（M55）
#-
HumanMethylation450_annotats[which(HumanMethylation450_annotats$CHR==17),]->chr17_annots;
chr17_annots[intersect(which(chr17_annots$MAPINFO>=43325292),which(chr17_annots$MAPINFO<=43345997)),]#MAP3K14-AS1（M58）

#---------------------------------------------------------------------
#-----------------TCGA CRC:
process_TCGA_pancancer_data(TCGA_CRC_methy,TCGA_CRC_clinical)->TCGA_CRC_methy_objRef;
TCGA_CRC_methy_objRef$getSubSet("SampleTypeCode",c("11","01"))->TCGA_CRC_methy_objRef;
change_values(TCGA_CRC_methy_objRef$CliniInfo,26,c("Cancer","Normal"))->TCGA_CRC_methy_objRef$CliniInfo;
#---------------------------------------------------------------------
#-----------------GEO CRC: GSE48684
#-for GSE48684
new("SampleObj",DATA=GSE48684_methy,CliniInfo=GSE48684_targets)->GSE48684_methy_objRef;
#---------------------------------------------------------------------
#-----------------GEO CRC: GSE40279/GSE122126
#--for GSE40279: 
new("SampleObj",DATA=GSE40279_methy,CliniInfo=GSE40279_targets)->GSE40279_methy_objRef;
#----------------------------------------------------------------------
#--GSE122126：
new("SampleObj",DATA=GSE122126_GPL21145_methy,CliniInfo=GSE122126_GPL21145_targets)->GSE122126_GPL21145_methy_objRef

######################################################################################################################################################
###############################################################################################################################
#------------------
TCGA_CRC_methy_objRef$getConditionGroup("SampleTypeCode",c("Normal","Cancer"))->TCGA_CRC_methy_group;
TCGA_CRC_methy_objRef$getGroupMatrix("SampleTypeCode",c("Normal","Cancer"))->TCGA_CRC_methy_matrix;
do_ttest_diff(TCGA_CRC_methy_matrix,TCGA_CRC_methy_group,"unpair",c("Normal","Cancer"))->TCGA_CRC_methy_ttest;
#--filter:
TCGA_CRC_methy_ttest$AveExpr-TCGA_CRC_methy_ttest$B->TCGA_CRC_methy_ttest$DeltaBeta;
TCGA_CRC_methy_ttest[which(TCGA_CRC_methy_ttest$DeltaBeta>=0.3),]->TCGA_CRC_methy_ttest_filter;
TCGA_CRC_methy_ttest_filter[TCGA_CRC_methy_ttest_filter$B<0.2,]->TCGA_CRC_methy_ttest_filter;

#--GEO tissue:
GSE48684_methy_objRef$getConditionGroup("SampleType",c("Normal","Cancer"))->GSE48684_CRC_methy_group;
GSE48684_methy_objRef$getGroupMatrix("SampleType",c("Normal","Cancer"))->GSE48684_CRC_methy_matrix;
do_ttest_diff(GSE48684_CRC_methy_matrix,GSE48684_CRC_methy_group,"unpair",c("Normal","Cancer"))->GSE48684_CRC_methy_ttest;
#--filter:
GSE48684_CRC_methy_ttest$AveExpr-GSE48684_CRC_methy_ttest$B->GSE48684_CRC_methy_ttest$DeltaBeta;
GSE48684_CRC_methy_ttest[which(GSE48684_CRC_methy_ttest$DeltaBeta>=0.3),]->GSE48684_CRC_methy_ttest_filter;
GSE48684_CRC_methy_ttest_filter[TCGA_CRC_methy_ttest_filter$B<0.2,]->GSE48684_CRC_methy_ttest_filter;
#---shared:
intersect(TCGA_CRC_methy_ttest_filter$Probe,GSE48684_CRC_methy_ttest_filter$Probe)->TCGA_GSE48684_shared_diff_probes;

#---blood samples:
GSE40279_methy_objRef$getConditionGroup("SampleType",c("whole blood"))->GSE40279_CRC_methy_group;
GSE40279_methy_objRef$getGroupMatrix("SampleType",c("whole blood"))->GSE40279_CRC_methy_matrix;
lapply(rownames(GSE40279_CRC_methy_matrix),function(rx){unlist(strsplit(rx,split="\\|"))[2]})->rownames(GSE40279_CRC_methy_matrix)
apply(GSE40279_CRC_methy_matrix,1,mean,na.rm=T)->GSE40279_CRC_methy_probe_values;
names(which(GSE40279_CRC_methy_probe_values<0.1))->GSE40279_CRC_methy_filtered_probes;

#############################################################################################
#-----------------------------figure 2： methylation values of "C9orf50|cg14015706"      "LOC100133991|cg08247376"
change_GEObjRef_to_expd(TCGA_CRC_methy_objRef,select_f="probe")->TCGA_CRC_methy_factors;
change_GEObjRef_to_expd(GSE48684_methy_objRef,select_f="probe")->GSE48684_methy_factors;
change_values(TCGA_CRC_methy_factors,7,c("Un","Stage I","Stage I","Stage II","Stage II","Stage II","Stage II","Stage III","Stage III","Stage III","Stage III","Stage IV","Stage IV","Stage IV"))->TCGA_CRC_methy_factors;
"Normal"->TCGA_CRC_methy_factors$AJCC_Stage[which(TCGA_CRC_methy_factors$SampleTypeCode=="Normal")]


#--sort by normal vs cancer
sort_by_factor_order<-function(inputd,f,f_orders,cg_id=NULL){	
	if(!is.null(cg_id)){
		inputd[order(inputd[,cg_id]),]->inputd;
	}
	lapply(f_orders,function(fx){
		which(inputd[,f]==fx)
	})->f_index;
	unlist(f_index)->f_index;
	inputd[f_index,]->inputd;
	1:nrow(inputd)->inputd$RowIndex;
	#---
	if(!is.null(cg_id)){
		inputd[,c(f,cg_id,"RowIndex")]->inputd;
		"Value"->names(inputd)[2];
		cg_id->inputd$CpG;
	}
	return(inputd)
}
sort_by_factor_order(TCGA_CRC_methy_factors,"SampleTypeCode",c("Normal","Cancer"),"cg14015706")->TCGA_CRC_cg14015706_methy;
sort_by_factor_order(TCGA_CRC_methy_factors,"SampleTypeCode",c("Normal","Cancer"),"cg08247376")->TCGA_CRC_cg08247376_methy;
rbind(TCGA_CRC_cg14015706_methy,TCGA_CRC_cg08247376_methy)->test_;
ggplot(test_,aes(x=RowIndex))+geom_area(aes(y=Value,fill=SampleTypeCode))+facet_wrap(~CpG,nrow=2)+scale_y_continuous(breaks=seq(0,1,0.2))+scale_fill_manual(values=pal_jco()(2))+theme_minimal()+scale_fill_manual(values=c("blue","red"))+theme(legend.position="top")->p1;

#-GSE48684: tissue validation
melt(GSE48684_methy_factors,value.name="Methylation",id.vars=c("SampleType"),as.is=T,measure.vars=c("cg14015706","cg08247376"))->ggboxplot_data;
factor(ggboxplot_data$SampleType,levels=c("Normal","Cancer"))->ggboxplot_data$SampleType;
ggboxplot(ggboxplot_data, x = "variable", y = "Methylation",col="SampleType",bxp.errorbar=T,add="jitter")+scale_y_continuous(breaks=seq(0,1,0.2))+theme_minimal()+scale_color_manual(values=c("blue","red"))->p2;

#--GSE40279: blood samples
library(ComplexHeatmap)
library(circlize)
colorRamp2(c(0,0.5,1), c("blue", "white", "red"))->fill_colors
Heatmap(GSE40279_CRC_methy_matrix[c("cg14015706","cg08247376"),],col=fill_colors,show_column_names=F,show_row_names=F,show_heatmap_legend=F)->p3;
#---
plot_grid(p1,p2,as.grob(p3),nrow=2,rel_widths=c(3,2))
#------------ age vs methylation 
ggscatter(test_,x="Age",y="cg14015706",add="reg.line",size=0.5,conf.int = T,cor.coef = T,xlab="Patient age (year)",ylab="beta value (cg14015706)")+theme_minimal()->p1;
ggscatter(test_,x="Age",y="cg08247376",add="reg.line",size=0.5,conf.int = T,cor.coef = T,xlab="Patient age (year)",ylab="beta value (cg08247376)")+theme_minimal()->p2
plot_grid(p1,p2)

#########################################################################################################################################
#################################################---- GSE122126: cfDNA
#-: #"C9orf50|cg14015706"      "LOC100133991|cg08247376"
#-----------------------------heatmap plot:
c("cg14015706","cg03048083")->target_850k_probes;#,"cg01627847","cg08124910","cg26742995"
GSE122126_GPL21145_methy_objRef$getSubSet("SampleType",c("cfDNA"))->test_subset;
test_subset$getGroupMatrix("DiseaseStatus",c("Colon adenocarcinoma","Normal"))->GSE122126_CRC_matrix;
lapply(rownames(GSE122126_CRC_matrix),function(rx){unlist(strsplit(rx,split="\\|"))[2]})->rownames(GSE122126_CRC_matrix)
pheatmap(GSE122126_CRC_matrix[target_850k_probes,],display_numbers=T,number_color="white")

###################################################################################################################################
##--sanger sequencing for NTMT1（M55）and MAP3K14-AS1（M58）
library(sangerseqR)
ggmsa(M76_methy_pa[c("Ref_seq",target_samples)], color = "Shapely_NT", char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights)
#-------------------------------------------------- 处理M76正常组织的测序结果
ggmsa(M76_unmethy_pa[c("Ref_seq",target_samples)],char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights,color = "Shapely_NT")
#------------------------------------------------------------------------------
#----for NTMT1 sense-strand:
ggmsa(M55_2_methy_pa[c("Ref_seq",target_samples)], color = "Shapely_NT", char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights)
#------------
ggmsa(M55_2_unmethy_pa[c("Ref_seq",target_samples)], color = "Shapely_NT", char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights)
#----for NTMT1 antisense-strand 
ggmsa(M55_0_methy_pa[c("Ref_seq",target_samples)], color = "Shapely_NT", char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights)
#------------
ggmsa(M55_0_unmethy_pa[c("Ref_seq",target_samples)], color = "Shapely_NT", char_width = 0.8, seq_name = TRUE,position_highlight=position_highlights)
###############################################################################################################################
#---- figure 4
library(Superheat);
#--------------------------------------------------------------------------------------------------------------------------------------------
superheat(t(M55_0_CpG_count[,-1]),yt=apply(M55_0_CpG_count[,-1],1,sum,na.rm=T),yt.plot.type="bar",left.label.col = "white",left.label.text.size=2,scale=FALSE,membership.cols=M55_0_CpG_count$SampleType,heat.pal = c("white","black"),grid.hline=F)
#-
superheat(t(M55_2_CpG_count[,-1]),yt=apply(M55_2_CpG_count[,-1],1,sum,na.rm=T),yt.plot.type="bar",left.label.col = "white",left.label.text.size=2,scale=FALSE,membership.cols=M55_2_CpG_count$SampleType,heat.pal = c("white","black"),grid.hline=F)
#-
superheat(t(M76_CpG_count[,-1]),yt=apply(M76_CpG_count[,-1],1,sum,na.rm=T),yt.plot.type="bar",left.label.col = "white",left.label.text.size=2,scale=FALSE,membership.cols=M76_CpG_count$SampleType,heat.pal = c("white","black"),grid.hline=F)




#############################################################################################
#-----------------------------figure 5： 
read.table("training-set Ct.txt",header=T, stringsAsFactors=F, encoding="UTF-8",sep="\t")->CRC_blood_training_ct;
read.table("test-set Ct.txt",header=T, stringsAsFactors=F, encoding="UTF-8",sep="\t")->CRC_blood_test_ct;
CRC_blood_test_ct[-which(CRC_blood_test_ct$ACTB_Ct>38),]->CRC_blood_test_ct;
#--ct values
prepare_graphprim_column(CRC_blood_training_ct,c("M55_Ct"),f="DiseaseStatus")->test_;
write.table(test_,"results/CRC_blood_training_M55_Ct.txt",quote=F,sep="\t",row.names=F);
prepare_graphprim_column(CRC_blood_training_ct,c("M58_Ct"),f="DiseaseStatus")->test_;
#-ROC:
prepare_graphprim_column(CRC_blood_training_ct,c("M55_Ct"),f="SampleType")->test_;
write.table(test_,"results/CRC_blood_training_M55_Ct-ROC.txt",quote=F,sep="\t",row.names=F);
prepare_graphprim_column(CRC_blood_training_ct,c("M58_Ct"),f="SampleType")->test_;

#----
#--M55: CRC vs NonCRC
add_Y_labels(CRC_blood_training_ct,"SampleType",c("NonCRC","Cancer"))->CRC_blood_training_ct;
do_logistic_fit(CRC_blood_training_ct,"Y_label","M55_Ct")->CRC_blood_training_M55_logit;

#---------------
#-- CRC vs healthy
#do_logistic_fit(CRC_blood_training_ct_subset,"Y_label","M58_Ct")->CRC_blood_training_M58_subset_logit;
#--M58: CRC vs NonCRC
do_logistic_fit(CRC_blood_training_ct,"Y_label","M58_Ct")->CRC_blood_training_M58_logit;
#sigmoid(CRC_blood_training_ct$M58_Ct,coef(CRC_blood_training_M58_logit$model))
get_x_cutoff(coef(CRC_blood_training_M55_logit$model),CRC_blood_training_M55_logit$Coords[1,])->M55_Ct_cutoff;
get_x_cutoff(coef(CRC_blood_training_M58_logit$model),CRC_blood_training_M58_logit$Coords[1,])->M58_Ct_cutoff;
#---------------------------combination of two makers ---------------------------------------------------------------------------------
#-----
dichotomy_samples(CRC_blood_training_ct,"M55_Ct",M55_Ct_cutoff,larger_true=F)->CRC_blood_training_ct;
dichotomy_samples(CRC_blood_training_ct,"M58_Ct",M58_Ct_cutoff,larger_true=F)->CRC_blood_training_ct;
#-
batch_calculate_multiP(CRC_blood_training_ct,c("M55_Ct_dichotomy","M58_Ct_dichotomy"),"Y_label")

#----- logistic model 
do_logistic_fit(CRC_blood_training_ct,"Y_label","M55_Ct+M58_Ct")->CRC_blood_training_combin_logit;
CRC_blood_training_combin_logit$Score->CRC_blood_training_ct$M55_M58_combined_score;
#--
1->CRC_blood_training_ct$M55_M58_combined;
0->CRC_blood_training_ct$M55_M58_combined[which((CRC_blood_training_ct$M55_Ct_dichotomy+CRC_blood_training_ct$M58_Ct_dichotomy)==0)];
do_logistic_fit(CRC_blood_training_ct,"Y_label","M55_M58_combined")->CRC_blood_training_combin_logit2;

#############################################################################################################################
#------------------------------ patient age/stage in training cohort
change_values(CRC_blood_training_ct,6,c("Un","I","II","III","IV"))->CRC_blood_training_ct;
subset(CRC_blood_training_ct,DiseaseStatus=="Cancer")->CRC_blood_training_ct_subset;
table(CRC_blood_training_ct_subset$Stage,CRC_blood_training_ct_subset$M55_M58_combined)->test_;
test_/apply(test_,1,sum)
#--
table(CRC_blood_training_ct$DiseaseStatus,CRC_blood_training_ct$M55_M58_combined)->test_;
test_/apply(test_,1,sum)
#--------------------------------#--plot
ggboxplot(subset(CRC_blood_training_ct,DiseaseStatus=="Cancer"),x="Age_group",y="M55_Ct",add=c("jitter","median_q1q3"),order=c("a40","a40_a60","a60"),bxp.errorbar=T,color="Age_group",palette="Set1")+stat_compare_means()+theme_minimal()+theme(legend.position = "none")->p1
ggboxplot(subset(CRC_blood_training_ct,DiseaseStatus=="Cancer"),x="Age_group",y="M58_Ct",add=c("jitter","median_q1q3"),order=c("a40","a40_a60","a60"),bxp.errorbar=T,color="Age_group",palette="Set1")+stat_compare_means()+theme_minimal()+theme(legend.position = "none")->p2
#-
plot_grid(p1,p2)
#------------------------------------------------------------------------------------------
#-----------by sex 
subset(CRC_blood_training_ct,SampleType=="Cancer")->test_;
table(test_$Gender,test_$M55_M58_combined)



#############################################################################################################################
#------------------------------for validation cohort
add_Y_labels(CRC_blood_test_ct,"SampleType",c("NonCRC","Cancer"))->CRC_blood_test_ct;
#----- 
dichotomy_samples(CRC_blood_test_ct,"M55_Ct",M55_Ct_cutoff,larger_true=F)->CRC_blood_test_ct;
dichotomy_samples(CRC_blood_test_ct,"M58_Ct",M58_Ct_cutoff,larger_true=F)->CRC_blood_test_ct;
#-
batch_calculate_multiP(CRC_blood_test_ct,c("M55_Ct_dichotomy","M58_Ct_dichotomy"),"Y_label")
#----- logistic model 
do_logistic_predict(CRC_blood_test_ct,CRC_blood_training_M55_logit$model)->test_;
#############################################################################################################################
#------------------------------ by patient age, stage sex in validation cohort
1->CRC_blood_test_ct$M55_M58_combined;
0->CRC_blood_test_ct$M55_M58_combined[which((CRC_blood_test_ct$M55_Ct_dichotomy+CRC_blood_test_ct$M58_Ct_dichotomy)==0)];
change_values(CRC_blood_test_ct,5,c("Un","I","I","II","III","III","IV"))->CRC_blood_test_ct;
subset(CRC_blood_test_ct,DiseaseStatus=="CRC")->CRC_blood_test_ct_subset;
table(CRC_blood_test_ct_subset$Stage,CRC_blood_test_ct_subset$M55_M58_combined)->test_;
round(test_/apply(test_,1,sum)*100,2)
#------------ for different disease status
#-- CRC, Adenoma, BenignDisease, OtherCancer, OtherDisease, Polyps vs Healthy
table(CRC_blood_test_ct$DiseaseStatus,CRC_blood_test_ct$M55_M58_combined)->test_;
round(test_/apply(test_,1,sum)*100,2)
#----------- 
ggboxplot(subset(CRC_blood_test_ct,DiseaseStatus=="CRC"),x="Age_group",y="M55_Ct",add=c("jitter","median_q1q3"),order=c("a40","a40_a60","a60"),bxp.errorbar=T,color="Age_group",palette="Set1")+stat_compare_means()+theme_minimal()+theme(legend.position = "none")->p1
ggboxplot(subset(CRC_blood_test_ct,DiseaseStatus=="CRC"),x="Age_group",y="M58_Ct",add=c("jitter","median_q1q3"),order=c("a40","a40_a60","a60"),bxp.errorbar=T,color="Age_group",palette="Set1")+stat_compare_means()+theme_minimal()+theme(legend.position = "none")->p2
#-
plot_grid(p1,p2)
#------------------------------------------------------------------------------------------
#----------
subset(CRC_blood_test_ct,SampleType=="Cancer")->test_;
table(test_$Gender,test_$M55_M58_combined)


###################################################################################################################################
#===----- cfDNA copies for training set
ggboxplot(CRC_blood_training_ct,x="Stage",y="M55_copies",yscale="log2",order=c("Un","I","II","III","IV"),bxp.errorbar = T,fill="Stage",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p1
ggboxplot(CRC_blood_training_ct,x="Stage",y="M76_copies",yscale="log2",order=c("Un","I","II","III","IV"),bxp.errorbar = T,fill="Stage",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p2
plot_grid(p1,p2)
#---validation set
exp((CRC_blood_test_ct$M55_Ct-41.46)/(-3.080))*40->CRC_blood_test_ct$M55_copies;
exp((CRC_blood_test_ct$M58_Ct-39.98)/(-3.135))*40->CRC_blood_test_ct$M76_copies;
exp((CRC_blood_test_ct$ACTB_Ct-40.62)/(-3.077))*40->CRC_blood_test_ct$ACTB_copies;
#-
ggboxplot(CRC_blood_test_ct,x="DiseaseStatus",y="M55_copies",yscale="log2",order=c("Healthy","Polyps","BenignDisease","Adenoma","CRC","OtherCancer","OtherDisease"),bxp.errorbar = T,fill="DiseaseStatus",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p1
ggboxplot(CRC_blood_test_ct,x="DiseaseStatus",y="M76_copies",yscale="log2",order=c("Healthy","Polyps","BenignDisease","Adenoma","CRC","OtherCancer","OtherDisease"),bxp.errorbar = T,fill="DiseaseStatus",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p2
ggboxplot(CRC_blood_test_ct,x="Stage",y="M55_copies",yscale="log2",order=c("Un","I","II","III","IV"),bxp.errorbar = T,fill="Stage",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p3
ggboxplot(CRC_blood_test_ct,x="Stage",y="M76_copies",yscale="log2",order=c("Un","I","II","III","IV"),bxp.errorbar = T,fill="Stage",palette="Set1",add=c("jitter","median_q1q3"),legend="")+stat_compare_means()->p4
plot_grid(p1,NULL,p2,NULL,p3,p4,nrow=3,ncol=2)
#------------------------------summary:
draw_multiClass_roc_v3(CRC_blood_test_ct,"Y_label",c("M55_Ct_dichotomy","M58_Ct_dichotomy","M55_M58_combined"),get_palette("aaas",5),"NonCRC vs CRC")->test_;
"CRC vs NonCRC"->test_$Comparisons;
#-
draw_multiClass_roc_v3(subset(CRC_blood_test_ct,DiseaseStatus%in%c("OtherCancer","OtherDisease","Polyps","Adenoma","BenignDisease","CRC")),"Y_label",c("M55_Ct_dichotomy","M58_Ct_dichotomy","M55_M58_combined"),get_palette("aaas",5),"CRC vs interfering")->test_1;
"CRC vs interfering"->test_1$Comparisons;
rbind(test_,test_1)->test_;
#-
draw_multiClass_roc_v3(subset(CRC_blood_test_ct,DiseaseStatus%in%c("Healthy","CRC")),"Y_label",c("M55_Ct_dichotomy","M58_Ct_dichotomy","M55_M58_combined"),get_palette("aaas",5),"CRC vs Healthy")->test_1;
"CRC vs Healthy"->test_1$Comparisons;
rbind(test_,test_1)->test_;
















