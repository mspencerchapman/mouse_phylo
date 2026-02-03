library(stringr)
library(ape)
library(seqinr)
library(data.table)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)

my_working_directory = ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/mouse_phylo/","/lustre/scratch126/casm/team154pc/ms56/mouse_phylo/")
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch126/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
setwd(tree_mut_dir); source("treemut.R");setwd(my_working_directory)


bbprob.calculator <- function(x_i, n_i, alpha_error, beta_error, alpha_mut, beta_mut, prior_mut) {
  # x_i is vector of counts of mutant reads across different mutations
  # n_i is vector of total depth across mutations
  # alpha_error, beta_error are vectors of the parameters of the error distribution across mutations
  # alpha_mut, beta_mut are vectors of the parameters of the mutation distribution across mutations
  # prior_mut is the prior probability that the mutation is present in the sample
  
  log_Pr_D_M0 <- lbeta(x_i + alpha_error, n_i - x_i + beta_error) - lbeta(alpha_error, beta_error)
  log_Pr_D_M1 <- lbeta(x_i + alpha_mut, n_i - x_i + beta_mut) - lbeta(alpha_mut, beta_mut)
  
  log_bayes_factor <- log_Pr_D_M0 - log_Pr_D_M1
  log_odds_prior <- log(prior_mut) - log(1-prior_mut)
  
  post_prob <- 1 / (1 + exp(log_bayes_factor - log_odds_prior))
  return(post_prob)
}


moment.est <- function(x_i, n_i, w_i) {
  # Function to calculate mu_hat and gamma_hat from data and a given set of weights
  # Used in estimation of parameters of beta distribution
  w <- sum(w_i)
  p_i_hat <- x_i / n_i
  p_hat <- sum(w_i * p_i_hat) / w
  S <- sum(w_i * (p_i_hat - p_hat)^2) * (length(n_i)-1) / length(n_i)
  temp_sum <- sum(w_i * (1 - w_i/w) / n_i)
  temp_sum_2 <- sum(w_i * (1 - w_i/w))
  gamma_hat <- (S - p_hat * (1-p_hat) * temp_sum) / (p_hat * (1-p_hat) * (temp_sum_2 - temp_sum))
  if (is.nan(gamma_hat) | gamma_hat < 0) {gamma_hat <- 0}
  return(c(p_hat, gamma_hat))
}

generate_targ_seq_plots=function(samples,
                                 tree,
                                 details_targ,
                                 matrices,
                                 post.prob,
                                 info_type=c("post.prob","cell_frac","log_cell_frac"),
                                 prob_threshold_to_include=0.5, #Probability threshold from the post.prob matrix for plotting
                                 plot_cell_frac=TRUE,
                                 plot_donut=TRUE,
                                 donut_info="cell_frac", #other option is "lineages_lost"
                                 CI=0.8,  #Confidence intervals on the pie chart, default = 80% CI
                                 radius=3.5,  #Radius of the pie charts on the plot
                                 scale_muts_to_branch=FALSE,
                                 colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                                 vaf_lwd=5,
                                 overlay=F,
                                 title=NULL) {
  require(plotrix)
  post.prob.mat <-post.prob[details_targ$mut_ref,samples,drop=F]
  if(info_type=="post.prob") {
    post.prob.mat[post.prob.mat<prob_threshold_to_include] <- 0
    details_targ_full=cbind(details_targ,post.prob.mat)
  } else {
    cell_frac_present=calculate_cell_frac(matrices$NV[details_targ$mut_ref,gsub("_comb","",samples),drop=F],matrices$NR[details_targ$mut_ref,gsub("_comb","",samples),drop=F])
    colnames(cell_frac_present)<-samples
    cell_frac_present[post.prob.mat<prob_threshold_to_include]<-0
    if(info_type=="cell_frac") {
      details_targ_full=cbind(details_targ,cell_frac_present)
    } else if(info_type=="log_cell_frac") {
      log_cell_frac_present=cell_frac_present
      log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])#change all the non 0 vafs to the log of their vaf
      log_cell_frac_present_scaled=log_cell_frac_present
      scale_range=c(0.01,1)
      log_cell_frac_present_scaled[cell_frac_present!=0] = plotrix::rescale(log_cell_frac_present[cell_frac_present!=0],newrange = scale_range) #scale these figures between 0 and 1
      details_targ_full=cbind(details_targ,log_cell_frac_present_scaled)
    }
  }
  colnames(details_targ_full)<-gsub("_comb","",colnames(details_targ_full))
  lims=par("usr")
  sapply(samples, function(sample) {
    sample_stripped=gsub("_comb","",sample)
    if(!overlay){
      tree=plot_tree(tree, cex.label = 0,lwd=0.5,plot_axis = TRUE,default_edge_color="lightgrey")
    }
    lims=par("usr")
    if(is.null(title)){
      text(0,lims[4]-(0.05*lims[4]),paste0(sample,": Mean depth is ",round(mean(matrices$NR[,sample_stripped]),digits = 2)),cex=1,pos=4)
    } else {
      text(0,lims[4]-(0.05*lims[4]),paste0(title,": Mean depth is ",round(mean(matrices$NR[,sample_stripped]),digits = 2)),cex=1,pos=4)
    }
    
    
    
    add_annotation(tree=tree,
                   details=details_targ_full,
                   matrices,
                   annot_function=function(tree,details,matrices,node) {
                     add_var_col(tree,
                                 details,
                                 matrices,
                                 node,
                                 var_field = sample_stripped,
                                 pval_based=FALSE,
                                 lwd = vaf_lwd,
                                 colours=colour.scale,
                                 scale_muts_to_branch=scale_muts_to_branch)
                   }
    )
    if(plot_cell_frac) {
      print("Plotting cell fraction")
      add_annotation_targeted(sample,
                              tree=tree,
                              details=details_targ_full,
                              matrices=matrices,
                              annot_function=function(node,sample,tree,details,matrices,cex=0.6) {
                                node_cell_frac=get_node_cell_frac(node,sample_stripped,tree,details,matrices)
                                node_cell_frac<-min(node_cell_frac,1)
                                info=get_edge_info(tree,details,node)
                                if(!is.na(node_cell_frac) & any(post.prob.mat[info$idx.in.details,sample]>prob_threshold_to_include)) {
                                  text(info$x,info$yb,round(node_cell_frac,digits=3),cex = cex,col="black",font=2)
                                }
                              })
    }
    if(plot_donut) {
      #Detect which nodes to plot donuts for - do it for branches with a mean clean.post.prob of >0.5
      nodes_to_check=unique(tree$edge[,2])[!unique(tree$edge[,2])%in%1:length(tree$tip.label)]
      nodes_to_include=sapply(nodes_to_check,function(node) {
        if(sum(details_targ$node==node)==0) {
          return(NA)
        } else if(mean(post.prob.mat[details_targ$node==node,sample,drop=F])>prob_threshold_to_include){
          return(node)
        }else{
          return(NA)
        }
      })
      nodes_to_include<-nodes_to_include[!is.na(nodes_to_include)]
      
      if(donut_info=="cell_frac") {
        print("Starting to plot the cell fraction pie charts")
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = TRUE)
          #Make the pie chart for plotting on the node
          df2<-data%>%dplyr::select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          plotDonut(info$x,mean(c(info$yb,info$yt)),median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black")
          
          median_only_prop=c(data$median,1-data$median)
          
          #plotDonut(info$x,mean(c(info$yb,info$yt)),radius=radius,col=c( "#08306B" ,"#FFFFFF"),prop=median_only_prop,border="black",plotPie = TRUE)
        }
      } else if(donut_info=="lineages_lost") {
        
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = FALSE)
          #Make the pie chart for plotting on the node
          df2<-data%>%select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          plotDonut(info$x,info$yb,median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black",plotPie = TRUE)
        }
      }
    }
  }
  )
}

#Define function to import allele counter data similarly to cgpVAF NV and NR matrices
import_allele_counter_data=function(allele_counter_dir,suffix=".allelecounts.txt",bed_file=SNV_mut_refs,verbose=F){
  require(dplyr)
  require(stringr)
  output_files=list.files(path=allele_counter_dir,pattern=suffix,full.names = T)
  Samples=unlist(lapply(stringr::str_split(output_files,pattern="/"),function(x) gsub(suffix,"",tail(x,n=1))))
  if(verbose) {cat(paste("Importing",length(Samples),"allele counter files"),sep="\n")}
  NV_and_NR=lapply(output_files,function(file) {
    if(verbose) {cat(file,sep="\n")}
    dat=read.delim(file,stringsAsFactors = F)%>%
      dplyr::rename("Chrom"=X.CHR,"Pos"=POS)%>%
      left_join(bed_file,by=c("Chrom","Pos"))
    sample_NR=dat$Good_depth
    sample_NV=sapply(1:nrow(dat),function(i) {return(dat[i,paste0("Count_",dat$Alt[i])])})
    return(list(NV=sample_NV,NR=sample_NR))
  })
  
  #Get the standard set of mut_refs (in the same order as allele counter output) by looking at the first allele counter file
  if(verbose) {cat("Getting the 'mut_refs' in the correct order",sep="\n")}
  mut_refs=read.delim(output_files[1],stringsAsFactors = F)%>%
    dplyr::rename("Chrom"=X.CHR,"Pos"=POS)%>%
    left_join(bed_file,by=c("Chrom","Pos"))%>%
    tidyr::unite(col="mut_ref",Chrom,Pos,Ref,Alt,sep="-")%>%
    pull(mut_ref)
  
  #Bind the NV output and NR outputs into separate matrices
  if(verbose) {cat("Binding the NV and NR files into separate matrices",sep="\n")}
  NV=as.matrix(dplyr::bind_cols(lapply(NV_and_NR,function(list) list$NV)))
  NR=as.matrix(dplyr::bind_cols(lapply(NV_and_NR,function(list) list$NR)))
  
  #Name the columns as the samples and rows as the mutation references ('mut_refs')
  rownames(NV)=rownames(NR)<-mut_refs
  colnames(NV)=colnames(NR)<-Samples
  
  return(list(NV=NV,NR=NR))
}

calculate_cell_frac=function(NV,NR,sex="male") {
  #Remove 0's from the depth, to avoid dividing by 0
  NR[NR==0]<-1
  #For autosomal mutations, cell frac = NV/ (NR/2)
  cell_frac=NV/(NR/2)
  
  if(sex=="Male"|sex=="male"|sex=="M") {
    #Get vectors to select out the autosomal and XY chromosomal mutations
    XY_muts=grepl("X",rownames(NR))|grepl("Y",rownames(NR))
    #For XY mutations, cell frac = NV/NR -> replace these accordingly
    cell_frac[XY_muts,]<-(NV[XY_muts,])/(NR[XY_muts,])
  }
  #Cell frac cannot be greater than 1, therefore if comes out as > 1 (which can happen when dividing the NR by 2), coerce to 1
  cell_frac[cell_frac>1]<-1
  return(cell_frac)
}

get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

aggregate_cols=function(NV,NR,metadata,ID=NULL,Tissue_IDs=list(GranMono=c("CD11b","GM"),Bcells=c("B220"),Tcells=c("CD3e"))) {
  Samples=colnames(NV)
  Samples<-Samples[!Samples=="MDGRCm38is"]
  conv_samples=sapply(Samples,function(samp) {metadata%>%filter(INTERNAL_CASM_SAMPLE_NAME==samp)%>%pull(SUPPLIER_SAMPLE_ID)})
  
  sample_lists<-lapply(Tissue_IDs,function(labels){
    temp_list<-lapply(labels,function(label) grepl(label,conv_samples))
    select_vec<-Reduce(f=function(a,b) {a|b},temp_list)
    return(conv_samples[select_vec])
  })
  
  NV_agg=Map(labels=Tissue_IDs,matching_samples=sample_lists,function(labels,matching_samples){
    Tissue_NV<-rowSums(NV[,names(matching_samples),drop=F])
    return(Tissue_NV)
  })%>%dplyr::bind_cols()%>%
    as.matrix()
  NR_agg=Map(labels=Tissue_IDs,matching_samples=sample_lists,function(labels,matching_samples){
    Tissue_NR<-rowSums(NR[,names(matching_samples),drop=F])
    return(Tissue_NR)
  })%>%dplyr::bind_cols()%>%
    as.matrix()
  rownames(NV_agg)=rownames(NR_agg)=rownames(NV)
  
  if(!is.null(ID)) {
    colnames(NV_agg)=colnames(NR_agg)=paste(ID,names(Tissue_IDs),sep="_")
  } else {
    olnames(NV_agg)=colnames(NR_agg)=names(Tissue_IDs)
  }
  
  return(list(NV=NV_agg,NR=NR_agg,samples_included=sample_lists))
}

read_posterior_VAFs_return_cell_frac=function(posterior_VAFs_file) {
  require(readr)
  require(dplyr)
  tissue_post<-readr::read_delim(posterior_VAFs_file,delim = "\t",show_col_types = FALSE)
  
  VAF_distributions=tissue_post%>%
    dplyr::select(Posterior_VAFs)%>%
    separate(col="Posterior_VAFs",into=paste("n",1:100,sep="_"),sep=",")%>%
    mutate_all(as.numeric)
  
  cell_frac_distributions=2*VAF_distributions
  cell_frac_post<-bind_cols((tissue_post%>%dplyr::select(Node_assignment,mutation_ID)),cell_frac_distributions)
  return(cell_frac_post)
}

get_95percentCI_cellfrac<-function(post.df,CI_level=0.95) {
  require(dplyr)
  min_quantile=(1-CI_level)/2
  max_quantile=1-((1-CI_level)/2)
  quantile_cell_fracs=post.df%>%
    dplyr::select(-Node_assignment,-mutation_ID)%>%
    as.matrix()%>%
    apply(1,function(x) {quantile(x,c(min_quantile,0.5,max_quantile))})%>%t()
  return(cbind(post.df%>%dplyr::select(mutation_ID),as.data.frame(quantile_cell_fracs)))
}

get_median_cellfracs=function(post.df) {
  require(dplyr)
  median_cell_fracs=post.df%>%
    dplyr::select(-Node_assignment,-mutation_ID)%>%
    as.matrix()%>%
    apply(1,median)
  return(data.frame(mutation_ID=post.df$mutation_ID,median_cell_frac=median_cell_fracs))
}

#Get VAFs of all clones at 100 mutations of molecular time
get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

get_sum_of_frac=function(tissue_post,tree,cut_off,verbose=F) {
  
  if(verbose) {cat(paste("Clone cutoff of",cut_off,"being used"),sep="\n")}
  
  #Define the cutoff branches
  cutoff_branches=get_cutoff_branches(tree,cut_off = cut_off)
  
  if(verbose) {cat(paste(length(cutoff_branches)," branches at the cutoff level"),sep="\n")}
  
  #If the cutoff branches have no mutations, need to replace with the daughter branches of that node
  number_of_covered_muts=function(node,tissue_post) {sum(tissue_post$Node_assignment==node)} #convenience function
  
  while(any(sapply(cutoff_branches,number_of_covered_muts,tissue_post=tissue_post)==0)){
    updated_nodes<-unlist(lapply(cutoff_branches,function(node) {
      if(number_of_covered_muts(node,tissue_post=tissue_post)==0) {
        if(verbose) {cat(paste("No mutations on node",node,"=> replacing with daughter nodes."),sep="\n")}
        daughter_nodes<-tree$edge[,2][tree$edge[,1]==node]
        return(daughter_nodes)
      } else {
        return(node)
      }
    }))
    cutoff_branches<-updated_nodes
  }
  
  #Define the clone posteriors (as a list)
  clone_posteriors<-lapply(cutoff_branches,function(node) {
    if(verbose) {cat(node,sep="\n")}
    branch_heights=nodeHeights(tree)[tree$edge[,2]==node,]
    min_height=branch_heights[1]
    max_height=branch_heights[2]
    
    if(min_height>cut_off) {
      prop<-0
    } else {
      prop=(cut_off-min_height)/(max_height-min_height)
    }
    
    n_targseq_muts=sum(tissue_post$Node_assignment==node)
    
    which_branch_rank=max(round(prop*n_targseq_muts),1)
    node_cell_frac_distributions=tissue_post%>%
      filter(Node_assignment==node)%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()
    
    #Sort each iteration of posterior by cell fraction & choose that which corresponds to the rank
    if(nrow(node_cell_frac_distributions)>1) {
      #Original method - sorts cell fracs of each Gibbs sample, and takes the 10th rank of each. ?significant overestimates
      #node_cell_frac_distributions_sorted<-apply(node_cell_frac_distributions,2,sort,decreasing=T)
      
      #Second method - keeps values of specific mutations linked together & just orders them by their median values
      rank=order(apply(node_cell_frac_distributions,1,median),decreasing = T)
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions[rank,]
      
    } else {
      node_cell_frac_distributions_sorted<-node_cell_frac_distributions
    }
    
    clone_post=node_cell_frac_distributions_sorted[which_branch_rank,]
    return(clone_post)
  })
  names(clone_posteriors)<-paste("node",cutoff_branches,sep="_")
  
  #Use these to get the sum of VAF
  n_iter=length(clone_posteriors[[1]])
  sum_of_frac_post<-sapply(1:100,function(j) {
    sum(sapply(clone_posteriors,function(x) x[j]))
  })
  return(sum_of_frac_post)
}


Gibbs_targ_seq_plots=function(SampleID,
                              tree,
                              details_targ,
                              pair_cell_fracs,
                              scale_muts_to_branch=TRUE,
                              colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                              log_min=-5,
                              vaf_lwd=5,
                              title=NULL) {
  
  details_targ_full<-left_join(details_targ,
                               get_median_cellfracs(post.df=pair_cell_fracs),
                               by=c("mut_ref"="mutation_ID"))
  
  #Generate the rescaled log cell fraction for plotting with contrast
  details_targ_full$log_median_cell_frac<-log(details_targ_full$median_cell_frac)
  details_targ_full$log_median_cell_frac[details_targ_full$log_median_cell_frac<log_min]<-log_min
  details_targ_full$log_median_cell_frac=plotrix::rescale(details_targ_full$log_median_cell_frac,newrange = c(0,1))
  
  ##Generate the plot
  tree=plot_tree(tree,cex.label=F,title=title)
  add_annotation(tree=tree,
                 details=details_targ_full,
                 matrices=NULL,
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,
                               details,
                               node,
                               var_field = "log_median_cell_frac",
                               lwd = vaf_lwd,
                               colours=colour.scale,
                               scale_muts_to_branch=scale_muts_to_branch)
                 }
  )
}


Gibbs_targ_seq_comparison_plots=function(post1,
                                         post2,
                                         details_targ,
                                         title=NULL,
                                         vaf_lwd=5,
                                         log_min=-4) {
  
  require(dplyr)
  require(RColorBrewer)
  
  #Ensure post1 and post2 are in the same order (in terms of mutations)
  post2<-post2[match(post1$mutation_ID,post2$mutation_ID),]
  
  #Remove the mutation_ID/ node assignnment columns so that you only have the posterior cell fractions
  post1_mat<-dplyr::select(post1,-mutation_ID,-Node_assignment)%>%as.matrix()
  post2_mat<-dplyr::select(post2,-mutation_ID,-Node_assignment)%>%as.matrix()
  post_comp_mat<-post1_mat/post2_mat
  
  post_comp<-dplyr::bind_cols(post1%>%dplyr::select(Node_assignment,mutation_ID),as.data.frame(post_comp_mat))
  
  CI_df<-get_95percentCI_cellfrac(post_comp)
  
  #Get the median fold change
  details_targ_full<-left_join(details_targ,
                               get_median_cellfracs(post.df=post_comp),
                               by=c("mut_ref"="mutation_ID"))
  
  #Take the log of this, such that 0 = no difference
  details_targ_full$log_median_cell_frac<-log(details_targ_full$median_cell_frac)
  
  #Take the absolute values, to allow the same scaling for positive and negative values but keep 0s the same
  details_targ_full$abs_log_median_cell_frac<-abs(details_targ_full$log_median_cell_frac)
  
  #Rescale such that the maximum absolute log change = 1
  details_targ_full$abs_log_median_cell_frac_scaled<-plotrix::rescale(details_targ_full$abs_log_median_cell_frac,newrange = c(0,1))
  
  #Now convert back to whether the change is positive/ negative be taking the sign of original log value
  details_targ_full$log_median_cell_frac_scaled<-details_targ_full$abs_log_median_cell_frac_scaled*sign(details_targ_full$log_median_cell_frac)
  
  #Now centre the values on 0.5 (rather than 0) as this fits the colour scale used and the way that the add_var_col function works
  #(the value for no change should be 0.5, and values should be between minimum of 0 and maximum of 1)
  details_targ_full$final<-sapply(details_targ_full$log_median_cell_frac_scaled,function(x) {return((x+1)/2)})
  
  #Censor results if not significant (i.e. the 95% confidence intervals span 1 - no difference)
  details_targ_full$final[CI_df$`2.5%`<1&CI_df$`97.5%`>1]<-0.5
  
  tree=plot_tree(tree,cex.label=F,title=title)
  temp=add_annotation(tree=tree,
                      details=details_targ_full,
                      matrices=NULL,
                      annot_function=function(tree,details,matrices,node) {
                        add_var_col(tree,
                                    details,
                                    node,
                                    var_field = "final",
                                    lwd = 5,
                                    colours=RColorBrewer::brewer.pal(11,"PRGn"),
                                    scale_muts_to_branch=T)
                      }
  )
}


add_var_col=function(tree, ##<< enhanced phylo returned from plot_tree
                     details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                     node,
                     var_field,
                     b.add.line=TRUE,
                     colours = c("black","green","red"),
                     scale_muts_to_branch=TRUE,
                     ...){
  
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  muts_on_edge=length(info$idx.in.details)
  edge_length=tree$edge.length[tree$edge[,2]==node]
  
  if(muts_on_edge > 0 & edge_length>0) {
    bdat=details[[var_field]][info$idx]
    if(is.null(bdat) || class(bdat)!="numeric"){
      stop("Error in provided bfield (does it exist and is it numeric?)")
    }
    
    bdat = sort(bdat, decreasing = TRUE)
    if(scale_muts_to_branch) {
      mut_unit_of_edge=edge_length/muts_on_edge
    } else {
      mut_unit_of_edge=1
    }
    ##Could add in a third category NA
    #missing=sum(is.na(bdat))
    if(b.add.line){
      y0_next = info$yt
      for(i in 1:muts_on_edge) {
        arrows(y0=y0_next,y1=(y0_next - mut_unit_of_edge),x0=info$x,length = 0,col=col.scale[ceiling(100*bdat[i])],lend=1,...)
        y0_next = y0_next - mut_unit_of_edge
      }
    }
  }
}


root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/mouse_phylo/targeted/","/lustre/scratch126/casm/team154pc/ms56/mouse_phylo/targeted/")

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

metadata<-readxl::read_excel(paste0("sample_metadata/targeted_metadata.xlsx"))%>%
  mutate(flow_marker=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,1],
         age=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,3],
         PDID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))

baitset_SNV_bed_file_path=paste0(root_dir,"baitset_SNVs.bed")
baitset_INDEL_bed_file_path=paste0(root_dir,"baitset_INDELs.bed")
allele_counter_dir=paste0(root_dir,"allele_counter/")
cgpvaf_SNV_file=paste0(root_dir,"merged_SNVs_targeted.tsv")
cgpvaf_INDEL_file=paste0(root_dir,"merged_indels_targeted.tsv")
plots_dir=paste0(root_dir,"plots/")
matrices_files=paste0(root_dir,"full_matrices.RDS")
matrices_files_cgpVAF=paste0(root_dir,"full_cgpvaf_matrices.RDS")
output_dir=paste0(root_dir,"output/")

if(Sys.info()['sysname'] == "Darwin"){
  tree_folder=paste0(my_working_directory,"/tree_files/")
  annotated_muts_folder=paste0(my_working_directory,"/annotated_muts/")
} else {
  tree_folder=paste0(my_working_directory,"filtering_runs/tree_files")
  annotated_muts_folder=paste0(my_working_directory,"filtering_runs/annotated_muts")
}

#There are two sets of MDIDs as the targeted sequencing samples were mistakenly given new IDs
Phylo_MDIDs=c("MD7634","MD7635")
MDIDs=c("MD7816","MD7817")
names(MDIDs)<-Phylo_MDIDs

all_targeted_res=vector(mode = "list",length = length(MDIDs))
names(all_targeted_res)<-MDIDs

#Set up the all.muts and all.trees objects - if previously run, read in the RDS object, otherwise import the details matrices from the annotated muts objects
cat("Loading the individual annotated mutation files",sep = "\n")
annotated_muts_paths=list.files(annotated_muts_folder,pattern="_vaf_post_mix_post_dup",full.names = T)
all.muts<-lapply(Phylo_MDIDs,function(ID){
  file_paths<-grep(ID,annotated_muts_paths,value = T)
  load(grep("_m40_",file_paths,invert=T,value = T));return(filtered_muts$COMB_mats.tree.build$mat)
})

tree_paths=list.files(tree_folder,pattern="_vaf_post_mix_post_dup.tree",full.names = T)
all.trees<-lapply(Phylo_MDIDs,function(ID){
  tree_paths<-grep(ID,tree_paths,value = T)
  tree<-read.tree(grep("_m40_",tree_paths,invert=T,value = T))
  return(tree)
})
names(all.muts)<-names(all.trees)<-Phylo_MDIDs

#Import the allele counter targeted sequencing data
if(!file.exists(matrices_files)){
  #Import the SNV bedfile and name the columns appropriately
  SNV_mut_refs=read.delim(baitset_SNV_bed_file_path,stringsAsFactors = F,header=F)
  colnames(SNV_mut_refs)<-c("Chrom","Pos","Ref","Alt")
  
  matrices=import_allele_counter_data(allele_counter_dir = allele_counter_dir, bed_file=SNV_mut_refs)
  saveRDS(matrices,file=matrices_files)
} else {
  matrices<-readRDS(matrices_files)
}

# Import the cgpVAF targeted sequencing data ----

if(T|!file.exists(matrices_files_cgpVAF)){
  cat("Importing cgpVAF files",sep="\n")
  matrices_cgpvaf=import_cgpvaf_SNV_and_INDEL(SNV_output_file = cgpvaf_SNV_file,INDEL_output_file = cgpvaf_INDEL_file)
  saveRDS(matrices_cgpvaf,file=matrices_files_cgpVAF)
} else {
  cat("Reading in previously saved cgpVAF files",sep="\n")
  matrices_cgpvaf<-readRDS(matrices_files_cgpVAF)
}

for(ID in Phylo_MDIDs){
  cat(ID,sep="\n")
  details=all.muts[[ID]]
  tree=all.trees[[ID]]
  
  #In details file, create column for whether it is included in the targeted sequencing baits
  details$targeted_tree = ifelse(details$mut_ref %in% rownames(matrices_cgpvaf$NV),"YES","NO")
  
  #Get the details matrix only including mutations in the targeted panel
  details_targ = details[details$targeted_tree=="YES",]
  
  #Get the tissue names from this mouse
  this_ID_tissues=metadata%>%filter(PDID==MDIDs[ID])%>%pull(INTERNAL_CASM_SAMPLE_NAME)%>%unique()
  
  #Remove tissue names that failed targeted sequencing
  this_ID_tissues<-this_ID_tissues[this_ID_tissues%in%colnames(matrices_cgpvaf[['NV']])]
  
  matrices_ID=lapply(matrices_cgpvaf[c("NV","NR")],function(mat) {
    mat[details_targ$mut_ref,this_ID_tissues]
  })
  
  all_targeted_res[[ID]]<-list(details_targ=details_targ,matrices=matrices_ID,smry=(metadata%>%filter(PDID==MDIDs[ID])))
}

#--------------------------------------------------------------------------------------------------#
## Plot the coverage in the two mice ----
#--------------------------------------------------------------------------------------------------#

mean_cov_df<-lapply(c("MD7634","MD7635"),function(MDID) {
  res<-all_targeted_res[[MDID]]
  mean_cov<-apply(res$matrices$NR,2,mean)
  median_cov<-apply(res$matrices$NR,2,median)
  data.frame(tissueID=names(mean_cov),mean_cov=mean_cov,median_cov=median_cov)
})%>%dplyr::bind_rows()

left_join(metadata,mean_cov_df,by=c("INTERNAL_CASM_SAMPLE_NAME"="tissueID"))%>%
  mutate(pheno_simple=ifelse(TISSUE_PHENOTYPE=="Blood","Blood","Other"))%>%
  ggplot(aes(x=median_cov,fill=pheno_simple))+
  geom_histogram(col="black",linewidth=0.5,bins = 20)+
  #scale_x_log10(breaks=c(20,30,50,100,200,300,500))+
  facet_grid(rows=vars(PDID))+
  theme_classic()+
  labs(x="Median coverage (X)",fill="")

#--------------------------------------------------------------------------------------------------#
## -----Get output from the aggregated info in format for the Gibbs sampler to run ------
#--------------------------------------------------------------------------------------------------#

ID=Phylo_MDIDs[2]

#Use one of the samples to get the tree & details matrices
temp_samp<-metadata%>%filter(PDID==MDIDs[ID])%>%pull(INTERNAL_CASM_SAMPLE_NAME)%>%.[3]
dat<-readRDS(paste0(root_dir,"/GibbsSampler/data/",temp_samp,"_gibbs_info.RDS"))
details_targ=dat$details%>%tidyr::unite(col=mut_ref,Chrom,Pos,Ref,Alt,sep="-",remove=F)
tree=squash_tree(dat$tree)

#Subset the counts matrices to include only the mutations and samples from that mouse - the "test" matrix
matrices_test=lapply(matrices_cgpvaf,function(mat) {
  mat[details_targ$mut_ref,grepl(MDIDs[ID],colnames(mat))]
})

#Subset the counts matrices to include only the mutations from that mouse, but the samples from the OTHER mouse - the "control" matrix (these mutations are not expected to be present in the other mouse)
matrices_ctrl=lapply(matrices_cgpvaf,function(mat){
  mat[details_targ$mut_ref,!grepl(MDIDs[ID],colnames(mat))]
})

#This is a list of all the different ways the data will be aggregated
aggregation_info<-list(across_time=list(GranMono=c("CD11b\\+_\\d","GM"),Bcells=c("B220"),Tcells=c("CD3e"),BandT=c("B220","CD3e")),
                       by_time_point=list(nine_week_blood=c("_9w"),fifteen_week_blood=c("_15w"),twentyone_week_blood="_21w",twentyseven_week_blood="_27w"),
                       by_tissue=list(LymphNode=c("_LN_"),Spleen=c("_Sp_"),Lung=c("_Lu_"),Thymus="_Thy_",Skin="Sk_",Peritoneal="_PC_",Liver="_Li_"),
                       by_cell_type=list(iNKT=c("iNKT_"),ILC2="ILC2_",ILC3="ILC3_",allILC="ILC"))


SNVs_only=T
resave=F
all_aggregated_matrices=lapply(aggregation_info,function(Tissue_IDs) {
  matrices_aggregated<-aggregate_cols(NV=matrices_test$NV,NR=matrices_test$NR,metadata=metadata,ID=ID,Tissue_IDs = Tissue_IDs)
  
  ##Now restructure the data in the format for Nick's version of the Gibbs sampler and save for each individual sample
  details_targ<-all_targeted_res[[ID]]$details_targ
  
  if(SNVs_only) {
    details_targ<-details_targ%>%filter(nchar(Ref)==1 & nchar(Alt)==1)
  }
  
  this_mouse_IDs<-colnames(matrices_aggregated$NV)
  
  #Get the control counts (these will be the same for all samples from this mouse)
  other_mouse_IDs<-grep(MDIDs[names(MDIDs)!=ID],colnames(matrices_cgpvaf$NV),value=T)
  mtr_other<-rowSums(matrices_cgpvaf[["NV"]][details_targ$mut_ref,other_mouse_IDs])
  depth_other<-rowSums(matrices_cgpvaf[["NR"]][details_targ$mut_ref,other_mouse_IDs])
  
  
  for(tissueID in this_mouse_IDs) {
    cat(tissueID,sep="\n")
    details_tissue<-cbind(details_targ[,c("Chrom","Pos","Ref","Alt","node")],data.frame(mtr=matrices_aggregated[["NV"]][details_targ$mut_ref,tissueID],
                                                                                        depth=matrices_aggregated[["NR"]][details_targ$mut_ref,tissueID],
                                                                                        mtr_other=mtr_other,
                                                                                        depth_other=depth_other))
    
    GS_data_file<-paste0(root_dir,"GibbsSampler/data/",tissueID,"_gibbs_info.RDS")
    if(!file.exists(GS_data_file)|resave) {
      saveRDS(object=list(details=details_tissue,tree=all.trees[[ID]]),file = GS_data_file)
    }
  }
  
  return(matrices_aggregated)
})

##Run the Gibbs sampler from R for these samples
setwd(paste0(root_dir,"GibbsSampler/"))

#Loop through the aggregated matrices
for(i in 1:length(all_aggregated_matrices)) {
  mat<-all_aggregated_matrices[[i]]
  samples<-colnames(mat$NV)
  for (ID in samples) {
    output_file=paste0("output/",ID,"_branch_VAFs.txt.gz")
    if(!file.exists(output_file)) {
      cat(paste("Output files for",ID,"not found. Running the Gibbs Sampler algorithm."),sep="\n")
      system(paste("julia wrap_gibbs.jl",ID))
    } else {
      cat(paste("Output files for",ID,"already exist. Not re-running"),sep="\n")
    }
  }
}

## AFTER RUNNING THE GIBBS SAMPLER IN JULIA ----

## This code imports the posteriors for two samples/ cells types to do comparisons
ID="MD7635"
comparison=c("GranMono","MD7817bt")

tissue_names<-sapply(comparison,function(x) {
  if(grepl("^MD",x)) {
    deets<-metadata%>%filter(INTERNAL_CASM_SAMPLE_NAME==x)%>%dplyr::select(SUPPLIER_SAMPLE_ID,TISSUE_PHENOTYPE)
    return(paste(unlist(deets),collapse = " "))
  } else {
    return(x)
  }
})

post<-lapply(comparison,function(tissue_type) {
  
  if(grepl("^MD",tissue_type)) {tissueID<-tissue_type} else {tissueID=paste(ID,tissue_type,sep="_")}
  
  posterior_VAFs_file=paste0(root_dir,"/GibbsSampler/output/",tissueID,"_posterior_VAFs.txt.gz")
  cell_frac_post<-read_posterior_VAFs_return_cell_frac(posterior_VAFs_file)
  return(cell_frac_post)
})

#Import the original data to get the tree & details objects
tissueID=paste(ID,comparison[1],sep="_")

dat<-readRDS(paste0(root_dir,"/GibbsSampler/data/",tissueID,"_gibbs_info.RDS"))
tree=dat$tree
tree=squash_tree(tree,cut_off = 50)
details=dat$details%>%tidyr::unite(col=mut_ref,Chrom,Pos,Ref,Alt,sep="-",remove=F)

#Plot the two posteriors separately
par(mfrow=c(1,2))
Gibbs_targ_seq_plots(SampleID=comparison[1],
                     tree=tree,
                     details_targ=details,
                     pair_cell_fracs = post[[1]],
                     log_min = -3,
                     title = tissue_names[1],
                     vaf_lwd=5)

Gibbs_targ_seq_plots(SampleID=comparison[2],
                     tree=tree,
                     details_targ=details,
                     pair_cell_fracs = post[[2]],
                     log_min = -3,
                     title = tissue_names[2],
                     vaf_lwd=5)

#Now plot the comparison (green = post1 is higher than post2; purple = post2 is higher than post1)
par(mfrow=c(1,1))
Gibbs_targ_seq_comparison_plots(post1=post[[1]],
                                post2=post[[2]],
                                details_targ = details,title = paste(tissue_names[1],"/",tissue_names[2],collapse=" "))

## This code imports all posteriors from both mice
posterior_cell_fracs=vector(mode="list",length = length(MDIDs))
names(posterior_cell_fracs)<-MDIDs
for(i in 1:length(MDIDs)) {
  posterior_VAFs_files=list.files(paste0(root_dir,"/GibbsSampler/output/"),pattern="_posterior_VAFs.txt.gz")
  this_mouse_posteriors=grep(paste0(MDIDs[i],"|",names(MDIDs)[i],collapse=""),posterior_VAFs_files,value=T)
  
  this_mouse_post_cell_fracs<-lapply(this_mouse_posteriors,function(file_name) {
    cat(file_name,sep="\n")
    posterior_VAFs_file=paste0(root_dir,"/GibbsSampler/output/",file_name)
    post_cell_fracs<-read_posterior_VAFs_return_cell_frac(posterior_VAFs_file)
    return(post_cell_fracs)
  })
  names(this_mouse_post_cell_fracs)<-gsub("_posterior_VAFs.txt.gz","",this_mouse_posteriors)
  posterior_cell_fracs[[i]]<-this_mouse_post_cell_fracs
}

#--------------------------------------------------------------------------------------------------#
## -----Set theme, import metadata ------
#--------------------------------------------------------------------------------------------------#

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

metadata<-readxl::read_excel(paste0(my_working_directory,"sample_metadata/targeted_metadata.xlsx"))%>%
  mutate(flow_marker=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,1],
         age=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,3],
         PDID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))

ID="MD7817"
pass_targeted_ids<-grep("test",invert=T,value=T,gsub("_gibbs_info.RDS","",list.files(paste0(root_dir,"/GibbsSampler/data"))))
completed_ids<-grep("test",invert=T,value=T,gsub("_posterior_VAFs.txt.gz","",list.files(paste0(root_dir,"/GibbsSampler/output"),pattern = "posterior_VAFs.txt")))
pass_targeted_ids[!pass_targeted_ids%in%completed_ids]

input_dat<-lapply(this_mouse_tissueIDs,function(tissueID) {
  input_dat_file=paste0(root_dir,"/targeted/GibbsSampler/data/",tissueID,"_gibbs_info.RDS")
  dat<-readRDS(input_dat_file)
  return(dat)
})

mean(input_dat[[1]]$details$depth)
hist(sapply(input_dat,function(dat) mean(dat$details$depth)))

this_mouse_tissueIDs<-metadata%>%filter(PDID==ID & INTERNAL_CASM_SAMPLE_NAME%in%completed_ids)%>%pull(INTERNAL_CASM_SAMPLE_NAME)
mouse_posteriors<-lapply(this_mouse_tissueIDs,function(tissueID) {
  posterior_VAFs_file=paste0(root_dir,"/GibbsSampler/output/",tissueID,"_posterior_VAFs.txt.gz")
  tissue_posteriors<-readr::read_delim(posterior_VAFs_file,delim = "\t",show_col_types = FALSE)
  return(tissue_posteriors)
})

post_fracs<-lapply(mouse_posteriors,function(tissue_post) {
  VAF_distributions=tissue_post%>%
    dplyr::select(Posterior_VAFs)%>%
    separate(col="Posterior_VAFs",into=paste("n",1:100,sep="_"),sep=",")%>%
    mutate_all(as.numeric)
  
  cell_frac_distributions=2*VAF_distributions
  
  cell_frac_post<-bind_cols((tissue_post%>%dplyr::select(Node_assignment,mutation_ID)),cell_frac_distributions)
  return(cell_frac_post)
})

i=14
tissueID<-this_mouse_tissueIDs[i]
dat<-readRDS(paste0(root_dir,"/GibbsSampler/data/",tissueID,"_gibbs_info.RDS"))
tree=dat$tree
tree=squash_tree(tree,cut_off = 50)
details=dat$details%>%tidyr::unite(col=mut_ref,Chrom,Pos,Ref,Alt,sep="-",remove=F)

Gibbs_targ_seq_plots(SampleID=tissueID,
                     tree=tree,
                     details_targ=details,
                     pair_cell_fracs = post_fracs[[i]],
                     log_min = -4,
                     title = tissueID,
                     vaf_lwd=5)

#--------------------------------------------------------------------------------------------------#
# Sum of cell fraction ----
#--------------------------------------------------------------------------------------------------#

both_trees<-lapply(MDIDs,function(ID) {
  this_mouse_files<-grep(ID,list.files(paste0(root_dir,"/GibbsSampler/data/")),value=T)
  readRDS(paste0(root_dir,"/GibbsSampler/data/",this_mouse_files[1]))$tree
})

#Calculate how much of the tissues cells are accounted for by the captured lineages
sum_of_frac_list_file=paste0(root_dir,"sum_of_frac.RDS")
cutoffs_to_use<-seq(2,50,2)
if(file.exists(sum_of_frac_list_file)){
  sum_of_frac_list<-readRDS(sum_of_frac_list_file)
} else {
  #This is slow
  sum_of_frac_list<-Map(tree=both_trees,post=posterior_cell_fracs,ID=MDIDs,function(tree,post,ID) {
    cat(ID,sep="\n")
    tree.ultra<-make.ultrametric.tree(tree)
    tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(drop.tip(tree,"Ancestral")))
    tree.ultra$edge.length[tree.ultra$edge[,2]==which(tree.ultra$tip.label=="Ancestral")]<-0
    
    tissueIDs=names(post)
    this_mouse_sum_of_frac_by_cutoff<-lapply(tissueIDs,function(tissueID){
      cat(tissueID,sep="\n")
      tissue_post<-post[[tissueID]]
      tissue_sum_of_frac_by_cutoff<-lapply(cutoffs_to_use,function(clone_cutoff) {
        cat(clone_cutoff,sep="\n")
        get_sum_of_frac(tissue_post = tissue_post,cut_off = clone_cutoff,tree = tree.ultra,verbose=F)
      })
      return(tissue_sum_of_frac_by_cutoff)
    })
    names(this_mouse_sum_of_frac_by_cutoff)<-tissueIDs
    return(this_mouse_sum_of_frac_by_cutoff)
  })
  saveRDS(sum_of_frac_list,file=sum_of_frac_list_file)
}

##Extract the median and 95% CI of the sum of frac statistic
all_sof<-Map(ID=MDIDs,this_mouse_sum_of_frac_by_cutoff=sum_of_frac_list,f=function(ID,this_mouse_sum_of_frac_by_cutoff) {
  mouse_cell_fracs<-Map(tissue=this_mouse_sum_of_frac_by_cutoff,tissueID=names(this_mouse_sum_of_frac_by_cutoff),f=function(tissue,tissueID) {
    cat(tissueID,sep="\n")
    tissue_df<-Map(cutoff=cutoffs_to_use,cutoff_sof=tissue,function(cutoff,cutoff_sof) {
      data.frame(cutoff=cutoff,lowerCI=quantile(cutoff_sof,0.025),median=median(cutoff_sof),upperCI=quantile(cutoff_sof,0.975))
    })%>%dplyr::bind_rows()%>%
      mutate(tissueID=tissueID,.before=1)
    return(tissue_df)
  })%>%dplyr::bind_rows()%>%
    mutate(MDID=ID,.before=1)
  return(mouse_cell_fracs)
})

##Plot the results
SOF_plot<-all_sof%>%
  dplyr::bind_rows()%>%
  dplyr::filter(grepl("_GranMono|_BandT",tissueID))%>%
  dplyr::mutate(tissueID=stringr::str_split(tissueID,"_",simplify=T)[,2])%>%
  ggplot(aes(ymin=lowerCI,ymax=upperCI,y=median,x=cutoff,fill=tissueID,col=tissueID))+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(alpha=0.2,col=NA)+
  geom_line(size=0.3)+
  facet_grid(~MDID)+
  scale_y_continuous(breaks=seq(0,1,0.5),limits=c(0,1))+
  theme_bw()+
  my_theme+
  labs(x="Molecular time",y="Captured cell fraction")

ggsave(filename=paste0(plots_dir,"Sum_of_frac_through_time.pdf"),SOF_plot,width=6,height=8)


#========================================#
# ASSESS SIMILARITY OF CLONAL COMPOSITION BETWEEN TISSUES ####
#========================================#
##One way of assessing similarity of clonal composition between tissues is using the
##soft cosine similarity of the targeted sequencing data between tissues/ individuals

##Define custom functions for this analysis
generate_median_fracs_across_tissues=function(post_ind,selected_nodes) {
  tissues=names(post_ind)
  mutation_order=post_ind[[1]]%>%filter(Node_assignment%in%selected_nodes)%>%pull(mutation_ID)
  temp<-lapply(post_ind,function(post) {
    post_emb<-post%>%filter(Node_assignment%in%selected_nodes)
    cell_frac_distributions=post_emb%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()
    median_cell_fracs<-apply(cell_frac_distributions,1,median)
    names(median_cell_fracs)<-post_emb$mutation_ID
    return(median_cell_fracs)
  })
  all_tissue_median_fracs<-Reduce(cbind,lapply(temp,function(x) x[mutation_order]))
  dimnames(all_tissue_median_fracs)[[2]]<-tissues
  return(all_tissue_median_fracs)
}

##Function to characterise the shortest phylogenetic distance between mutations
#Uses a reference object for distances between nodes
shortest_dist = function(mut_a,mut_b, tree, ref,node_comb_dist) { #the function to look up the distance
  node_a=ref$Node_assignment[ref$mutation_ID==mut_a]
  node_b=ref$Node_assignment[ref$mutation_ID==mut_b]
  genetic_dist<-node_comb_dist[[paste(node_a,node_b,sep="-")]]
  return(genetic_dist)
}

##Define functions to calculate the cosine similarity for tissue pairs
sum_tissue_scores = function(data,selected_muts,tissue_1,tissue_2,mut_pairs_sim) {
  tissue_pairs_df=as.data.frame(data.table::CJ(a=data[tissue_1,selected_muts],b=data[tissue_2,selected_muts],unique=FALSE,sorted = FALSE))
  tissue_pairs_df$sim=mut_pairs_sim
  tissue_pairs_df$numerator=apply(tissue_pairs_df[,1:3],1,prod)
  return(sum(tissue_pairs_df$numerator)/2)
}

get_soft_cosim=function(data,selected_muts,tissue_1,tissue_2,mut_pairs_sim) {
  num=sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_1,tissue_2,mut_pairs_sim = mut_pairs_sim)
  denom=sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_1,tissue_1,mut_pairs_sim = mut_pairs_sim)^0.5 * sum_tissue_scores(data=data,selected_muts=selected_muts,tissue_2,tissue_2,mut_pairs_sim = mut_pairs_sim)^0.5
  soft_cosim=num/denom
  return(soft_cosim)
}

## Generate the similarity matrices -----
## For this analysis we focus on the early embryonic mutations (<30 mutations of molecular time, corresponds to ~1st trimester)
# Most of the cell compartments will have most information on these early mutations
# The analysis also becomes very costly if start to have too many mutations (>1,000) - therefore if there are to many, need to reduce
# Could do this by random subsampling, but might miss key early branches, therefore reduce by gradually reducing the cut-off time until the number of mutations is <1,000

#This step takes some time so we have pre-saved the output to this
similarity_matrices_file=paste0(root_dir,"/data/similarity_matrices.Rds")
if(!file.exists(similarity_matrices_file)) {
  similarity_matrices=list()
  for(j in 1:length(posterior_cell_fracs)) {
    embryonic_mutation_cutoff=30
    tree=all.trees.cc.nodups[[j]]
    pair=names(all.trees.cc.nodups)[j]
    post_ind=posterior_cell_fracs[[j]]
    cat(pair,sep="\n")
    
    #Get vector of the tissue IDs that relate to the initial time point only
    initial_time_point_samples<-bulk_smry_all%>%filter(Pair==pair & time_point==0)%>%pull(tissueID)%>%unique()
    if(!all(names(post_ind)%in%initial_time_point_samples)) {
      post_ind<-post_ind[-which(!names(post_ind)%in%initial_time_point_samples)]
    }
    
    ##Select which branches are definitively 'embryonic' - i.e. all mutations are from early life
    embryonic_branches=tree$edge[,2][which(nodeHeights(tree)[,2]<embryonic_mutation_cutoff)]
    
    #Rename the tissue IDs with their donor/ recipient and cell type identity
    embryonic_median_cell_fracs<-generate_median_fracs_across_tissues(post_ind,selected_nodes=embryonic_branches)
    tissues=sapply(colnames(embryonic_median_cell_fracs),function(id) {
      bulk_smry_all%>%filter(tissueID==id)%>%mutate(new_id=paste(individual_type,cell_type,sep="_"))%>%pull(new_id)%>%.[1]
    })
    colnames(embryonic_median_cell_fracs)=tissues
    
    #(2) Get the names of the mutations that are on these embryonic branches
    embryonic_muts=rownames(embryonic_median_cell_fracs) #get a vector of mut_refs for these differing mutations
    nmut_max=1000
    cat(paste("There are",length(embryonic_muts),"embryonic mutations"),sep="\n")
    while(length(embryonic_muts)>nmut_max) {
      cat(paste("Reducing mutation set by reducing cut-off molecular time to",embryonic_mutation_cutoff-1),sep="\n")
      embryonic_mutation_cutoff<-embryonic_mutation_cutoff-1
      embryonic_branches=tree$edge[,2][which(nodeHeights(tree)[,2]<embryonic_mutation_cutoff)]
      embryonic_median_cell_fracs<-generate_median_fracs_across_tissues(post_ind,selected_nodes=embryonic_branches)
      embryonic_muts=rownames(embryonic_median_cell_fracs) #get a vector of mut_refs for these differing mutations
    }
    cat(paste("There are",length(embryonic_muts),"embryonic mutations included"),sep="\n")
    colnames(embryonic_median_cell_fracs)=tissues
    
    #(3) Calculate the phylogenetic distances for these mutations
    mut_pairs=as.data.frame(data.table::CJ(a=embryonic_muts,b=embryonic_muts,unique=FALSE,sorted=FALSE)) #create a df of all possible pairs
    mut_pairs$uid=apply(mut_pairs[,1:2],1,paste,collapse="-") #create a unique id for each pair
    cat(paste("This results in",nrow(mut_pairs),"mutation pairs"),sep="\n") #How many are included - gives an idea of how long it will take
    
    node_combs=as.data.frame(data.table::CJ(a=embryonic_branches,b=embryonic_branches,unique=TRUE))
    
    #Calculate genetic distances for all of these node pairs - save as a list
    cat("Calculating genetic distances between node pairs",sep="\n")
    nodeheights=nodeHeights(tree) #The lapply function needs the nodeheights object
    node_comb_dist=lapply(1:nrow(node_combs),function(i) {
      node_a=node_combs[i,1]
      node_b=node_combs[i,2]
      node_height_a=nodeheights[tree$edge[,2]==node_a,2]
      node_height_b=nodeheights[tree$edge[,2]==node_b,2]
      if(node_a==node_b) {stop(return(0))}
      ancestral_nodes_a=get_ancestral_nodes(node_a,tree$edge) #get all the ancestral nodes of each node & include the node itself
      ancestral_nodes_b=get_ancestral_nodes(node_b,tree$edge)
      common_ancestors=dplyr::intersect(c(node_a,ancestral_nodes_a),c(node_b,ancestral_nodes_b)) #pull out those that are common ancestors
      if(length(common_ancestors)==0) {
        genetic_dist <-(node_height_a+node_height_b) #If there are no common ancestors listed, then the closest common ancestor is the tree root & the genetic distance is the sum of the node heights
      } else {
        common_ancestors_heights=sapply(common_ancestors, function(x) {
          nodeheights[tree$edge[,2]==x,2] #otherwise find the heights of all common ancestors
        })
        genetic_dist = node_height_a+node_height_b-2*max(common_ancestors_heights) #the common ancestor with the maximum nodeheight is the most recent
      }
      return(genetic_dist)
    })
    
    #Name the list with the node pair in the format "node_1-node_2"
    node_comb_names=apply(node_combs,1,paste,collapse="-")
    names(node_comb_dist) <- node_comb_names
    
    #Now using this list as a reference, get the distance for all individual mutation pairs
    cat("Calculating the distance between mutation pairs",sep="\n")
    mut_pairs$dist<-NA
    #mut_pairs$dist[mut_pairs$uid %in% rownames(ref_dist)]=ref_dist[mut_pairs$uid[mut_pairs$uid %in% rownames(ref_dist)],"dist"]
    empties=which(is.na(mut_pairs$dist))
    for(i in empties) {
      mut_pairs$dist[i] <-shortest_dist(mut_a=mut_pairs[i,1],mut_b=mut_pairs[i,2],tree=tree,ref=post_ind[[1]][1:2],node_comb_dist=node_comb_dist)
      if(i%%1000==0) {print(i)}
    }
    
    ##(4) Convert the "distance" into a "similarity", by taking the inverse.
    mut_pairs$sim=1/(mut_pairs$dist+1) #Add one to denominator to avoid dividing by 0, and make max of 1.
    
    
    ##Apply this across all possible tissue pairs
    cat("Calculating the cosine similarity across all tissue pairs",sep="\n")
    tissue_pairs=as.data.frame(data.table::CJ(tissues,tissues))
    colnames(tissue_pairs) <- c("tissue_1","tissue_2")
    tissue_pairs$soft_cosim = NA
    for(i in 1:nrow(tissue_pairs)) {
      if(i%%10==0) {print(i)}
      tissue_pairs$soft_cosim[i] <-get_soft_cosim(data=t(embryonic_median_cell_fracs),selected_muts=embryonic_muts,tissue_1 = tissue_pairs[i,1],tissue_2=tissue_pairs[i,2],mut_pairs_sim = mut_pairs$sim)
    }
    
    #(6) Convert to wide format for the heatmap
    sim_mat<-tissue_pairs%>%
      tidyr::pivot_wider(names_from = tissue_2,values_from = soft_cosim) %>%
      dplyr::select(-1) %>%
      as.matrix()
    rownames(sim_mat) <- colnames(sim_mat)
    
    cat("Completed similarity matrix",sep="\n")
    similarity_matrices[[j]]<-sim_mat
  }
  
  saveRDS(object=similarity_matrices,file = similarity_matrices_file)
  
} else {
  similarity_matrices=readRDS(similarity_matrices_file)
}

#========================================#
# VISUALIZE THE RESULTS AS HEATMAPS OF SIMILARITY ####
#========================================#

#First bind the results into a single tidy dataframe for easy visualization with ggplot2
all_tidy_df<-Map(sim_mat=similarity_matrices,pair=names(all.trees.cc.nodups),function(sim_mat,pair) {
  cat(pair,sep="\n")
  if(!is.null(sim_mat)){
    #diag(sim_mat)<-NA
    tidy_df<-as.data.frame(sim_mat)%>%
      tibble::rownames_to_column(var = "Cell_type1")%>%
      tidyr::gather(-Cell_type1,key="Cell_type2",value="similarity")%>%
      mutate(Pair=pair,.before=1)
    return(tidy_df)
  } else {
    return(NULL)
  }
})%>%dplyr::bind_rows()


#--------------------------------------------------------------------------------------------------#
## -----BAYESIAN CLASSIFIER APPROACH (not phylo aware) ------
#--------------------------------------------------------------------------------------------------#

NV=matrices_test$NV[details_targ$mut_ref,]
NR=matrices_test$NR[details_targ$mut_ref,]
depth.cols=matrices_ctrl$NR[details_targ$mut_ref,]; mtr.cols=matrices_ctrl$NV[details_targ$mut_ref,] #The matrices for calculating background error (single cell colonies)

# Check matrix structure is identical
print(all(row.names(NR) == row.names(NV)))
print(all(names(NR) == names(NV)))
print(all(row.names(NR) == details_targ$mut_ref))
print(all(row.names(depth.cols) == row.names(NR)))
print(all(row.names(mtr.cols) == row.names(NR)))

# Set parameters
alpha_mut <- rep(1.5, nrow(NR))
beta_mut <- rep(10, nrow(NR))
min_theta <- 0.001
max_theta <- 0.1

# Estimate alpha and beta parameters for the error distribution
# Uses empirically weighted method of moments, as described by Keinman, JASA 1973
# Loop through each mutation -> find colonies not carrying that variant -> estimate alpha / beta from those colonies

alpha_error <- beta_error <- rep(0,nrow(NR))

for (i in 1:nrow(NR)) {
  
  curr.n <- unlist(depth.cols[i,])
  curr.x <- unlist(mtr.cols[i,])
  curr.x <- curr.x[curr.n > 0]
  curr.n <- curr.n[curr.n > 0]
  
  min_error <- max((sum(curr.x) + 1) / (sum(curr.n) + 2), 0.0005)
  moment.ests.1 <- moment.est(curr.x, curr.n, curr.n)
  if (moment.ests.1[1] < min_error) {
    p_hat <- min_error
    gamma_hat <- 0
  } else {
    new.wts <- curr.n / (1 + moment.ests.1[2] * (curr.n-1))
    moment.ests.2 <- moment.est(curr.x, curr.n, new.wts)
    p_hat <- moment.ests.2[1]
    gamma_hat <- moment.ests.2[2]
  }
  theta_hat <- min(max(gamma_hat / (1-gamma_hat), min_theta), max_theta)
  
  alpha_error[i] <- p_hat / theta_hat
  beta_error[i] <- (1 - p_hat) / theta_hat
}

# Define prior prob as fraction of colonies in the haematopoietic tree with the mutation
prior_mut <- sapply(1:nrow(NR), function(i) {max(length(Descendants(tree, details_targ$node[i], "tips")[[1]]) / length(tree$tip.label),0.01)})
names(alpha_error)=names(beta_error)=names(alpha_mut)=names(beta_mut)=names(prior_mut)=row.names(NR) #these need to be named for the clean_up_post_2 function

# Calculate posterior probabilities
post.prob <- cbind(sapply(1:ncol(NR), function(i) {bbprob.calculator(x_i = NV[,i], n_i = NR[,i], alpha_error = alpha_error, beta_error = beta_error, alpha_mut = alpha_mut, beta_mut = beta_mut, prior_mut = prior_mut)}))
colnames(post.prob) <- colnames(NR); row.names(post.prob) <- row.names(NR)


#--------------------------------------------------------------------------------------------------#
## -----Plot using Bayesian classifier approach (not phylogeny aware) ------
#--------------------------------------------------------------------------------------------------#


pdf(file=paste0(root_dir,"plots/Bayesian_classifier_plots_",MDIDs[ID],".pdf"),width=12,height=8)
temp=lapply(this_mouse_tissueIDs,function(sampleID) {
  dat<-readRDS(paste0(root_dir,"/GibbsSampler/data/",sampleID,"_gibbs_info.RDS"))
  details_targ=dat$details%>%tidyr::unite(col=mut_ref,Chrom,Pos,Ref,Alt,sep="-",remove=F)
  tree=squash_tree(dat$tree)
  
  details_targ<-details_targ%>%mutate(Mut_type=ifelse(nchar(Ref)==1 & nchar(Alt)==1,"SNV","indel"))
  generate_targ_seq_plots(samples=sampleID,
                          tree=tree,
                          details_targ=details_targ%>%filter(Mut_type=="SNV"),
                          matrices = matrices_test,
                          post.prob=post.prob,
                          info_type="cell_frac",
                          prob_threshold_to_include=0.1, #Probability threshold from the post.prob matrix for plotting
                          plot_cell_frac=TRUE,
                          plot_donut=TRUE,
                          donut_info="cell_frac", #other option is "lineages_lost"
                          CI=0.95,  #Confidence intervals on the pie chart, default = 80% CI
                          radius=2.8,  #Radius of the pie charts on the plot
                          scale_muts_to_branch=T,
                          colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                          vaf_lwd=5,
                          overlay=F,
                          title=metadata%>%filter(INTERNAL_CASM_SAMPLE_NAME==sampleID)%>%pull(SUPPLIER_SAMPLE_ID)) 
  
})
dev.off()


