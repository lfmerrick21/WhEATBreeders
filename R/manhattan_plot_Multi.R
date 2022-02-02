manhattan_plot_Multi <- function(GWAS,cutoff=0.05,QTN_index2 = NULL,QTN_index3=NULL,ZZ_label=NULL,ZZ_circle=NULL,Sig_L=NULL,Sig_Multi=NULL, trait = NULL,labels=NULL,model="unknown",Nchr=22,Multi=NULL,Multi_labels=NULL,correction="FDR",model_cor="BLINK",Third_Labels=NULL,Trait_one="Trait 1",Trait_two="Trait 2",Comparison="Population",Multi2=NULL,Multi2_labels=NULL){

  #library(compiler)
  #source("http://zzlab.net/GAPIT/GAPIT.library.R")
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  library(GAPIT3)
  #find cutoff values using FDR or Bonferonni Cutfoff
  if(!is.null(Sig_L)){
    if(model_cor=="BLINK"){
      if(correction=="FDR"){
        cutoff.final=max(Sig_L$P.value)
      }else{
        cutoff.final=cutoff/length(GWAS$P.value)
      }}

    if(model_cor=="MLM"){
      if(correction=="FDR"){
        cutoff.final=max(Sig_L$P.value)
      }else{
        cutoff.final=cutoff/length(GWAS$P.value)
      }}

    SG_pval <- -log10(cutoff.final)
  }

  if(!is.null(Sig_Multi)){
    if(model_cor=="BLINK"){
      if(correction=="FDR"){
        cutoff.final.Multi=max(Sig_Multi$P.value)
      }else{
        cutoff.final.Multi=cutoff/length(Multi$P.value)
      }}

    if(model_cor=="MLM"){
      if(correction=="FDR"){
        cutoff.final.Multi=max(Sig_Multi$P.value)
      }else{
        cutoff.final.Multi=cutoff/length(Multi$P.value)
      }}

    SG_pval_Multi <- -log10(cutoff.final.Multi)
  }

  if(is.null(Sig_L)){
    cutoff.final=cutoff/length(GWAS$P.value)
    SG_pval <- -log10(cutoff.final)
  }
  if(is.null(Sig_Multi)){
    cutoff.final.Multi=cutoff/length(Multi$P.value)
    SG_pval_Multi <- -log10(cutoff.final.Multi)

  }


  col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(Nchr/5))

  #GWAS$Chromosome=as.character(GWAS$Chromosome)
  #GWAS$Position=as.numeric(GWAS$Chromosome)

  don <- GWAS %>%
    # Compute chromosome size
    group_by(Chromosome) %>%
    summarise(chr_len=max(Position)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(GWAS, ., by=c("Chromosome"="Chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot)

  if(!is.null(labels)){
    don=don%>% mutate( is_annotate=ifelse(SNP %in% labels, "yes", "no"))
    don=don%>%mutate( is_highlight=ifelse(P.value<=cutoff.final, "yes", "no"))
    SNP_don= don %>% filter(SNP %in% labels$SNP)
  }else{
    don$is_annotate=NA
    don$is_highlight=NA
  }

  #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  don$Trait=rep(Trait_one,length(don$SNP))


  if(!is.null(Multi)){
    #Multi$Chromosome=as.character(Multi$Chromosome)
    #Multi$Position=as.numeric(Multi$Chromosome)
    don_Multi <- Multi %>%
      # Compute chromosome size
      group_by(Chromosome) %>%
      summarise(chr_len=max(Position)) %>%
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(Multi, ., by=c("Chromosome"="Chromosome")) %>%
      # Add a cumulative position of each SNP
      arrange(Chromosome, Position) %>%
      mutate( BPcum=Position+tot)
    if(!is.null(Multi_labels)){
      don_Multi=don_Multi%>% mutate( is_annotate=ifelse(SNP %in% Multi_labels$SNP, "yes", "no"))
      don_Multi=don_Multi%>%mutate( is_highlight=ifelse(P.value<=cutoff.final.Multi, "yes", "no"))
      SNP_don_Multi= don_Multi %>% filter(SNP %in% Multi_labels$SNP)

    }else{
      don_Multi$is_annotate=NA
      don_Multi$is_highlight=NA
    }
  }

  if(!is.null(Multi)){
    #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
    axisdf = don_Multi %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    don_Multi$Trait=rep(Trait_two,length(don_Multi$SNP))


    don_Both=rbind(don,don_Multi)
    don_Both$Chr_Trait=paste(don_Both$Chromosome,don_Both$Trait)
    #don_Both=mutate(don_Both,Chr_Trait,paste(chromosome,Trait))
  }

  if(!is.null(Multi_labels)){
    SNP_both=rbind(SNP_don,SNP_don_Multi)
    both_labels=intersect(labels$SNP,Multi_labels$SNP)
    SNP_Multi= don_Both %>% filter(SNP %in% both_labels)
  }else{
    SNP_Multi=data.frame()
  }

  if(is.null(Multi)){

    manhattan_plot <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
      # Show all points
      geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=1.3) +
      scale_color_manual(values = col.Oceanic) +

      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +
      scale_shape_manual(values=c(21))+# remove space between plot area and x axis
      #scales::number_format(accuracy = 0)
      #ylim(0,  max(-log10(don$P.value)))+
      # Custom the theme:
      geom_hline(yintercept= SG_pval,color="black")+
      labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
        y = "-log10(p)",
        x = "Chromosome",
        color = "Chromosome",
        tag=model)+
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag=element_text(angle=-90,size=10),
        plot.tag.position=c("right")
      )
    if(!is.null(QTN_index3)){
      QTN=don%>%filter(SNP %in% QTN_index3)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red",linetype="solid") }
    if(!is.null(QTN_index2)){
      QTN=don%>%filter(SNP %in% QTN_index2)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "green",linetype="dashed") }
    if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
    if(!is.null(labels)){
      manhattan_plot=manhattan_plot+
        # Add highlighted points
        geom_point(data=subset(don, is_highlight=="yes"),size=3,shape=21) +
        # Add label using ggrepel to avoid overlapping
        geom_text( data=SNP_don, aes(label=SNP, size=3), nudge_y = 5)
    }

    if(!is.null(Third_Labels)){

      don_both_third=don%>%filter(SNP %in% Third_Labels)
      manhattan_plot=manhattan_plot+ geom_label_repel(data=don_both_third, aes(label=SNP), size=3) }


    if(!is.null(ZZ_label)){

      ZZ=don%>%filter(SNP %in% ZZ_label)
      manhattan_plot=manhattan_plot+ geom_text_repel( data=ZZ, aes(label=SNP, size=3,y = 10))  }

  }else{
    col.dual=rep(c(  '#194CEC','#B10913', '#349EF7','#F5414C'),ceiling(Nchr*2/4))
    col.dual.shape=rep(c(1,2),ceiling(Nchr*2/2))
    levels(as.factor(don_Both$Chr_Trait))
    manhattan_plot <- ggplot(don_Both, aes(x=BPcum, y=-log10(P.value))) +
      # Show all points
      geom_point( aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), alpha=0.8, size=1.3) +
      scale_color_manual(values = col.dual) +
      scale_shape_manual(values=c(21, 24))+
      #scale_shape_manual(values=col.dual.shape)+

      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +     # remove space between plot area and x axis
      #scales::number_format(accuracy = 0)
      #ylim(0,  max(-log10(don$P.value)))+
      # Custom the theme:
      geom_hline(yintercept= SG_pval,color="black")+
      geom_hline(yintercept= SG_pval_Multi,color="black",linetype="dashed")+

      labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
        y = "-log10(p)",
        x = "Chromosome",
        color = Comparison,
        tag=model)+
      theme_bw() +
      theme(
        #legend.position=c(1,1), legend.justification=c(1,1),legend.background=element_rect(fill="white", colour="black"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag=element_text(angle=-90,size=10),
        plot.tag.position=c("right"))
    if(!is.null(QTN_index3)){
      QTN=don%>%filter(SNP %in% QTN_index3)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red",linetype="solid") }
    if(!is.null(QTN_index2)){
      QTN=don%>%filter(SNP %in% QTN_index2)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "green",linetype="dashed") }
    if(!is.null(ZZ_label)){
      ZZ=don%>%filter(SNP %in% ZZ_label)
      manhattan_plot=manhattan_plot+ geom_text_repel( data=ZZ, aes(label=SNP,y = 10)) }
    if(!is.null(ZZ_circle)){
      ZZ=don%>%filter(SNP %in% ZZ_circle)
      manhattan_plot=manhattan_plot+ geom_point(data=ZZ,shape=1, size=5,color="green") }
    if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
    if(!is.null(labels)){
      if(!is.null(Multi_labels)){
        manhattan_plot=manhattan_plot+
          # Add highlighted points
          geom_point(data=subset(don_Both, is_highlight=="yes"),aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), size=3) +
          # Add label using ggrepel to avoid overlapping
          geom_text_repel( data=SNP_both, aes(label=SNP, size=3))
      }else{
        manhattan_plot=manhattan_plot+
          # Add highlighted points
          geom_point(data=subset(don_Both, is_highlight=="yes"),aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), size=3) +
          # Add label using ggrepel to avoid overlapping
          geom_text( data=SNP_don, aes(label=SNP, size=3), nudge_y = 5)
      }

      # Add highlighted points
      #geom_point(data=subset(don, is_highlight=="yes"), size=3) +
      # Add label using ggrepel to avoid overlapping
      #geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3)+

      #geom_point(data=subset(don_Multi, is_highlight=="yes"), size=3) +
      # Add label using ggrepel to avoid overlapping
      #geom_text_repel( data=subset(don_Multi, is_annotate=="yes"), aes(label=SNP), size=3)
    }
    if(!is.null(Third_Labels)){

      don_both_third=don%>%filter(SNP %in% Third_Labels)
      manhattan_plot=manhattan_plot+ geom_label_repel(data=don_both_third, aes(label=SNP), size=3) }
    if(!nrow(SNP_Multi)==0){manhattan_plot=manhattan_plot+geom_point(data=SNP_Multi,shape=1, size=5,color="green")


    }
  }

  return(manhattan_plot)
}
