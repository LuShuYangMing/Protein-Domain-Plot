#library(UniProt.ws)
#up <- UniProt.ws(taxId=83332)
#pkng <- select(up, keys=c("P9WI73"), columns=c("REFSEQ_NUCLEOTIDE","REFSEQ_PROTEIN","LENGTH","DATABASE(PFAM)"), keytype="UNIPROTKB")
library(httr)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(ggthemes)

get_json <- function(protein){
  # accession numbers need to be combined with "%2C" for Uniprot API
  proteins_acc_url <- gsub(" ", "%2C", proteins_acc)
  baseurl<-"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
  # make url to GET features for multiple accession numbers
  url <- paste0(baseurl, proteins_acc_url)
  # accept_json() argument gives a JSON object
  prots_feat <- GET(url, accept_json())
  code <- status_code(prots_feat)  
  # if it returns a 200 - that's good
  if(code == 200){
    print("Download has worked")
  } else {
    print(paste("An error has occured. Code:", code))
  }
  protein_json <- content(prots_feat)
  return(protein_json)
}

extract_names <- function(protein_json){
  # Create a variable of protein_json[[1]] to prevent repitition
  prot_info <- protein_json[[1]]
  # Extract list of names...
  names <- list(
    accession = prot_info$accession,
    name = prot_info$id,
    protein.recommendedName.fullName =
      prot_info$protein$recommendedName$fullName$value,
    gene.name.primary = prot_info$gene[[1]]$name$value,
    gene.name.synonym = prot_info$gene[[1]]$synonyms[[1]]$value,
    organism.name.scientific = prot_info$organism$names[[1]]$value
  )
  return(names)
}

get_features <- function(protein){
  feat_api_url <- c("https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=")
  comb_acc_api <- paste0(feat_api_url, protein)
  # basic function is GET() which accesses the API
  prot_feat <- GET(comb_acc_api, accept_json())
  code <- status_code(prot_feat)  
  # returns a 200 is good
  if(code==200){
    print("Download has worked")
  }else{
    print(paste("An error has occurred. Code:", code))
  }
  prot_feat %>%
    content() %>%
    flatten() ->
    features_list
  # create the data.frame object called features
  features <- NULL
  for(i in 1:length(features_list$features)){
    if(is.null(features_list$features[[i]]$description) == TRUE){
      Type <- c(features_list$features[[i]]$type,
                        "NONE",
                        as.numeric(features_list$features[[i]]$begin),
                        as.numeric(features_list$features[[i]]$end))
    } else{
      Type<- c(features_list$features[[i]]$type,
                        as.character(features_list$features[[i]]$description),
                        as.numeric(features_list$features[[i]]$begin),
                        as.numeric(features_list$features[[i]]$end))
    }
    features <- rbind(features, Type)
  }
  features_dataframe <- as.data.frame(features, stringsAsFactors = FALSE)
  features_dataframe$order <- 0
  colnames(features_dataframe) <- c("type", "description", "begin", "end", "order")
  features_dataframe$begin <- as.numeric(features_dataframe$begin)
  features_dataframe$end <- as.numeric(features_dataframe$end)
  features_dataframe$length <-
    features_dataframe$end - features_dataframe$begin
  
  # add accession number to each row of dataframe
  features_dataframe$accession <- rep(features_list$accession,
                                      times = nrow(features_dataframe))
  
  # add entryName (e.g. p65_HUMAN) to each row of dataframe
  features_dataframe$entryName <- rep(features_list$entryName,
                                      times = nrow(features_dataframe))
  
  # add taxid to each row of datafame
  features_dataframe$taxid <- rep(features_list$taxid,
                                  times = nrow(features_dataframe))
  
  return(list(features_dataframe,features_list$sequence))
}

draw_canvas <- function(data = data){
  begin = end = NULL
  p <- ggplot()
  p <- p + labs(x = "Amino acid") # label x-axis
  p <- p + labs(y = "Number of mutations") # label y-axis
  return(p)
}

draw_chains <- function(p, data = data, outline = "#FFCA05", fill = "#FFCA05", label_chains = TRUE, 
                        labels = data[data$type == "CHAIN",]$entryName, size = 0.5, label_size = 4, ... ,expand_value = 0.2){
  begin = end = NULL
  data$expand_value <- expand_value
  p <- p + geom_rect(data = data[data$type == "CHAIN",], mapping= aes(xmin = begin,xmax = end,
                                                   ymin = order - expand_value, ymax = order + expand_value),
                     colour = outline, fill = fill, size = size)
  if(label_chains == TRUE){
    p <- p + annotate("text", x = data[data$type == "CHAIN",]$length+100, y = data[data$type == "CHAIN",]$order,
                        label = labels, hjust = 1, size = label_size)
  }
  return(p)
}

draw_domains <- function(p,color=color,data = data, label_domains = TRUE, label_size = 4, ... , expand_value = 0.25, show.legend=FALSE){
  arg <- match.call()
  begin = end = NULL
  data$expand_value <- expand_value
  data$domaincolor <- "white"
  data[data$type == "DOMAIN","domaincolor"] <- color
  p <- p + geom_rect(data = data[data$type == "DOMAIN",], mapping = aes(xmin = begin, xmax = end, ymin = order-expand_value,
                                                                        ymax = order + expand_value,colour=NA,fill=domaincolor),
                     show.legend = show.legend)+scale_fill_identity()
  if(label_domains == TRUE){
    p <- p + geom_text(data = data[data$type == "DOMAIN", ], aes(x = begin + (end-begin)/2, y = order, label = description), 
                       size = label_size,color=NA,fill=domaincolor)
  }
  return(p)
}

draw_phospho <- function(p, data = data, size = 2, fill = "yellow", expand_value = 0.25){
  begin = end = NULL
  data$expand_value <- expand_value
  features <- data[data$type == "MOD_RES",]
  phospho_list <- grep("Phospho", features$description)
  phospho_features <- features[phospho_list,]
  p <- p + geom_point(data = phospho_features, aes(x = begin, y = order+expand_value),
                               shape = 21, colour = "black", fill = fill,
                               size = size)
  return(p)
}

draw_motif <- function(p, data = data, expand_value = 0.25){
  begin = end = NULL
  data$expand_value <- expand_value
  p <- p + geom_rect(data= data[data$type == "MOTIF",], mapping = aes(xmin = begin, xmax = end, ymin=order - expand_value,
                                                                      ymax=order + expand_value, fill=description))
  return(p)
}

draw_repeat <- function(p, data = data, label_size = 4, outline = "dimgrey", 
                        fill = "dimgrey", label_repeats = TRUE, expand_value = 0.25){
  begin = end = NULL
  data$expand_value <- expand_value
  p <- p + geom_rect(data= data[data$type == "REPEAT",], mapping= aes(xmin = begin, xmax = end, ymin=order - expand_value, 
                                                                      ymax=order + expand_value), colour = outline, fill = fill)
  if(label_repeats == TRUE){
    # label repeats (for this they are ANK but remove digits)
    p <- p + geom_text(data = data[data$type == "REPEAT",], aes(x = begin + (end-begin)/2, y = order,
                                                                label = gsub("\\d", "", description)), size = label_size)
  }
  return(p)
}

draw_recept_dom <- function(p, data = data, label_domains = FALSE, label_size = 4, expand_value = 0.25){
  begin = end = NULL
  data$expand_value <- expand_value
  p <- p +geom_rect(data= data[data$type == "TOPO_DOM",], mapping = aes(xmin = begin, xmax = end, ymin = order - expand_value,
                                                                        ymax=order + expand_value, fill = description))
  p <- p + geom_rect(data= data[data$type == "TRANSMEM",], mapping = aes(xmin = begin, xmax = end, ymin = order - expand_value,
                                                                         ymax=order + expand_value, fill = description))
  if(label_domains == TRUE){
    p <- p + geom_label(data = data[data$type == "TOPO_DOM", ], aes(x = begin + (end-begin)/2, y = order, label = description), size = label_size)
    p <- p + geom_label(data = data[data$type == "TRANSMEM", ], aes(x = begin + (end-begin)/2, y = order, label = "TM"), size = label_size)
  }
  return(p)
}

cal_number <- function(data){
  data %>% 
    group_by(Codon) %>%
    mutate( y= row_number()) %>%
    separate(col= Amino, c("Old_Amino", "new_Amino"), sep="/" ) ->
    df
  df$text <- paste(df$Old_Amino,df$Codon,df$new_Amino,sep="")
  df <- df[,c("Codon", "y", "text", "Effect")]
  df$Effect <- as.factor(df$Effect)
  return(df)
}

draw_snp <- function(p, data, point_size = 1 ,text_point= TRUE, text_size = 2){
  p <- p + geom_segment(data = data, aes(x = Codon, y = y, yend = 0, xend = Codon), colour = "grey")
  p <- p + geom_point(data = data ,aes(x = Codon, y = y, colour=Effect), size = point_size)#+
    #scale_color_manual(name = "Syn&Nonsyn", values = c("Syn"="#3CB200", "Nonsyn"="red"))
  if (text_point == TRUE){
    p <- p + geom_text(data = data, mapping = aes(x = Codon, y = y, label=text), size=text_size, check_overlap= T)
  }
  return(p)
}

get_xlabel <- function(sequence_data, point, normal=TRUE){
  new  <- sequence_data[sequence_data$index%in%point,]
  seq <- c(1,seq(0,nrow(sequence_data),200)[-1])
  if (!(nrow(sequence_data) %in% seq)){
    seq <- c(seq, nrow(sequence_data))
  }
  if (normal){
    return(seq)
  }
  if (!all(new$index %in% seq)){
    for (i in new$index){
      if (i %in% seq){
        seq <- seq[-which(seq==i)]
      }
    }
  }
  rest <- data.frame(index=seq,sequence=NA,label=seq)
  new <- rbind(new,rest)
  new <- new[order(new$index),]
  xlabel <- new[,c("index","label")]
  return(xlabel)
}

draw_theme <- function(p, data,xlabel ,max_snp_number=100, negetive_value=-1, ymin=-1){
  p <- p + theme_bw()
  #a = data.frame(max_snp_number=max_snp_number, max=max(data$end, na.rm=TRUE))
  #p <- p + geom_segment(data = a, mapping = aes(y = 0, yend = max_snp_number, x = -max*0.2, xend = -max*0.2), lwd = 1)
  #b = data.frame(negetive_value=negetive_value,max=max(data$end))
  #p <- p + geom_segment(data = b, mapping = aes(y = negetive_value, yend = negetive_value, x = 1, xend= max), lwd = 1)
  if (is.vector(xlabel)){
    p <- p + scale_x_continuous(breaks = xlabel, labels = xlabel)
  }else{
    p <- p + scale_x_continuous(breaks = xlabel$index, labels = xlabel$label)
  }
  #p <- p + expand_limits(x = 0)
  #p <- p + xlim(-max(data$end, na.rm=TRUE)*0.2,
  #              max(data$end, na.rm=TRUE) + max(data$end, na.rm=TRUE)*0.1)
  p <- p + ylim(ymin,max_snp_number)
  oo <- ggplot_build(p)
  xrange <- range(oo$layout$panel_ranges[[1]]$x.major_source)
  yrange <- range(oo$layout$panel_ranges[[1]]$y.major_source)
  a = data.frame(x1=xrange[1], x2=xrange[2],y1=yrange[1], y2 =yrange[2])
  p <- p + geom_segment(data = a, aes(x=x1, xend= x2, y= -Inf, yend= -Inf))
  p <- p + geom_segment(data = a, aes(y=y1, yend= y2, x= -Inf, xend= -Inf))
  p <- p + theme_classic() + theme(axis.line = element_blank(),
                                   axis.ticks.length = unit(5,"pt"))
  return(p)
}

#protein_json <- get_features("P9WI73")
# use content() function from httr to give us a list gives a Large list
# with 14 primary parts and lots of bits inside
#names <- extract_names(protein_json)



feature_list <- get_features("Q79FP0")
sequence <- feature_list[[2]]
sequence_frame <- data.frame(index=c(1:length(strsplit(sequence,"")[[1]])),sequence=strsplit(sequence,"")[[1]])
sequence_frame$label <- paste(sequence_frame$sequence, sequence_frame$index,sep = "")
feature <- feature_list[[1]]
features_plot <- feature[feature$length > 0, ]
#features_plot <- features_plot[complete.cases(features_plot),]
#features_plot <- features_plot[features_plot$type!= "CONFLICT",]
library(readxl)

setwd("E:\\Institute_of_Microbiology\\Mycobacterium_tuberculosis\\Rv1468c\\SNP\\")
snp <- read_xlsx("gmtv-Rv1468c.xlsx")
wanted <- snp[,c("Codon number","Amin acid", "Effect")]
colnames(wanted) <- c("Codon","Amino","Effect")

want2 <- wanted[as.character(wanted$Effect)=="Nonsyn",]

snp <- cal_number(wanted)

#snp2 <- cal_number(want2)



labels <- get_xlabel(sequence_data = sequence_frame, point = c(32,66,93), normal = F)
features_plot$order <- rep(-2.1,3)
#special for UBA domain
features_plot["type.3",] <- c("DOMAIN","UBA",32,66,-2.1,34,"Q79FP0","Q79FP0_MYCTU","83332")
features_plot$order <- as.numeric(features_plot$order)
features_plot$begin <- as.numeric(features_plot$begin)
features_plot$end <- as.numeric(features_plot$end)
features_plot$length <- as.numeric(features_plot$length)
p <- draw_canvas(features_plot)
p <- draw_snp(p ,snp, text_point = F,point_size = 1.2)
p <- draw_chains(p, features_plot, expand_value = 2)
p <- draw_domains(p,color=c("#6699CC","#996633"), data=features_plot,label_domains = F, label_size =3, expand_value = 2.1)
p <- draw_repeat(p, features_plot, label_size =2.5, expand_value = 1.1)
p <- draw_theme(p, features_plot,labels,ymin = -4.5)
p+theme(axis.line = element_line(size = 1),axis.ticks = element_line(size=1),
        axis.text = element_text(face = "bold"),axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),legend.title = element_text(face = "bold"))
p
#p <- draw_canvas(features_plot)
#p <- draw_snp(p ,snp2, text_point = F)
#p <- draw_chains(p, features_plot, expand_value = 0.75)
#p <- draw_domains(p, features_plot, label_size =3, expand_value = 0.9)
#p <- draw_repeat(p, features_plot, label_size =2.5, expand_value = 0.9)
#p <- draw_theme(p, features_plot,xlabel = labels ,max_snp_number = 35)
#p 

