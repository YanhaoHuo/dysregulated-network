regulateDS<-function(BL_regulate,d){
  mean_c<-c()
  for(j in 1:length(d)){
    w<-which(BL_regulate$Gene1==d[j])#|BL_regulate$Gene2==drivergenelist[[i]][j])#gene1????????????gene2
    if(length(w)!=0){mean_c<-c(mean_c,sum(abs(BL_regulate$logFC[w])))}#abs logFC
    else{mean_c<-c(mean_c,0)}
  }
  DS<-data.frame(DS=mean_c[order(-mean_c)])
  rownames(DS)<-d[order(-mean_c)]
  return(DS)
}
coverDS<-function(BL_regulate,ds){
  d<-rownames(ds)
  cover_gene<-c()
  cover_gene_length<-c()
  cover_gene_pre<-c()
  for(i in 1:length(d)){
    w<-which(BL_regulate$Gene1==d[i])
    if(length(w)!=0){
      cover_gene<-union(cover_gene,BL_regulate$Gene2[w])
      cover_gene<-union(cover_gene,d[i])
      cover_gene_length<-c(cover_gene_length,length(cover_gene))
      cover_gene_pre<-c(cover_gene_pre,length(cover_gene)/length(d))
    }else{
      cover_gene<-union(cover_gene,cover_gene)
      cover_gene<-union(cover_gene,d[i])
      cover_gene_length<-c(cover_gene_length,length(cover_gene))
      cover_gene_pre<-c(cover_gene_pre,length(cover_gene)/length(d))
      }
  }
  ds$cover=cover_gene_length
  ds$percentage=cover_gene_pre
  return(ds)
}
