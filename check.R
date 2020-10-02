
list<-list()
for(i in c(1:4)){
  file=paste("./", i, "/density.csv", sep="")
  ex<-as.data.frame(file.exists(file))
  ex$perm<-paste(i)
  list[[length(list)+1]]= ex
}
exs<-do.call(rbind.data.frame, list)
non<-exs[(exs$`file.exists(file)`=="FALSE"),
write.csv(non, file="missing.csv", quote)
