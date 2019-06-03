out <- 1.2


classify <- function(data){
  k <- dist(t(data))
  k <- hclust(k,method = "ward.D2")
  k <- cutree(k,k = 2)
  k$cluster <- k
  #Create the definition of UC1 (lower expression) and UC2 (higher expression)
  if( mean(colMeans(data)[k$cluster==2]) < mean(colMeans(data)[k$cluster==1]) ){ UC_subtype <- factor(paste0("UC",k$cluster)) } else {UC_subtype <- factor(ifelse(k$cluster==1,"UC2","UC1"))}
 return(UC_subtype)
}


my_violins <- function(Data,group,main="",color,xlab="",ylab="",ylim=c(min(Data), max(Data)*1.1)){
#boxplot(Data ~ group,ylim=c(min(Data), max(Data)*1.1),las=1,xlab=xlab,ylab=ylab,main=main,border="white")
plot(Data ~ group,ylim=ylim,las=1,xlab=xlab,ylab=ylab,main=main,border="white",type="n",yaxs="i")
for (i in 1:length(levels(group))){
  a = Data[group == levels(group)[i]]
  vioplot(at = i, a ,add = T, col = paste0(color[i],"80"),drawRect = F,border = F,range = c(0,15))
  a = Data[group == levels(group)[i]]
  lines(c(i-0.3,i+0.3),c(median(a),median(a)),lwd=2)}
  beeswarm(Data~group,col=color[1:length(levels(group))],pch=16,cex=0.6,add=T,priority = "density")}



plot_series <- function(data){
  mypar(1,5)
  PC <- prcomp(apply(data,1,function(x)scale(x,center = T,scale = T)))
plot(PC$x[,1],PC$x[,2],col="grey60",pch=16,xlim=c(min(PC$x[,1])*out,max(PC$x[,1])*out),ylim=c(min(PC$x[,2])*out,max(PC$x[,2])*out),las=1,xlab="PC1",ylab="PC2",cex=1.5,main="PCA")

set.seed(123)
tSNE <- Rtsne(t(data),perplexity = 30,pca = T,max_iter = 1000,theta = 0.0)
plot(tSNE$Y[,1],tSNE$Y[,2],cex=1.5,las=1,xaxt="n",yaxt="n",col="grey60",pch=16,line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out),main="tSNE")
axis(side = 1, labels = FALSE, tck = -0.01);axis(side = 2, labels = FALSE, tck = -0.01)


z <- kde2d(tSNE$Y[,1], tSNE$Y[,2], n=100, lims = c(min(tSNE$Y[,1])*out, max(tSNE$Y[,1])*out, min(tSNE$Y[,2])*out, max(tSNE$Y[,2])*out))

plot(tSNE$Y[,1],tSNE$Y[,2],type="n",xaxt="n",yaxt="n",line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out),main="contour tSNE")
contour(z, drawlabels=FALSE, nlevels=5, col=c("red","grey","blue2","green4","orange3","firebrick"), add=TRUE)
points(tSNE$Y[,1], tSNE$Y[,2],cex=.5,pch=16)

UC_subtype <- classify(data)

plot(tSNE$Y[,1],tSNE$Y[,2],cex=1.5,las=1,xaxt="n",yaxt="n",col=pallete[UC_subtype],line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out),main="tSNE per cluster",pch=as.numeric(UC_subtype)+15)
axis(side = 1, labels = FALSE, tck = -0.01);axis(side = 2, labels = FALSE, tck = -0.01)

for(i in 1:length(unique(UC_subtype))){
  temp <- levels(UC_subtype)[i]
  for (j in (1:length(UC_subtype))[UC_subtype == temp]){
    lines(c(tSNE$Y[j,1],mean(tSNE$Y[UC_subtype == temp,1])),c(tSNE$Y[j,2],mean(tSNE$Y[UC_subtype == temp,2])),col=paste0(pallete[i],"90"),lwd=0.5)}}

#points(k$centers,pch=22,bg="black",col=pallete[1:ncol(k$centers)],cex=2)
a <- dataEllipse(tSNE$Y[,1], tSNE$Y[,2],draw = F,levels = 0.8,plot.points = F,group.labels = NA,groups = factor(UC_subtype))
for(i in 1:2){
  polygon(a[[i]],lty=1,border=pallete[i],lwd=2)}

model <- glm((UC_subtype)~tSNE$Y,family="binomial")
slope <- coef(model)[2]/(-coef(model)[3])
intercept <- coef(model)[1]/(-coef(model)[3]) 
abline(intercept,slope,lty=2)

grad <- (colMedians(data) - min(colMedians(data))) / (max(colMedians(data)) - min(colMedians(data)))
plot(tSNE$Y[,1],tSNE$Y[,2],cex=1.5,las=1,xaxt="n",yaxt="n",col=colorRampPalette(c("grey90","grey80","navy"))(50)[grad*49+1],pch=as.numeric(UC_subtype)+15,line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out),main="tSNE median expr.")
axis(side = 1, labels = FALSE, tck = -0.01);axis(side = 2, labels = FALSE, tck = -0.01)
abline(intercept,slope,lty=2)

}