library(ggplot2)
library(Cairo)

df <- read.table('nxrepair.csv',col.names=c("scaffold","pos","zscore"),as.is=T)

# I want nice headings for each subplot
pl <- unlist(strsplit(df[,1],"_"))
df.nrows=dim(df)[1]
sf <- rep("scaffold",df.nrows)
idx <- seq(2,length(pl),by=8)
nums <- pl[idx]
pp <- paste(sf,nums,sep="_")
# Overwrite old scaffold names
df[,1] <- as.factor(pp)

CairoPNG(file="zscores.per.scaffold.png",width=600,height=600)
#png(file="zscores.per.scaffold.png",width=600,height=600)
ggplot(df,aes(x=pos,y=zscore)) + geom_line() + facet_wrap(~scaffold,scales="free")
dev.off()
