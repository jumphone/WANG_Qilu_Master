
output_dir <- "/home/zhangfeng/projects/QRLu/data/12.5/"
setwd(output_dir)

library(Seurat)
library(dplyr)
library(future)
getOption("future.globals.maxSize")
options(future.globals.maxSize = 100 * 1024^3)


seurat_obj=readRDS('seurat_obj_umap.rds')

xy=seurat_obj@images$image@coordinates
x=xy[,1]
y=xy[,2]
mat=seurat_obj@assays$SCT@data


draw_exp <-function(x,y,z, output){
    library(circlize) 
    this_col_fun=colorRamp2(c(min(z), max(z)), 
                            c("grey90","red"))
    this_col=this_col_fun(z)
    scale_bar_z=seq(from=min(z),to=max(z),length.out=1000)
    scale_bar_col=this_col_fun(scale_bar_z)
    pdf(output)
    plot(x,y,pch=16,cex=0.1,col=this_col)
    plot(scale_bar_z,rep(1,length(scale_bar_z)),col=scale_bar_col,type='h',lwd=2)
    dev.off()
    }

z=mat['SOX11',]
draw_exp(x,y,z,output='plot_exp_SOX11.pdf')

z=mat['HNRNPH1',]
draw_exp(x,y,z,output='plot_exp_HNRNPH1.pdf')



set.seed(123)
km=kmeans(cbind(x,y),centers=60)
clst=km$cluster

pdf('plot_kmeans.pdf')
plot(x,y,cex=0.1,pch=16,col=clst)
text(x=km$centers[,1],y=km$centers[,2],labels=rownames(km$centers),cex=1.2,col='black')
text(x=km$centers[,1],y=km$centers[,2],labels=rownames(km$centers),cex=1,col='white')
text(x=km$centers[,1],y=km$centers[,2],labels=rownames(km$centers),cex=0.8,col='black')
dev.off()

saveRDS(km,'km.rds')

label=rep(0,length(x))
label[which(clst %in% c(42))]=1

this_col=rep('grey70',length(x))
this_col[which(label==1)]='red'
pdf('plot_select_cell.pdf')
plot(x,y,cex=0.1,pch=16,col=this_col)
dev.off()



z=mat['SOX11',]
z1=z[which(label==1)]
z0=z[which(label==0)]

pdf('plot_den_SOX11.pdf')

dz1=density(z1,bw=0.2)
dz0=density(z0,bw=0.2)
#dz1_y=dz1$y/max(dz1$y)
#dz0_y=dz0$y/max(dz0$y)
minx=0 #min(dz1$x,dz0$x)
maxx=2.5#max(dz1$x,dz0$x)
miny=min(dz1$y,dz0$y)
maxy=max(dz1$y,dz0$y)

plot(dz1$x, dz1$y, xlim=c(minx,maxx),ylim=c(0,maxy),col='red',type='l',lwd=3)
points(dz0$x, dz0$y, xlim=c(minx,maxx),ylim=c(0,maxy),col='grey70',type='l',lwd=3)

abline(v=mean(z1),lwd=3,col='red')
print(mean(z1))
#0.7571536
abline(v=mean(z0),lwd=3,col='grey79')
print(mean(z0))
#0.639591
dev.off()




pdf('plot_box_SOX11.pdf')
boxplot(z1,z0,col=c('red','grey70'),pch='+',lwd=2)
dev.off()
print(wilcox.test(z1,z0))
#p-value < 2.2e-16
print(t.test(z1,z0))
#p-value =4.035e-16



z=mat['HNRNPH1',]
z1=z[which(label==1)]
z0=z[which(label==0)]

pdf('plot_box_HNRNPH1.pdf')
boxplot(z1,z0,col=c('red','grey70'),pch='+',lwd=2)
dev.off()
print(wilcox.test(z1,z0))
#p-value < 2.2e-16
print(t.test(z1,z0))
#p-value < 2.2e-16


mat1=mat[,which(label==1)]
mat2=mat[,which(label==0)]
mean1=apply(mat1,1,mean)
mean2=apply(mat2,1,mean)

bkg_fc= (mean1 +0.1)/(mean2+0.1)

saveRDS(mean1,file='mean1.rds')
saveRDS(mean2,file='mean2.rds')
saveRDS(bkg_fc,file='bkg_fc.rds')

print(bkg_fc['SOX11'])
#1.158956
ecdf(bkg_fc)(bkg_fc['SOX11'])
#0.9439625

print(bkg_fc['HNRNPH1'])
#1.361956
ecdf(bkg_fc)(bkg_fc['HNRNPH1'])
#0.9697876

pdf('plot_fc_den.pdf')
plot(density(bkg_fc),lwd=3,col='grey70',type='h')
abline(v=bkg_fc['SOX11'],col='red',lwd=2)
abline(v=bkg_fc['HNRNPH1'],col='blue',lwd=2)
dev.off()








z=mat['CTNNB1',]
z1=z[which(label==1)]
z0=z[which(label==0)]

pdf('plot_box_CTNNB1.pdf')
boxplot(z1,z0,col=c('red','grey70'),pch='+',lwd=2)
dev.off()
print(wilcox.test(z1,z0))
#p-value = 2.87e-09
print(t.test(z1,z0))
#p-value = 3.571e-16



















