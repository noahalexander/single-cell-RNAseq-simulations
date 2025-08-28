
#library(monocle3)
library(dplyr)
library(rafalib)
library(ComplexHeatmap)
library(RColorBrewer)
library(MASS)
#library(gcookbook)  
library(scales)
library(ggplot2)

################-------------------------------------------------------------------------------simulate data using params from m22 vs. by de testing in ypd 
#for fc, see betas in earlier versions of code 
#thetas=c(1, 5, 10, 20) 
#avg.exp.level.per.cell=c(.0001, .001, .0025, .05)
#ncells.per.exp = c(100, 200, 500, 1000)

thetas=c(1, 5, 10, 20) 
avg.exp.level.per.cell=c(.001, .002, .05, .5)
ncells.per.exp = c(50, 100, 500, 1000)

effect.size=c(0.025, 0.2, 0.4, .6)

list.of.arrays.g.pval = vector(mode = "list", length = 100)
list.of.arrays.g.estimate = vector(mode = "list", length = 100)
list.of.arrays.g.std.err = vector(mode = "list", length = 100)
list.of.arrays.theta = vector(mode = "list", length = 100)

t.1 = Sys.time()

#n simulations
#Each entry (here, n in 1:100) in the lists produced by the first for loop is a similar array of results for 
#one simulated dataset/summary of resulting neg binomial model

for (n in 1:100) { 
  
  #initialize arrays based on sets of parameters  
  array.g.pval <- array(rep(1, length(thetas)*length(avg.exp.level.per.cell)*length(ncells.per.exp)*length(effect.size)), dim=c(length(thetas), length(avg.exp.level.per.cell), length(ncells.per.exp), length(effect.size)))
  #array.g.estimate <- array(rep(1, length(thetas)*length(avg.exp.level.per.cell)*length(ncells.per.exp)*length(effect.size)), dim=c(length(thetas), length(avg.exp.level.per.cell), length(ncells.per.exp), length(effect.size)))
  #array.g.std.err <- array(rep(1, length(thetas)*length(avg.exp.level.per.cell)*length(ncells.per.exp)*length(effect.size)), dim=c(length(thetas), length(avg.exp.level.per.cell), length(ncells.per.exp), length(effect.size)))
  
  for (i in 1:length(thetas)) {
    for (j in 1:length(avg.exp.level.per.cell)) {
      for (k in 1:length(ncells.per.exp)) {
        for (l in 1:length(effect.size)) {
          nUMIs=1000+750*rexp(ncells.per.exp[k]*2,1) 
          # biological covariates of interest are two genotypes and two experimental conditions
          g1=rep(0,ncells.per.exp[k])
          g2=rep(1,ncells.per.exp[k])
          #one way to specify the experiment design matrix
          X=cbind(log(nUMIs),g=c(g1,g2))
          #simulate some effects
          betas=c(1,effect.size[l])
          #get the expected mean of expression  
          pMu=(X %*% betas)[,1]
          #push into count space with neg bin over-dispersion 
          ycounts=rnegbin(length(pMu),exp(pMu)*avg.exp.level.per.cell[j], theta=thetas[i])
          model_i_j_k_l=glm.nb(ycounts~X)
          # model_i_j_k_l=glm(ycounts ~X, family = "quasipoisson")
          psum.model_i_j_k_l = summary(model_i_j_k_l)
          array.g.pval[i,j,k,l] = psum.model_i_j_k_l$coefficients[3,4]

        }
      }
    }
  }
  
  print(n)
  #add each array to the list of arrays or do this in a better way
  list.of.arrays.g.pval[[n]] = array.g.pval
  #list.of.arrays.g.estimate[[n]] = array.g.estimate
  #list.of.arrays.g.std.err[[n]] = array.g.std.err
  #list.of.arrays.theta[[n]] = array.theta
  
}

t.2 = Sys.time()
time.100.iterations = t.2-t.1
time.100.iterations



#list.of.lists = list(list.of.arrays.g.pval, list.of.arrays.g.estimate, list.of.arrays.g.std.err, list.of.arrays.theta)

#consider packakging in an S4 object but for now this works

#saveRDS(list.of.lists, file = "list.of.lists.20221109.RDS")
#saveRDS(list.of.arrays.g.estimate, file="list.of.arrays.g.estimate.20221109.RDS")
#saveRDS(list.of.arrays.g.pval, file="list.of.arrays.g.pval.20221118.RDS")
#saveRDS(list.of.arrays.theta, file="list.of.arrays.theta.20221109.RDS")
#saveRDS(list.of.arrays.g.std.err, file="list.of.arrays.g.std.err.20221109.RDS")

########---------------------------loop through the list of 100 arrays (or however many sims n is) of 256 pvalues (i*J*k*l, here 4*4*4*4, values) to make a list of vectors 

#the lists of arrays are organized such that any element [i,j,k,l] of any of the arrays will correspond to the same set of parameters
#it will be necessary to extract all the values of each specific set of values [i,j,k,l] in each array and then see how many of the 100 simulations using that
#set of parameters has pvalues under a threshold of 0.01. 

iter.list = vector(mode = "list", length = 256)
avg.exp.level.per.cell.list = vector(mode = "list", length = 256)
effect.size.list = vector(mode = "list", length = 256)
theta.list = vector(mode = "list", length = 256)
ncells.per.exp.list = vector(mode = "list", length = 256)

iter.list.params = vector(mode = "list", length = 256)
count = 1
for (i in 1:length(thetas)) {
  for (j in 1:length(avg.exp.level.per.cell)) {
    for (k in 1:length(ncells.per.exp)) {
      for (l in 1:length(effect.size)) {
        iter.list[[count]] = vector()
        for (n in 1:100){
          #theta.list[[count]][n] = list.of.arrays.g.pval[[n]][i,j,k,l]
          #avg.exp.level.per.cell.list[[count]][n] = list.of.arrays.[[n]][i,j,k,l]
          #ncells.per.exp.list[[count]][n] = list.of.arrays.g.pval[[n]][i,j,k,l]
          #effect.size.list[[count]][n] = list.of.arrays.g.pval[[n]][i,j,k,l]
          iter.list[[count]][n] = list.of.arrays.g.pval[[n]][i,j,k,l]
          iter.list.params[[count]][n] = paste0("theta: ",thetas[i],"avg.exp.per.cell: ",avg.exp.level.per.cell[j],"ncells: ",ncells.per.exp[k],"effect.size:", effect.size[l])
        }
        count = count+1
      }
    }
  }
}


rows=vector()
#get vector of rownames for mat
for (i in 1:length(iter.list.params)){
  rows[i] = iter.list.params[[i]][1]
}

#organize matrix with pvals from each simulation 
mat = matrix(nrow = 256, ncol = 100)
rownames(mat) = rows

for (i in 1:length(iter.list)) {
  mat[i,] = iter.list[[i]]
  mat[i,] = p.adjust(mat[i,], method = "BH")
}



#dd = cbind.data.frame(pvals.under.0.01, pvals.under.0.05)
#rownames(dd) = rows
#how many p values are at or below some threshold (I'm using 0.01 and 0.05 for the moment)

pvals.under.0.01=vector()
pvals.under.0.05=vector()
for (i in 1:nrow(mat)) {
  pvals.under.0.01[i] = length(mat[i,][mat[i,]<0.01])
  pvals.under.0.05[i] = length(mat[i,][mat[i,]<0.05])
  
}


da = cbind.data.frame(pvals.under.0.01)
rownames(da) = rows
#pheatmap(mat, annotation_row = da, cluster_rows = F, cluster_cols = F)
rownames(da) = make.unique(substr(rownames(da), 1,5))
rownames(mat) = make.unique(substr(rownames(mat), 1,5))
#pheatmap(mat, annotation_row = da, cluster_rows = F, cluster_cols = F)

#generate vectors of parameter identities for plotting later 
param.theta = vector()
param.avg.exp.level.per.cell = vector()
param.ncells.per.exp = vector()
param.effect.size = vector()

count=1
for (i in 1:length(thetas)) {
  for (j in 1:length(avg.exp.level.per.cell)) {
    for (k in 1:length(ncells.per.exp)) {
      for (l in 1:length(effect.size)) {
        
        print(count)
        param.theta[count] = thetas[i]
        param.avg.exp.level.per.cell[count] = avg.exp.level.per.cell[j]
        param.ncells.per.exp[count] = ncells.per.exp[k]
        param.effect.size[count] = effect.size[l]
        count = count+1
        
        
      }
    }
  }
}

#organize vectors of params and pval counts for 100 sims per set of parameters 
#plot power af a function of effect size facet by various relevant parameters 
df = cbind.data.frame(da$pvals.under.0.01, param.theta , param.avg.exp.level.per.cell, param.ncells.per.exp, param.effect.size)
rownames(df) = rownames(da)



p <- ggplot(data = df, aes(x = param.effect.size, y = (df$`da$pvals.under.0.01`)/100, colour=param.theta, group=param.theta))

p  + facet_wrap(param.avg.exp.level.per.cell~param.ncells.per.exp) + 
    geom_point() +
    scale_colour_gradient2(low = muted("red"), mid = ("white"), high = muted("blue"), midpoint = 10) + 
    labs(x="effect size", y="power" ) + 
    geom_line() +
     theme( axis.title.x = element_text(size=20), axis.title.y.left =element_text(size=20) )




#########-----------------additional plotting 
p <- ggplot(data = df, aes(x = param.effect.size, y = (df$`da$pvals.under.0.01`)/100, colour=param.theta)) +geom_point()
p = p  + facet_wrap(~param.ncells.per.exp) + scale_colour_gradient2(low = muted("red"), mid = ("white"), high = muted("blue"), midpoint = 10) + geom_jitter(position = position_jitter(width = 0.05, height = 0))
p + labs(x="effect size", y="power")



#############

