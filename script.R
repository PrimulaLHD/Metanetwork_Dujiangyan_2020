###Data analysis for metanetwork in Dujiangyan,R2 according to reviewers' comments
###by Hai-Dong Li



##load packages and functions
rm(list=ls()) # befor start remove everything include functions in R
setwd("C:/Dataset/Metanetwork/Finnal") # set working directory
# All packages
library("devtools")
library("openxlsx")
library("bipartite")
#library("plyr")
#library("dplyr")
library("betalink")
library("ggplot2")
library("mgcv")
library("vegan")
library("igraph")
library("networkD3")
library("car")
library("betapart")
library("patchwork")

# Functions
sumbygroup <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, .drop=TRUE) {
    require(plyr)
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    # This does the summary. For each group's data frame, return a vector with
    # sum
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(sum = sum(xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the " sum" column    
    datac <- rename(datac, c("sum" = measurevar))
    
    return(datac)
}

right.matrix <- function(x) {
    m<-as.matrix(x[,-1])
    rownames(m)<-x[,1]
    m
}

m2f <- function(m){# make interaction matrix to dataframe
    df1 <- as.data.frame.table(m)
    n <- length(df1[,1])
    x=df1[,1]
    y=df1[,2]
    IU <- vector()
    for (i in 1:n) {
        IU[i] <- paste(x[i],y[i],sep = "-")
    }
    df <- data.frame(cbind(IU,df1$Freq))
    names(df) <- c("Interaction_Unit","Value")
    df$Value <-as.numeric(levels(df$Value))[df$Value]
    df
}





##Load dataset and prepare 

# Load data set of 2014
Net_DF_2014on <- read.xlsx("DispersalNet2014.xlsx",sheet = 1) # read data on the tree 2014
names(Net_DF_2014on) # Check data

Net_DF_2014under <- read.xlsx("DispersalNet2014.xlsx",sheet = 2)# read data under the tree 2014
names(Net_DF_2014under) # Check data

# exclude uneaten birds and mammals from network
Net_DF_2014on_Eat <- droplevels(Net_DF_2014on[which(Net_DF_2014on$Eaten=="Yes"),])
networkdata_2014_on <- Net_DF_2014on_Eat[,c(5,7,10,11)]
head(networkdata_2014_on)

Net_DF_2014under_Eat <- droplevels(Net_DF_2014under[which(Net_DF_2014under$Eaten=="Yes"),])
Net_DF_2014under_Bird <- droplevels(Net_DF_2014under_Eat[which(Net_DF_2014under_Eat$Bird=="Yes"),])
networkdata_2014_under <- Net_DF_2014under_Bird[,c(2,4,13,14)]
head(networkdata_2014_under)

names(networkdata_2014_on) <- c("Patch","TreeSpecies","AnimalSpecies","Animal_Number")
sort(unique(Net_DF_2014on$TreeSpecies))==sort(unique(Net_DF_2014under$TreeSpecies)) #should be true

df_network2014 <- rbind(networkdata_2014_on,networkdata_2014_under)

#convert to bipartite matrix , this is list object of 2014
frame2webs(df_network2014,varnames=c("TreeSpecies","AnimalSpecies","Patch",
                                     "Animal_Number"))->Networkdata2014
Patch <- names(Networkdata2014)

#bipartite network for all patch combined, full network of 2014
networkdata2014_2 <-df_network2014[,2:3]
networkdata2014_com<-table(networkdata2014_2)

#*****************************************************************************************#

#Load dataset of 2015

Net_DF_2015on <- read.xlsx("DispersalNet2015.xlsx",sheet = 1) # read data on the tree 2015
names(Net_DF_2015on) # Check data
unique(Net_DF_2015on$TreeSpecies)

Net_DF_2015under <- read.xlsx("DispersalNet2015.xlsx",sheet = 2)# read data under the tree 2015
names(Net_DF_2015under) # Check data
unique(Net_DF_2015under$TreeSpecies)

sort(unique(Net_DF_2015on$TreeSpecies))==sort(unique(Net_DF_2015under$TreeSpecies)) #should be true

# exclude uneaten birds and mammals from network of 2015
Net_DF_2015on_Eat <- droplevels(Net_DF_2015on[which(Net_DF_2015on$Eaten=="Yes"),])
networkdata_2015_on <- Net_DF_2015on_Eat[,c(5,7,10,11)]
head(networkdata_2015_on)


Net_DF_2015under_Eat <- droplevels(Net_DF_2015under[which(Net_DF_2015under$Eaten=="Yes"),])
Net_DF_2015under_Bird <- droplevels(Net_DF_2015under_Eat[which(Net_DF_2015under_Eat$Bird=="Yes"),])
networkdata_2015_under <- Net_DF_2015under_Bird[,c(2,4,13,14)]
head(networkdata_2015_under)

#names(networkdata_2015_on) <- c("Patch","TreeSpecies","AnimalSpecies","Animal_Number")
df_network2015 <- rbind(networkdata_2015_on,networkdata_2015_under)
#convert to bipartite matrix, this is a list of 2015
frame2webs(df_network2015,
           varnames=c("TreeSpecies","AnimalSpecies","Patch","Animal_Number"))->Networkdata2015
Patch <- names(Networkdata2015)

#bipartite network for all patch combined of 2015
networkdata2015_2 <-df_network2015[,2:3]
networkdata2015_com<-table(networkdata2015_2)


df_network_two <- rbind(df_network2014,df_network2015)

save(df_network_two,file="df_network_two.raw.Rdata")

frame2webs(df_network_two,varnames=c("TreeSpecies","AnimalSpecies","Patch",
                                     "Animal_Number"))->Networks

save(df_network_two,Networks,file="df_network_two0213.raw.Rdata")

# on-the-tree network
on_the_tree <- rbind(networkdata_2014_on,networkdata_2015_on)
on_the_tree$webID <- "on"
web_on <- frame2webs(on_the_tree,varnames=c("TreeSpecies","AnimalSpecies",
                                               "webID", "Animal_Number"))[[1]]

sum(web_on>0)

# under-the-tree network
under_the_tree <- rbind(networkdata_2014_under,networkdata_2015_under)
under_the_tree$webID <- "under"
web_under <- frame2webs(under_the_tree,varnames=c("TreeSpecies","AnimalSpecies",
                                                  "webID", "Animal_Number"))[[1]]
sum(web_under>0)

# full network
df_network_two$webID <- "all"
web_all <- frame2webs(df_network_two,varnames=c("TreeSpecies","AnimalSpecies",
                                                "webID", "Animal_Number"))[[1]]
sum(web_all>0)


## Build quatitative metanetwork and quanlitative metanetwork

plain_df <- function(df,patch="patch",plant="Plant",
                     animal="Animal",Value="feq"){# make dataframe to individual interactions
    
    require(dplyr)
    require(tidyr)
    df1 <- subset(df,df[,Value]>0)
    n <- length(df1[,1])
    x=df1[,plant]
    y=df1[,animal]
    z=df1[,patch]
    IU <- vector()
    for (i in 1:n) {
        IU[i] <- paste(x[i],y[i],z[i],sep = "-")
    }
    df2 <- data.frame(cbind(IU,df1[,c(Value)]))
    names(df2) <- c("Interaction_Unit","Value")
    df2$Value <-as.numeric(levels(df2$Value))[df2$Value]
    
    df3 <-sumbygroup(df2,measurevar = "Value",
                     groupvars = c("Interaction_Unit"))
    
    
    plain_df1 <- list()
    m=length(df3[,1])
    for (i in 1:m) {
        plain_df1[[i]] <- data.frame(Interaction_Unit=rep(df3[i,1],df3[i,2]))
    }
    plain_df2 <-do.call("rbind",plain_df1)
    plain_df2 %>% 
        separate(Interaction_Unit, c("plant","animal","patch"), "-") ->plain_df
    return(plain_df)
}


metanetwork_df <- plain_df(df=df_network_two,patch = "Patch",
               plant = "TreeSpecies",animal = "AnimalSpecies",Value = "Animal_Number")

#write.csv(metanetwork_df,"metanetwork_dataframe.csv")

df2metanetwork <- function(df,patch="patch",plant="Plant",
                           animal="Animal",Value="feq",quantitative=TRUE){
    df1 <- subset(df,df[,Value]>0)
    n <- length(df1[,1])
    x=df1[,plant]
    y=df1[,animal]
    IU <- vector()
    for (i in 1:n) {
        IU[i] <- paste(x[i],y[i],sep = "-")
    }
    df2 <- data.frame(cbind(IU,df1[,c(patch,Value)]))
    names(df2) <- c("Interaction_Unit","patch","Value")
    df3 <-sumbygroup(df2,measurevar = "Value",
                     groupvars = c("patch","Interaction_Unit"))
    xtabs(Value~Interaction_Unit + patch,data = df3)-> web
    attr(web,"class") <- NULL
    attr(web,"call") <- NULL
    if (quantitative=="TRUE"){
        metanetwork <- web
    }else{
        web[web>0]=1
        metanetwork <- web
    }
    return(metanetwork)
}



head(df_network_two)

q_metanetwork <-df2metanetwork(df_network_two,patch = "Patch",
               plant = "TreeSpecies",animal = "AnimalSpecies",Value = "Animal_Number")

b_metanetwork <-df2metanetwork(df_network_two,patch = "Patch",
                               plant = "TreeSpecies",animal = "AnimalSpecies",
                               Value = "Animal_Number",quantitative = FALSE)


metanetwork_df <- read.csv("metanetwork_dataframe.csv",header = T)
head(metanetwork_df)


metanet_df_drop_on <-droplevels(metanetwork_df[which(metanetwork_df[,"guild"]!="on"),])
metanet_df_drop_under <-droplevels(metanetwork_df[which(metanetwork_df[,"guild"]!="under"),])
metanet_df_drop_shared <-droplevels(metanetwork_df[which(metanetwork_df[,"guild"]!="shared"),])



q_metanetwork_drop_on <- df2metanetwork(df=metanet_df_drop_on,patch="patch",
                               plant="plant",animal="animal",
                               Value="feq",quantitative=TRUE)
q_metanetwork_drop_under <- df2metanetwork(df=metanet_df_drop_under,patch="patch",
                                  plant="plant",animal="animal",
                                  Value="feq",quantitative=TRUE)
q_metanetwork_drop_shared <- df2metanetwork(df=metanet_df_drop_shared,patch="patch",
                                   plant="plant",animal="animal",
                                   Value="feq",quantitative=TRUE)


b_metanetwork_drop_on <- df2metanetwork(df=metanet_df_drop_on,patch="patch",
                                        plant="plant",animal="animal",
                                        Value="feq",quantitative=FALSE)
b_metanetwork_drop_under <- df2metanetwork(df=metanet_df_drop_under,patch="patch",
                                           plant="plant",animal="animal",
                                           Value="feq",quantitative=FALSE)
b_metanetwork_drop_shared <- df2metanetwork(df=metanet_df_drop_shared,patch="patch",
                                            plant="plant",animal="animal",
                                            Value="feq",quantitative=FALSE)



## 1. metanetwork structure
### network visualization

plotweb(q_metanetwork) # bipartite network
plotweb(b_metanetwork) # bipartite network


metanetwork3d_q <- graph.incidence(q_metanetwork,weighted = T)#######Use 3D network

wc <- cluster_walktrap(metanetwork3d_q)
members <- membership(wc)
# Convert to object suitable for networkD3
metanetwork_d3 <- igraph_to_networkD3(metanetwork3d_q, group = members)

# Create force directed network plot
forceNetwork(Links = metanetwork_d3$links, Nodes = metanetwork_d3$nodes, 
             Source = 'source', Target = 'target', Value = "value",
             NodeID = 'name', Group = 'group', opacity = 1, bounded = T)


#------------------------------------------------------------------------------
metanetwork3d_q <- graph.incidence(q_metanetwork_drop_on,weighted = T)#######Use 3D network

wc <- cluster_walktrap(metanetwork3d_q)
members <- membership(wc)
# Convert to object suitable for networkD3
metanetwork_d3 <- igraph_to_networkD3(metanetwork3d_q, group = members)

# Create force directed network plot
forceNetwork(Links = metanetwork_d3$links, Nodes = metanetwork_d3$nodes, 
             Source = 'source', Target = 'target', Value = "value",
             NodeID = 'name', Group = 'group', opacity = 1, bounded = T)

#------------------------------------------------------------------------------
metanetwork3d_q <- graph.incidence(q_metanetwork_drop_under,weighted = T)#######Use 3D network

wc <- cluster_walktrap(metanetwork3d_q)
members <- membership(wc)
# Convert to object suitable for networkD3
metanetwork_d3 <- igraph_to_networkD3(metanetwork3d_q, group = members)

# Create force directed network plot
forceNetwork(Links = metanetwork_d3$links, Nodes = metanetwork_d3$nodes, 
             Source = 'source', Target = 'target', Value = "value",
             NodeID = 'name', Group = 'group', opacity = 1, bounded = T)

#------------------------------------------------------------------------------
metanetwork3d_q <- graph.incidence(q_metanetwork_drop_shared,weighted = T)#######Use 3D network

wc <- cluster_walktrap(metanetwork3d_q)
members <- membership(wc)
# Convert to object suitable for networkD3
metanetwork_d3 <- igraph_to_networkD3(metanetwork3d_q, group = members)

# Create force directed network plot
forceNetwork(Links = metanetwork_d3$links, Nodes = metanetwork_d3$nodes, 
             Source = 'source', Target = 'target', Value = "value",
             NodeID = 'name', Group = 'group', opacity = 1, bounded = T)





##1.1 quantitative analysis
myFunwithNull <- function(network,nulls,index){
    obs <- networklevel(network,index=index)
    null.value <-sapply(nulls,function(x)networklevel(x,index=index))
    mean.null <- mean(null.value)
    Z<-(obs-mean.null)/sd(null.value)
    return(cbind(obs,mean.null,Z))
}

MFunwithNull <- function(network,nulls){
    obs <- computeModules(network)@likelihood
    null.value <-sapply(nulls,function(x)computeModules(x)@likelihood)
    mean.null <- mean(null.value)
    Z<-(obs-mean.null)/sd(null.value)
    return(cbind(obs,mean.null,Z))
}

H2FunwithNull <- function(network,nulls){
    obs <- H2fun(network)[1]
    null.value <-sapply(nulls,function(x)H2fun(x)[1])
    mean.null <- mean(null.value)
    Z<-(obs-mean.null)/sd(null.value)
    return(cbind(obs,mean.null,Z))
}


set.seed(123)
nulls_q <- nullmodel(q_metanetwork, N=1000,method = "r2d")

(con_q <- networklevel(q_metanetwork, index="connectance")) # connectance
#0.1492537
(con_nulls_q <- myFunwithNull(q_metanetwork,nulls_q,index = "connectance"))

#                 obs    mean.null    Z
#connectance 0.1492537 0.3797865 -45.08766

M_q <- computeModules(t(q_metanetwork), method="Beckett") # Modularity
(M_q@likelihood)
# M = 0.45
plotModuleWeb(M_q)

(M_nulls_q <- MFunwithNull(q_metanetwork,nulls = nulls_q))
#         obs      mean.null  Z
#[1,] 0.4514764 0.0767809 105.6306



(nestedness_nulls_q <- myFunwithNull(q_metanetwork,
                                    nulls = nulls_q,index = "weighted NODF")) 

#obs mean.null        Z
#weighted NODF 10.85949    41.049 -23.0028



##drop on
set.seed(123)
nulls_q_drop_on <- nullmodel(q_metanetwork_drop_on, N=1000,method = "r2d")


(con_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,nulls_q_drop_on,
                                      index = "connectance"))

#              obs       mean.null   Z
#connectance 0.1407343 0.3546058 -39.35942


(M_nulls_q_drop_on <- MFunwithNull(q_metanetwork_drop_on,nulls = nulls_q_drop_on))

#obs mean.null       Z
#[1,] 0.6005207 0.1235814 80.1134




(nestedness_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,
                                     nulls = nulls_q_drop_on,index = "weighted NODF")) 
#                obs      mean.null   Z
#weighted NODF 10.0844  33.80326 -18.27791



##drop under
set.seed(123)
nulls_q_drop_under <- nullmodel(q_metanetwork_drop_under, N=1000,method = "r2d")


(con_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,nulls_q_drop_under,
                                      index = "connectance"))
#                obs      mean.null    Z
#connectance 0.1591973 0.3997706 -42.85222

(M_nulls_q_drop_under <- MFunwithNull(q_metanetwork_drop_under,nulls = nulls_q_drop_under))

#        obs   mean.null        Z
#[1,] 0.4469315 0.07285514 100.6145


(nestedness_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,
                                           nulls = nulls_q_drop_under,index = "weighted NODF")) 

#obs mean.null         Z
#weighted NODF 12.27524  44.19098 -21.77546


##drop shared
set.seed(123)
nulls_q_drop_shared <- nullmodel(q_metanetwork_drop_shared, N=1000,method = "r2d")


(con_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,nulls_q_drop_shared,
                                         index = "connectance"))
#              obs      mean.null      Z
#connectance 0.1431953 0.3688225 -30.87714

(M_nulls_q_drop_shared <- MFunwithNull(q_metanetwork_drop_shared,nulls = nulls_q_drop_shared))

#         obs     mean.null    Z
#[1,] 0.354539 0.06231726 65.51993

(nestedness_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,
                                              nulls = nulls_q_drop_shared,
                                              index = "weighted NODF")) 

#                obs     mean.null     Z
#weighted NODF 9.964156  45.81443 -16.5851



net_structure_q <- rbind(con_nulls_q,M_nulls_q,nestedness_nulls_q,
      con_nulls_q_drop_on,M_nulls_q_drop_on,nestedness_nulls_q_drop_on,
      con_nulls_q_drop_under,M_nulls_q_drop_under,nestedness_nulls_q_drop_under,
      con_nulls_q_drop_shared,M_nulls_q_drop_shared,nestedness_nulls_q_drop_shared)

net_structure_q
#                 obs    mean.null    Z
#connectance 0.1492537 0.3797865 -45.08766
#[1,] 0.4514764 0.0767809 105.6306
#weighted NODF 10.85949    41.049 -23.0028
#connectance 0.1407343 0.3546058 -39.35942
#[1,] 0.6005207 0.1235814 80.1134
#weighted NODF 10.0844  33.80326 -18.27791
#connectance 0.1591973 0.3997706 -42.85222
#[1,] 0.4469315 0.07285514 100.6145
#weighted NODF 12.27524  44.19098 -21.77546
#connectance 0.1431953 0.3688225 -30.87714
#[1,] 0.354539 0.06231726 65.51993
#weighted NODF 9.964156  45.81443 -16.5851


##redundancy
Sq <- function(web,level="lower"){
    if(level=="lower"){
        I <- nrow(web)
        A <- rowSums(web)
    }else{
        I <- ncol(web)
        A <- colSums(web)
    }
    m = sum(web)
    EH <- bipartite::specieslevel(web = web,level = level,
                                  index = "effective partners")
    Sq <-0
    for (i in 1:I) {
        Sq <- Sq + (A[i]/m)*EH[i,1]
    }
    names(Sq) <-NULL
    return(Sq)
}


Sq(q_metanetwork)#3.848041
Sq(q_metanetwork_drop_on) #3.1069
Sq(q_metanetwork_drop_under) #3.937121
Sq(q_metanetwork_drop_shared) #4.12475





## 1.2 Binary analysis
myFunwithNull_vegan <- function(network,nulls,index){
    obs <- networklevel(network,index=index)
    null.value <-apply(nulls,3,function(x)networklevel(x,index=index))
    mean.null <- mean(null.value)
    Z<-(obs-mean.null)/sd(null.value)
    return(cbind(obs,mean.null,Z))
}

MFunwithNull_vegan <- function(network,nulls){
    obs <- computeModules(network)@likelihood
    null.value <-apply(nulls,3,function(x)computeModules(x)@likelihood)
    mean.null <- mean(null.value)
    Z<-(obs-mean.null)/sd(null.value)
    return(cbind(obs,mean.null,Z))
}


set.seed(123)
#suggested by Carsten Dormann using fixed-fixed null model

nulls_b  <- (simulate(vegan::nullmodel(b_metanetwork, method="quasiswap"), nsim = 1000))

(con_nulls_b <- myFunwithNull_vegan(b_metanetwork,nulls_b,index = "connectance"))
#obs mean.null   Z
#connectance 0.1492537 0.1492537 NaN

M_b <- computeModules(t(b_metanetwork), method="Beckett") # Modularity
(M_b@likelihood)

plotModuleWeb(M_b)

(M_nulls_b <- MFunwithNull_vegan(b_metanetwork,nulls = nulls_b))

#       obs      mean.null      Z
#[1,] 0.461642 0.4537553 1.768266

(nestedness_nulls_b <- myFunwithNull_vegan(b_metanetwork,
                                     nulls = nulls_b,index = "NODF")) 
#        obs      mean.null   Z
#NODF 19.16523  19.68255 -2.495415



##drop on
set.seed(123)

nulls_b_drop_on <- (simulate(vegan::nullmodel(b_metanetwork_drop_on, method="quasiswap"),
                             nsim = 1000))

(con_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,nulls_b_drop_on,
                                            index = "connectance"))
#                  obs mean.null   Z
#connectance 0.1407343 0.1407343 NaN

(M_nulls_b_drop_on <- MFunwithNull_vegan(b_metanetwork_drop_on,nulls = nulls_b_drop_on))

#          obs    mean.null    Z
#[1,] 0.5021797  0.487332 2.233659

(nestedness_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,
                                           nulls = nulls_b_drop_on,index = "NODF")) 

#          obs    mean.null   Z
#NODF 16.79688  17.58016 -2.522972





##drop under
set.seed(123)
nulls_b_drop_under <- (simulate(vegan::nullmodel(b_metanetwork_drop_under, method="quasiswap"),
                             nsim = 1000))

(con_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,nulls_b_drop_under,
                                            index = "connectance"))
#                  obs mean.null   Z
#connectance 0.1591973 0.1591973 NaN

(M_nulls_b_drop_under <- MFunwithNull_vegan(b_metanetwork_drop_under,nulls = nulls_b_drop_under))

#          obs mean.null        Z
#[1,] 0.444619 0.4263839 3.784796

(nestedness_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,
                                                   nulls = nulls_b_drop_under,index = "NODF")) 

#        obs mean.null         Z
#NODF 21.823  22.51838 -2.937353




##drop shared
set.seed(123)
nulls_b_drop_shared <- (simulate(vegan::nullmodel(b_metanetwork_drop_shared, method="quasiswap"),
                                nsim = 1000))

(con_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,nulls_b_drop_shared,
                                               index = "connectance"))

#                  obs mean.null   Z
#connectance 0.1431953 0.1431953 NaN


(M_nulls_b_drop_shared <- MFunwithNull_vegan(b_metanetwork_drop_shared,nulls = nulls_b_drop_shared))

#         obs   mean.null  Z
#[1,] 0.4839833 0.4837922 0.023648

(nestedness_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,
                                                      nulls = nulls_b_drop_shared,index = "NODF")) 

#obs mean.null         Z
#NODF 18.19669  18.45566 -0.651131


net_structure_b <- rbind(con_nulls_b,M_nulls_b,nestedness_nulls_b,
                         con_nulls_b_drop_on,M_nulls_b_drop_on,nestedness_nulls_b_drop_on,
                         con_nulls_b_drop_under,M_nulls_b_drop_under,nestedness_nulls_b_drop_under,
                         con_nulls_b_drop_shared,M_nulls_b_drop_shared,nestedness_nulls_b_drop_shared)

net_structure_b

#obs mean.null   Z
#connectance 0.1492537 0.1492537 NaN
#[1,] 0.461642 0.4537553 1.768266
#NODF 19.16523  19.68255 -2.495415
#connectance 0.1407343 0.1407343 NaN
#[1,] 0.5021797  0.487332 2.233659
#NODF 16.79688  17.58016 -2.522972
#connectance 0.1591973 0.1591973 NaN
#[1,] 0.444619 0.4263839 3.784796
#NODF 21.823  22.51838 -2.937353
#connectance 0.1431953 0.1431953 NaN
#[1,] 0.4839833 0.4837922 0.023648
#NODF 18.19669  18.45566 -0.651131
#########----------------------------------------------------------------------------
##Other null models
set.seed(123)
nulls_q <- nullmodel(q_metanetwork, N=1000,method = "swap.web")

(con_q <- networklevel(q_metanetwork, index="connectance")) # connectance
#0.1492537
(con_nulls_q <- myFunwithNull(q_metanetwork,nulls_q,index = "connectance"))


M_q <- computeModules(t(q_metanetwork), method="Beckett") # Modularity
(M_q@likelihood)
# M = 0.45
plotModuleWeb(M_q)

(M_nulls_q <- MFunwithNull(q_metanetwork,nulls = nulls_q))


(nestedness_nulls_q <- myFunwithNull(q_metanetwork,
                                     nulls = nulls_q,index = "weighted NODF")) 


##drop on
set.seed(123)
nulls_q_drop_on <- nullmodel(q_metanetwork_drop_on, N=1000,method = "swap.web")


(con_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,nulls_q_drop_on,
                                      index = "connectance"))



(M_nulls_q_drop_on <- MFunwithNull(q_metanetwork_drop_on,nulls = nulls_q_drop_on))


(nestedness_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,
                                             nulls = nulls_q_drop_on,index = "weighted NODF")) 


##drop under
set.seed(123)
nulls_q_drop_under <- nullmodel(q_metanetwork_drop_under, N=1000,method = "swap.web")


(con_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,nulls_q_drop_under,
                                         index = "connectance"))

(M_nulls_q_drop_under <- MFunwithNull(q_metanetwork_drop_under,nulls = nulls_q_drop_under))


(nestedness_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,
                                                nulls = nulls_q_drop_under,index = "weighted NODF")) 


##drop shared
set.seed(123)
nulls_q_drop_shared <- nullmodel(q_metanetwork_drop_shared, N=1000,method = "swap.web")


(con_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,nulls_q_drop_shared,
                                          index = "connectance"))

(M_nulls_q_drop_shared <- MFunwithNull(q_metanetwork_drop_shared,nulls = nulls_q_drop_shared))


(nestedness_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,
                                                 nulls = nulls_q_drop_shared,
                                                 index = "weighted NODF")) 



net_structure_q_swap.web <- rbind(con_nulls_q,M_nulls_q,nestedness_nulls_q,
                                  con_nulls_q_drop_on,M_nulls_q_drop_on,nestedness_nulls_q_drop_on,
                                  con_nulls_q_drop_under,M_nulls_q_drop_under,nestedness_nulls_q_drop_under,
                                  con_nulls_q_drop_shared,M_nulls_q_drop_shared,nestedness_nulls_q_drop_shared)

net_structure_q_swap.web
write.csv(net_structure_q_swap.web,"./R2/output/net_structure_q_swap.web.csv")

#################

set.seed(123)
nulls_q <- nullmodel(q_metanetwork, N=1000,method = "vaznull")

(con_q <- networklevel(q_metanetwork, index="connectance")) # connectance
#0.1492537
(con_nulls_q <- myFunwithNull(q_metanetwork,nulls_q,index = "connectance"))


M_q <- computeModules(t(q_metanetwork), method="Beckett") # Modularity
(M_q@likelihood)
# M = 0.45
plotModuleWeb(M_q)

(M_nulls_q <- MFunwithNull(q_metanetwork,nulls = nulls_q))


(nestedness_nulls_q <- myFunwithNull(q_metanetwork,
                                     nulls = nulls_q,index = "weighted NODF")) 


##drop on
set.seed(123)
nulls_q_drop_on <- nullmodel(q_metanetwork_drop_on, N=1000,method = "vaznull")


(con_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,nulls_q_drop_on,
                                      index = "connectance"))


(M_nulls_q_drop_on <- MFunwithNull(q_metanetwork_drop_on,nulls = nulls_q_drop_on))



(nestedness_nulls_q_drop_on <- myFunwithNull(q_metanetwork_drop_on,
                                             nulls = nulls_q_drop_on,index = "weighted NODF")) 



##drop under
set.seed(123)
nulls_q_drop_under <- nullmodel(q_metanetwork_drop_under, N=1000,method = "vaznull")


(con_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,nulls_q_drop_under,
                                         index = "connectance"))

(M_nulls_q_drop_under <- MFunwithNull(q_metanetwork_drop_under,nulls = nulls_q_drop_under))



(nestedness_nulls_q_drop_under <- myFunwithNull(q_metanetwork_drop_under,
                                                nulls = nulls_q_drop_under,index = "weighted NODF")) 


##drop shared
set.seed(123)
nulls_q_drop_shared <- nullmodel(q_metanetwork_drop_shared, N=1000,method = "vaznull")


(con_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,nulls_q_drop_shared,
                                          index = "connectance"))

(M_nulls_q_drop_shared <- MFunwithNull(q_metanetwork_drop_shared,nulls = nulls_q_drop_shared))


(nestedness_nulls_q_drop_shared <- myFunwithNull(q_metanetwork_drop_shared,
                                                 nulls = nulls_q_drop_shared,
                                                 index = "weighted NODF")) 



net_structure_q_vaznull <- rbind(con_nulls_q,M_nulls_q,nestedness_nulls_q,
                                 con_nulls_q_drop_on,M_nulls_q_drop_on,nestedness_nulls_q_drop_on,
                                 con_nulls_q_drop_under,M_nulls_q_drop_under,nestedness_nulls_q_drop_under,
                                 con_nulls_q_drop_shared,M_nulls_q_drop_shared,nestedness_nulls_q_drop_shared)

net_structure_q_vaznull
write.csv(net_structure_q_vaznull,"./R2/output/net_structure_q_vaznull.csv")

###*******************************
b_metanetwork <-t(b_metanetwork)
b_metanetwork_drop_on <-t(b_metanetwork_drop_on)
b_metanetwork_drop_under <-t(b_metanetwork_drop_under)
b_metanetwork_drop_shared <-t(b_metanetwork_drop_shared)

set.seed(123)
nulls_b  <- (simulate(vegan::nullmodel(b_metanetwork, method="r1"), nsim = 1000))

(con_nulls_b <- myFunwithNull_vegan(b_metanetwork,nulls_b,index = "connectance"))

M_b <- computeModules(t(b_metanetwork), method="Beckett") # Modularity
(M_b@likelihood)

plotModuleWeb(M_b)

(M_nulls_b <- MFunwithNull_vegan(b_metanetwork,nulls = nulls_b))


(nestedness_nulls_b <- myFunwithNull_vegan(b_metanetwork,
                                           nulls = nulls_b,index = "NODF")) 

##drop on
set.seed(123)

nulls_b_drop_on <- (simulate(vegan::nullmodel(b_metanetwork_drop_on, method="r1"),
                             nsim = 1000))

(con_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,nulls_b_drop_on,
                                            index = "connectance"))

(M_nulls_b_drop_on <- MFunwithNull_vegan(b_metanetwork_drop_on,nulls = nulls_b_drop_on))


(nestedness_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,
                                                   nulls = nulls_b_drop_on,index = "NODF")) 


##drop under
set.seed(123)
nulls_b_drop_under <- (simulate(vegan::nullmodel(b_metanetwork_drop_under, method="r1"),
                                nsim = 1000))

(con_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,nulls_b_drop_under,
                                               index = "connectance"))

(M_nulls_b_drop_under <- MFunwithNull_vegan(b_metanetwork_drop_under,nulls = nulls_b_drop_under))


(nestedness_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,
                                                      nulls = nulls_b_drop_under,index = "NODF")) 


##drop shared
set.seed(123)
nulls_b_drop_shared <- (simulate(vegan::nullmodel(b_metanetwork_drop_shared, method="r1"),
                                 nsim = 1000))

(con_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,nulls_b_drop_shared,
                                                index = "connectance"))



(M_nulls_b_drop_shared <- MFunwithNull_vegan(b_metanetwork_drop_shared,nulls = nulls_b_drop_shared))


(nestedness_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,
                                                       nulls = nulls_b_drop_shared,index = "NODF")) 



net_structure_b_r1 <- rbind(con_nulls_b,M_nulls_b,nestedness_nulls_b,
                            con_nulls_b_drop_on,M_nulls_b_drop_on,nestedness_nulls_b_drop_on,
                            con_nulls_b_drop_under,M_nulls_b_drop_under,nestedness_nulls_b_drop_under,
                            con_nulls_b_drop_shared,M_nulls_b_drop_shared,nestedness_nulls_b_drop_shared)

net_structure_b_r1

write.csv(net_structure_b_r1,"./R2/output/net_structure_b_r1-2.csv")


#####
set.seed(123)
nulls_b  <- (simulate(vegan::nullmodel(b_metanetwork, method="r2"), nsim = 1000))

(con_nulls_b <- myFunwithNull_vegan(b_metanetwork,nulls_b,index = "connectance"))

M_b <- computeModules(t(b_metanetwork), method="Beckett") # Modularity
(M_b@likelihood)

plotModuleWeb(M_b)

(M_nulls_b <- MFunwithNull_vegan(b_metanetwork,nulls = nulls_b))


(nestedness_nulls_b <- myFunwithNull_vegan(b_metanetwork,
                                           nulls = nulls_b,index = "NODF")) 

##drop on
set.seed(123)

nulls_b_drop_on <- (simulate(vegan::nullmodel(b_metanetwork_drop_on, method="r2"),
                             nsim = 1000))

(con_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,nulls_b_drop_on,
                                            index = "connectance"))

(M_nulls_b_drop_on <- MFunwithNull_vegan(b_metanetwork_drop_on,nulls = nulls_b_drop_on))


(nestedness_nulls_b_drop_on <- myFunwithNull_vegan(b_metanetwork_drop_on,
                                                   nulls = nulls_b_drop_on,index = "NODF")) 


##drop under
set.seed(123)
nulls_b_drop_under <- (simulate(vegan::nullmodel(b_metanetwork_drop_under, method="r2"),
                                nsim = 1000))

(con_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,nulls_b_drop_under,
                                               index = "connectance"))

(M_nulls_b_drop_under <- MFunwithNull_vegan(b_metanetwork_drop_under,nulls = nulls_b_drop_under))


(nestedness_nulls_b_drop_under <- myFunwithNull_vegan(b_metanetwork_drop_under,
                                                      nulls = nulls_b_drop_under,index = "NODF")) 


##drop shared
set.seed(123)
nulls_b_drop_shared <- (simulate(vegan::nullmodel(b_metanetwork_drop_shared, method="r2"),
                                 nsim = 1000))

(con_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,nulls_b_drop_shared,
                                                index = "connectance"))



(M_nulls_b_drop_shared <- MFunwithNull_vegan(b_metanetwork_drop_shared,nulls = nulls_b_drop_shared))


(nestedness_nulls_b_drop_shared <- myFunwithNull_vegan(b_metanetwork_drop_shared,
                                                       nulls = nulls_b_drop_shared,index = "NODF")) 




net_structure_b_r2 <- rbind(con_nulls_b,M_nulls_b,nestedness_nulls_b,
                            con_nulls_b_drop_on,M_nulls_b_drop_on,nestedness_nulls_b_drop_on,
                            con_nulls_b_drop_under,M_nulls_b_drop_under,nestedness_nulls_b_drop_under,
                            con_nulls_b_drop_shared,M_nulls_b_drop_shared,nestedness_nulls_b_drop_shared)

net_structure_b_r2

write.csv(net_structure_b_r2,"./R2/output/net_structure_b_r2-2.csv")

######
#####
set.seed(123)
nulls_b  <- bipartite::nullmodel(b_metanetwork, method=5,N=1000)

(con_nulls_b <- myFunwithNull(b_metanetwork,nulls_b,index = "connectance"))

M_b <- computeModules(t(b_metanetwork), method="Beckett") # Modularity
(M_b@likelihood)

plotModuleWeb(M_b)

(M_nulls_b <- MFunwithNull(b_metanetwork,nulls = nulls_b))


(nestedness_nulls_b <- myFunwithNull(b_metanetwork,
                                           nulls = nulls_b,index = "NODF")) 

##drop on
set.seed(123)

nulls_b_drop_on <- bipartite::nullmodel(b_metanetwork_drop_on, method=5,
                             N = 1000)

(con_nulls_b_drop_on <- myFunwithNull(b_metanetwork_drop_on,nulls_b_drop_on,
                                            index = "connectance"))

(M_nulls_b_drop_on <- MFunwithNull(b_metanetwork_drop_on,nulls = nulls_b_drop_on))


(nestedness_nulls_b_drop_on <- myFunwithNull(b_metanetwork_drop_on,
                                                   nulls = nulls_b_drop_on,index = "NODF")) 


##drop under
set.seed(123)
nulls_b_drop_under <- bipartite::nullmodel(b_metanetwork_drop_under, method=5,
                                N = 1000)

(con_nulls_b_drop_under <- myFunwithNull(b_metanetwork_drop_under,nulls_b_drop_under,
                                               index = "connectance"))

(M_nulls_b_drop_under <- MFunwithNull(b_metanetwork_drop_under,nulls = nulls_b_drop_under))


(nestedness_nulls_b_drop_under <- myFunwithNull(b_metanetwork_drop_under,
                                                      nulls = nulls_b_drop_under,index = "NODF")) 


##drop shared
set.seed(123)
nulls_b_drop_shared <- bipartite::nullmodel(b_metanetwork_drop_shared, method=5,
                                 nsim = 1000)

(con_nulls_b_drop_shared <- myFunwithNull(b_metanetwork_drop_shared,nulls_b_drop_shared,
                                                index = "connectance"))



(M_nulls_b_drop_shared <- MFunwithNull(b_metanetwork_drop_shared,nulls = nulls_b_drop_shared))


(nestedness_nulls_b_drop_shared <- myFunwithNull(b_metanetwork_drop_shared,
                                                       nulls = nulls_b_drop_shared,index = "NODF")) 




net_structure_b_equip <- rbind(con_nulls_b,M_nulls_b,nestedness_nulls_b,
                            con_nulls_b_drop_on,M_nulls_b_drop_on,nestedness_nulls_b_drop_on,
                            con_nulls_b_drop_under,M_nulls_b_drop_under,nestedness_nulls_b_drop_under,
                            con_nulls_b_drop_shared,M_nulls_b_drop_shared,nestedness_nulls_b_drop_shared)

net_structure_b_equip

write.csv(net_structure_b_equip,"./R2/output/net_structure_b_equip.csv")




#-------------------------------------------------------------------------------------
## 2 Interaction level analysis

## 2.1 quantitative analysis
#q_metanetwork <- t(q_metanetwork)

degree.interaction <- specieslevel(q_metanetwork,level = "lower",index = "degree")
degree.patch <- specieslevel(q_metanetwork,level = "higher",index = "degree")

BC.interaction_q <- BC(q_metanetwork,weighted = T)[1]
BC.patch_q <- BC(q_metanetwork,weighted = T)[2]

## 2.2 binary analysis
(BC.interaction_b <- BC(b_metanetwork,weighted = F)[1])
(BC.patch_b <- BC(b_metanetwork,weighted = T)[2])


## 3. Species interaction contribution to metanework structure
## 3.1 quantitative analysis

set.seed(123)
mod.q <- computeModules(q_metanetwork)
cz <- czvalues(mod.q,weighted=T,level="lower")
c.interaction <- cz$c
z.interaction <- cz$z


cz.patch <- czvalues(mod.q,weighted=T,level="higher")
c.patch <- cz.patch$c
z.patch <- cz.patch$z # NA produced

nestedcontribution(q_metanetwork, nsimul = 999) -> nestcontribution
save(nestcontribution,file = "nestcontribution20200123.Rdata")
load("nestcontribution20200123.Rdata")
cn_interaction <- nestcontribution$`lower level`
cn_patch <- nestcontribution$`higher level`

##perform a PCA analysis
contribution2metanet <- data.frame(degree.interaction,c.interaction,
                                   z.interaction,cn_interaction)

str(contribution2metanet)

contribution2metanet.pca <- prcomp(as.matrix(contribution2metanet), 
                                   center = TRUE,scale. = TRUE)

summary(contribution2metanet.pca)

contribution2metanet.pca$x

contribution2metanet.pca$rotation

write.csv(contribution2metanet,"./R2/output/contribution2metanet_p.csv")
write.csv(contribution2metanet.pca$x,"./R2/output/PCA_p.csv")



library(corrplot)
Mq <- cor(contribution2metanet,method = "spearman")
resq <- cor.mtest(contribution2metanet,method = "spearman",exact=FALSE)

corrplot(Mq, p.mat = resq$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white",
         type = "upper")





##3.2 binary analysis

set.seed(123)
mod.b <- computeModules(b_metanetwork)
cz.b <- czvalues(mod.b,weighted=F,level="lower")
c.interaction_b <- cz.b$c
z.interaction_b <- cz.b$z

nestedcontribution(b_metanetwork, nsimul = 999) -> nestcontribution #same as quantitative

##perform a PCA analysis
contribution2metanet.b <- data.frame(degree.interaction,c.interaction_b,
                                   z.interaction_b,cn_interaction)
pairs(contribution2metanet.b)

str(contribution2metanet.b)
library(corrplot)
Mb <- cor(contribution2metanet.b,method = "spearman")
resb <- cor.mtest(contribution2metanet.b,method = "spearman",exact=FALSE)

corrplot(Mb, p.mat = resb$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white",
         type = "lower")




contribution2metanet.b.pca <- prcomp(as.matrix(contribution2metanet.b), 
                                   center = TRUE,scale. = TRUE)

summary(contribution2metanet.b.pca)

contribution2metanet.b.pca$x


write.csv(contribution2metanet.b,"./R2/output/contribution2metanet_b.csv")
write.csv(contribution2metanet.b.pca$x,"./R2/output/PCA_b.csv")





## 5. Statistic analysis
##5.1 data visualization
Contribute_Bird <- read.xlsx("Outpute/Contribution.xlsx",sheet = 2)
names(Contribute_Bird)
Contribute_Plant <- read.xlsx("Outpute/Contribution.xlsx",sheet = 1)
names(Contribute_Plant)
Contribute_Interaction <- read.xlsx("Outpute/Contribution.xlsx",sheet = 3)
names(Contribute_Interaction)


aggregate(nestedcontribution ~ Methods, data = Contribute_Interaction, mean)
aggregate(nestedcontribution ~ Methods, data = Contribute_Interaction, sd)

aggregate(BSD_degree ~ Methods, data = Contribute_Interaction, mean)
aggregate(BSD_degree ~ Methods, data = Contribute_Interaction, sd)

aggregate(BC.Inter ~ Methods, data = Contribute_Interaction, mean)
aggregate(BC.Inter ~ Methods, data = Contribute_Interaction, sd)

summary(aov(log(BSD_degree) ~ Methods, data = Contribute_Interaction))
summary(aov(asin(sqrt(BC.Inter)) ~ Methods, data = Contribute_Interaction))

plot(aov(asin(sqrt(BC.Inter)) ~ Methods, data = Contribute_Interaction))

kruskal.test(log(BSD_degree) ~ factor(Methods), data = Contribute_Interaction)
kruskal.test(BC.Inter ~ factor(Methods), data = Contribute_Interaction)


###Interactions
Contribute_Interaction <- read.xlsx("./R2/output/metanet_analysis.xlsx",sheet = 3)
names(Contribute_Interaction)
Contribute_Interaction$Methods <- as.factor(Contribute_Interaction$Methods)

ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_degree),
                                  y=BSD_degree,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("BP-interactions")+
  ylab("degree of BP-interactions, k")


ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BC.Inter),
                                  y=BC.Inter,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("BP-interactions")+
  ylab("betweenness of BP-interactions, Bc")

ggplot(Contribute_Interaction,aes(x=PB_Interaction,
                                  y=BC.Inter,group=Methods,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("BP-interactions")+
  ylab("betweenness of BP-interactions, Bc")


##########For PC1
library(multcomp)
aggregate(PC1_q~Methods,data=Contribute_Interaction, mean)
aggregate(PC1_q~Methods,data=Contribute_Interaction, sd)

aggregate(PC1_b~Methods,data=Contribute_Interaction, mean)
aggregate(PC1_b~Methods,data=Contribute_Interaction, sd)


summary(aov(PC1_q~Methods,data=Contribute_Interaction))

fit <- aov(PC1_q~Methods,data=Contribute_Interaction)
summary(glht(fit, linfct=mcp(Methods="Tukey")),adjusted(type = "bonferroni"))

tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)

summary(aov(PC1_b~Methods,data=Contribute_Interaction))

fit <- aov(PC1_b~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)




PC1_q_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-PC1_q),
                          y=PC1_q,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's contribution \n to metanetwork structure")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
        axis.title.x = element_text(size = 12,family = "serif"),
        axis.title.y = element_text(size = 12,family = "serif"),
        legend.title = element_text(size = 12,family = "serif"),
        legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(A)")

PC1_q_bar




PC1_b_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-PC1_b),
                                               y=PC1_b,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's contribution \n to metanetwork structure")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(B)")

PC1_b_bar

PC1_q_bar+PC1_b_bar+plot_layout(ncol = 1)



PC1_q_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=PC1_q,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's contribution \n to metanetwork structure")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

PC1_q_boxplot

PC1_b_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=PC1_b,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's contribution \n to metanetwork structure")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

PC1_q_boxplot+PC1_b_boxplot+plot_layout(ncol = 2)
PC1_b_boxplot
#-------------------------------------------------------------------------------

##########For Ci
aggregate(BSD_C_q~Methods,data=Contribute_Interaction, mean)
aggregate(BSD_C_q~Methods,data=Contribute_Interaction, sd)

ggplot(Contribute_Interaction,aes(x=BSD_C_q,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()

aggregate(BSD_C_b~Methods,data=Contribute_Interaction, mean)
aggregate(BSD_C_b~Methods,data=Contribute_Interaction, sd)

summary(aov(log(BSD_C_q+0.00001)~Methods,data=Contribute_Interaction))

fit <- aov(BSD_C_q~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)

BSD_C_q_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_C_q),
                                               y=BSD_C_q,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's C value")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(A)")

BSD_C_q_bar

BSD_C_b_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_C_b),
                                               y=BSD_C_b,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's C value")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(B)")

BSD_C_b_bar

BSD_C_q_bar+BSD_C_b_bar+plot_layout(ncol = 1)



BSD_C_q_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=BSD_C_q,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's C value")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

BSD_C_q_boxplot

BSD_C_b_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=BSD_C_b,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's C value")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

BSD_C_b_boxplot

BSD_C_q_boxplot+BSD_C_b_boxplot+plot_layout(ncol = 2)


#-------------------------------------------------------------------------------

##########For Zi
aggregate(BSD_Z_q~Methods,data=Contribute_Interaction, mean)
aggregate(BSD_Z_q~Methods,data=Contribute_Interaction, sd)

ggplot(Contribute_Interaction,aes(x=BSD_Z_q,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()

aggregate(BSD_Z_b~Methods,data=Contribute_Interaction, mean)
aggregate(BSD_Z_b~Methods,data=Contribute_Interaction, sd)

summary(aov(BSD_Z_q~Methods,data=Contribute_Interaction))

fit <- aov(BSD_Z_q~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)

BSD_Z_q_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_Z_q),
                                                 y=BSD_Z_q,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's Z value")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(A)")

BSD_Z_q_bar

BSD_Z_b_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_Z_b),
                                                 y=BSD_Z_b,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's Z value")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))+
  labs(tag = "(B)")

BSD_Z_b_bar

BSD_Z_q_bar+BSD_Z_b_bar+plot_layout(ncol = 1)



BSD_Z_q_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=BSD_Z_q,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's Z value")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

BSD_Z_q_boxplot

BSD_Z_b_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=BSD_Z_b,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's Z value")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

BSD_Z_b_boxplot

BSD_Z_q_boxplot+BSD_Z_b_boxplot+plot_layout(ncol = 2)

BSD_C_q_boxplot+BSD_C_b_boxplot+
  BSD_Z_q_boxplot+BSD_Z_b_boxplot+
  plot_layout(ncol = 2)



#-------------------------------------------------------------------------------

##########For Cni
aggregate(BSD_Cn~Methods,data=Contribute_Interaction, mean)
aggregate(BSD_Cn~Methods,data=Contribute_Interaction, sd)

ggplot(Contribute_Interaction,aes(x=BSD_Cn,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()


summary(aov(BSD_Cn~Methods,data=Contribute_Interaction))

fit <- aov(BSD_Cn~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)

BSD_Cn_bar <- ggplot(Contribute_Interaction,aes(x=reorder(PB_Interaction,-BSD_Cn),
                                                 y=BSD_Cn,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Plant-Bird interactions")+
  ylab("Bird-plant interaction's contribution to metanetwork nestedness")+
  theme( axis.text.y= element_text(size = 12,family = "serif"), 
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))


BSD_Cn_bar



BSD_Cn_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=BSD_Cn,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("Bird-plant interaction's contribution \n to metanetwork nestedness")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))

BSD_Cn_boxplot


################## For eH and d'
aggregate(eH~Methods,data=Contribute_Interaction, mean)
aggregate(eH~Methods,data=Contribute_Interaction, sd)

ggplot(Contribute_Interaction,aes(x=eH,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()


summary(aov(eH~Methods,data=Contribute_Interaction))

fit <- aov(eH~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)


BSD_eH_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=eH,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("eH")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))
BSD_eH_boxplot


aggregate(d~Methods,data=Contribute_Interaction, mean)
aggregate(d~Methods,data=Contribute_Interaction, sd)

ggplot(Contribute_Interaction,aes(x=d,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()


summary(aov(d~Methods,data=Contribute_Interaction))

fit <- aov(d~Methods,data=Contribute_Interaction)
tuk <- glht(fit, linfct=mcp(Methods="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)


BSD_d_boxplot <- ggplot(Contribute_Interaction,aes(x=Methods,y=d,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab("d'")+
  xlab("")+
  theme( axis.text.x= element_text(size = 12,family = "serif"), 
         axis.text.y= element_text(size = 12,family = "serif"),
         axis.title.x = element_text(size = 12,family = "serif"),
         axis.title.y = element_text(size = 12,family = "serif"),
         legend.title = element_text(size = 12,family = "serif"),
         legend.text =element_text(size = 12,family = "serif"))
BSD_d_boxplot

BSD_eH_boxplot+BSD_d_boxplot+plot_layout(ncol = 2)



######Bird role vs. degree and BC
library(car)
library(multcomp)
Contribute_Interaction$BirdRole <- as.factor(Contribute_Interaction$BirdRole)
summary(aov(BSD_degree~BirdRole,data=Contribute_Interaction))
#cor.test(Contribute_Interaction$BSD_degree,Contribute_Interaction$BirdRole,
#         method = "spearman")
fit <- aov(log(BSD_degree)~BirdRole,data=Contribute_Interaction)
summary(fit)
summary(glht(fit,linfct=mcp(BirdRole="Tukey")),adjusted(type = "bonferroni"))


tuk <- glht(fit, linfct = mcp(BirdRole="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk, level=.05),col="lightgrey")
par(op)

op <- par(mfrow=c(2,2))
plot(fit)
par(op)

qqPlot(fit)
boxplot(BSD_degree~BirdRole,data=Contribute_Interaction)

summary(aov(BC.Inter~BirdRole,data=Contribute_Interaction))
#cor.test(Contribute_Interaction$BSD_degree,Contribute_Interaction$BirdRole,
#         method = "spearman")
fit.bc <- aov(log(BC.Inter+0.000001)~BirdRole,data=Contribute_Interaction)
summary(fit.bc)

summary(glht(fit.bc,linfct=mcp(BirdRole="Tukey")),adjusted(type = "bonferroni"))

tuk.bc <- glht(fit.bc, linfct = mcp(BirdRole="Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk.bc, level=.05),col="lightgrey")
par(op)

op <- par(mfrow=c(1,2))
plot(cld(tuk, level=.05),col="lightgrey")
plot(cld(tuk.bc, level=.05),col="lightgrey")
par(op)


op <- par(mfrow=c(2,2))
plot(fit)
par(op)

qqPlot(fit.bc)

boxplot(BC.Inter~BirdRole,data=Contribute_Interaction)


###Species traits~role
Birdtraits <- read.xlsx("./R2/output/metanet_analysis.xlsx",sheet = 2)
names(Birdtraits)
summary(aov(BodyMass~BirdRole,data = Birdtraits))
boxplot(BodyMass~BirdRole,data = Birdtraits)
summary(aov(degree~BirdRole,data = Birdtraits))
boxplot(degree~BirdRole,data = Birdtraits)


Planttraits <- read.xlsx("./R2/output/metanet_analysis.xlsx",sheet = 1)
names(Planttraits)
summary(aov(SeedMass2~PlantRole,data = Planttraits))
boxplot(SeedMass2~PlantRole,data = Planttraits)
summary(aov(degree~PlantRole,data = Planttraits))
boxplot(degree~PlantRole,data = Planttraits)
########Species

ggplot(Contribute_Bird,aes(x=reorder(Bird,-nestedcontribution),
                           y=nestedcontribution,fill=Methods))+
  geom_bar(stat = "identity")+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Bird")+
  ylab("Species' contribution to nestedness")

ggplot(Contribute_Bird,aes(y=nestedcontribution,fill=Methods))+
  geom_boxplot()+scale_fill_manual(values =c("blue","black","red"))+
  theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())+xlab("Bird")+
  ylab("Species' contribution to nestedness")


ggplot(Contribute_Bird,aes(x=nestedcontribution,color=Methods))+
  geom_density(alpha=.1)+theme_bw()+scale_color_manual(values =c("blue","black","red"))+
  theme_bw()+xlab("Species' contribution to nestedness")



############Due to rejected the SEM models in above, I just use spearman correlation

##In species level, calculating as following: 
## Bird degree~ Body mass
cor.test(Centrity_Bird$degree,Centrity_Bird$BodyMass,method = "spearman")
cor.test(Centrity_Bird$degree,Centrity_Bird$BillLength,method = "spearman")

## Bird Bc~ Body mass
cor.test(Centrity_Bird$BC.Bird,Centrity_Bird$BodyMass,method = "spearman")
cor.test(Centrity_Bird$BC.Bird,Centrity_Bird$BillLength,method = "spearman")


## Plant degree~ seed mass
cor.test(Centrity_Plant$degree,Centrity_Plant$FruitMass,method = "spearman")
cor.test(Centrity_Plant$degree,Centrity_Plant$FruitDiameter,method = "spearman")
cor.test(Centrity_Plant$degree,Centrity_Plant$SeedMass,method = "spearman")


## Plant Bc~ fruit mass
cor.test(Centrity_Plant$BC.Plant,Centrity_Plant$FruitMass,method = "spearman")
cor.test(Centrity_Plant$BC.Plant,Centrity_Plant$FruitDiameter,method = "spearman")


###

Planttraits
cor.test(Planttraits$FruitDiameter,Planttraits$SeedMass2)
cor.test(Planttraits$SeedDiameter,Planttraits$SeedMass2)
sort(Planttraits$SeedMass2)

## Interaction level
rand<-rnorm(134,sd=1e-10)

library(corrplot)
## corrplot 0.84 loaded

cordat <- Contribute_Interaction[,c("BSD_degree","BC.Inter","PC1_q","PC1_b","BSD_Cn")]
M <- cor(cordat,method = "spearman")
res1 <- cor.mtest(cordat,method = "spearman",exact=FALSE)

corrplot(M, p.mat = res1$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white",
         type = "upper")

corrplot(M, type = "upper")


names(Contribute_Interaction)
cordat2 <-Contribute_Interaction[,c("BSD_degree","BC.Inter","PC1_q","PC1_b","BSD_Cn",
                                    "Plant_degree","SeedMass2","Bird_degree","BodyMass")]
cordat2 <- na.omit(cordat2)

M2 <- cor(cordat2,method = "spearman")
res2 <- cor.mtest(cordat2,method = "spearman",exact=FALSE)

cor(cordat2+rnorm(nrow(cordat2),sd=1e-10),method = "spearman")
cor.mtest(cordat2+rnorm(nrow(cordat2),sd=1e-10),method = "spearman")[[1]]

corrplot(M2, p.mat = res2$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white",
         type = "upper")

## add all p-values
corrplot(M2, p.mat = res2$p, insig = "p-value", sig.level = -1)

corrplot(M2, method = "number")

cor.test(Contribute_Interaction$BSD_degree,Contribute_Interaction$BodyMass,
         method = "spearman")
cor.test(Centrity_Interaction$BSD_degree,Centrity_Interaction$BillLength,
         method = "spearman")
cor.test(Centrity_Interaction$BSD_degree,Centrity_Interaction$FruitMass,
         method = "spearman")

cor.test(Contribute_Interaction$BC.Inter+rand,Contribute_Interaction$BodyMass+rand,
         method = "spearman")
cor.test(Centrity_Interaction$BC.Inter,Centrity_Interaction$BillLength,
         method = "spearman")
cor.test(Centrity_Interaction$BC.Inter,Centrity_Interaction$FruitDiameter,
         method = "spearman")

cor.test(Contribute_Interaction$BC.Inter+rand,(Contribute_Interaction$SeedMass+rand),
         method = "spearman")
cor.test(Contribute_Interaction$BC.Inter,Contribute_Interaction$SeedMass,
         method = "spearman",exact = FALSE)



plot(Centrity_Interaction$BC.Inter~log(Centrity_Interaction$SeedMass))


## 6. beta-diversity analysis
#beta diversity analysis
# get betapart objects
metanet_matrix.core <- betapart.core(b_metanetwork)

# multiple site measures
metanet.multi <- beta.multi(metanet_matrix.core)
metanet.multi

# sampling across equal sites
metanet.samp <- beta.sample(metanet_matrix.core,
                            sites=10, samples=100)

# plotting the distributions of components
dist.metanet  <- metanet.samp$sampled.values

plot(density(dist.metanet$beta.SOR),
     xlim=c(0,1), ylim=c(0, 100), xlab='Beta diversity', main='', lwd=3)
lines(density(dist.metanet$beta.SNE), lty=1, lwd=2)
lines(density(dist.metanet$beta.SIM), lty=2, lwd=2)

# pairwise for all patch
pair.metanet <- beta.pair(metanetwork)

as.matrix(pair.metanet$beta.sor)
as.matrix(pair.metanet$beta.sim)
as.matrix(pair.metanet$beta.sne)

write.csv(as.matrix(pair.metanet$beta.sor),"Interaction_beta_SOR.csv")
write.csv(as.matrix(pair.metanet$beta.sim),"Interaction_beta_SIM.csv")
write.csv(as.matrix(pair.metanet$beta.sne),"Interaction_beta_SNE.csv")




dist2f <- function(dist){
    require(reshape2)
    df <-melt(as.matrix(dist), varnames = c("i", "j"))
    df <- df[as.numeric(df$i) > as.numeric(df$j),]
    return(df)
}

dist2f(pair.metanet$beta.sne)

interaction_beta <- merge(dist2f(pair.metanet$beta.sne),
                          dist2f(pair.metanet$beta.sim),by=c("i","j"))

interaction_beta <- merge(dist2f(pair.metanet$beta.sor),
                          interaction_beta,by=c("i","j"))

names(interaction_beta)[3:5]<-c("beta_SOR","beta_SNE","beta_SIM")



Area=c(BRS=80.45, 
       NNC=41.89, 
       SJHS=6.05, 
       SZZ=17.63, 
       WJQS=22.99, 
       XJHS=2.68, 
       XJY=57.51, 
       XE1=20.23, 
       XE3=16.18, 
       XE4=3.75, 
       XE5=2.85, 
       XiaoJHS=4.20, 
       XXHS=5.23) 

dist(Area)
area_dist_df <- dist2f(dist(Area))

MRM(pair.metanet$beta.sne~dist(Area))
Geo_dist <- read.csv("patchdistance_center.csv",row.names = 1)
Geo_dist <- Geo_dist[names(Area),names(Area)]
geo_dis_df <- dist2f(as.dist(Geo_dist))


dist_data <- merge(area_dist_df,geo_dis_df,by=c("i","j"))
names(dist_data)[3:4]<- c("area","geo_dist")

interaction_beta <- merge(interaction_beta,dist_data)
write.csv(interaction_beta,"interaction_beta_betapart.csv",row.names = F)

interaction_beta <- read.csv("interaction_beta_betapart.csv")

mean(interaction_beta$beta_SOR)
#0.7599523
sd(interaction_beta$beta_SOR)
#0.1206006

mean(interaction_beta$beta_SIM)
#0.6872141
sd(interaction_beta$beta_SIM)
#0.1421138

mean(interaction_beta$beta_SNE)
#0.07273826
sd(interaction_beta$beta_SNE)
#0.0643046


library(reshape2)
interaction_beta_2 <-melt(interaction_beta[,c(1:5)],id=c("i","j"))
names(interaction_beta_2)[3:4] <- c("Partition","Turnover")

Fig.E <- ggplot(interaction_beta,aes(x=log(geo_dist),y=beta_SOR))+
    geom_point(size=2)+
    geom_point(aes(x=log(geo_dist),y=beta_SIM),shape=0,size=2)+
    geom_point(aes(x=log(geo_dist),y=beta_SNE),shape=2,size=2)+
    geom_smooth(aes(x=log(geo_dist),y=beta_SOR),method = "lm",col="black")+
    geom_smooth(aes(x=log(geo_dist),y=beta_SIM),
                method = "lm",col="blue")+
    geom_smooth(aes(x=log(geo_dist),y=beta_SNE),
                method = "lm",col="red",linetype="dashed")+
    theme_bw()+ theme(panel.grid=element_blank())+
    xlab("Log(geographic distance)")+ylab("Interaction beta diversity")+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    labs(tag="(E)")

Fig.E



Fig.F=ggplot(data = interaction_beta_2,
             aes(x=factor(Partition,levels = c("beta_SOR","beta_SIM","beta_SNE")),
                 y=Turnover,fill=Partition))+
    geom_boxplot()+
    scale_color_manual(values = c("black","red","blue"))+
    scale_fill_manual(values = c("black","red","blue"))+
    xlab("")+ylab("Interaction beta diversity")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    labs(tag = "(F)")

Fig.F






summary(lm(beta_SOR~log(area)+log(geo_dist),data = interaction_beta))
summary(lm(beta_SIM~log(area)+log(geo_dist),data = interaction_beta))
summary(lm(beta_SNE~log(area)+log(geo_dist),data = interaction_beta))


set.seed(123)
summary(lmp(beta_SOR~area+geo_dist,data = interaction_beta, perm="Prob"))
summary(lmp(beta_SIM~area+geo_dist,data = interaction_beta, perm="Prob"))
summary(lmp(beta_SNE~area+geo_dist,data = interaction_beta, perm="Prob"))


set.seed(123)
MRM(pair.metanet$beta.sor~dist(Area)+as.dist(Geo_dist),nperm = 1000)
MRM(pair.metanet$beta.sim~dist(Area)+as.dist(Geo_dist),nperm = 1000)
MRM(pair.metanet$beta.sne~dist(Area)+as.dist(Geo_dist),nperm = 1000)


plot(mgram(pair.metanet$beta.sor, as.dist(Geo_dist), nclass=8))

log(as.dist(Geo_dist))

x=dist(Area,diag = T,upper = T)

x<- as.matrix(x)

rownames(x)==names(Area)



Beta_Species <- read.xlsx("Beta_Species.xlsx",sheet=1)
mean(Beta_Species$Beta_SOR_Plant);sd(Beta_Species$Beta_SOR_Plant)
#[1] 0.3385385
#[1] 0.1129542
mean(Beta_Species$Beta_SIM_Plant);sd(Beta_Species$Beta_SIM_Plant)
#[1] 0.2542125
#[1] 0.136165
mean(Beta_Species$Beta_SNE_Plant);sd(Beta_Species$Beta_SNE_Plant)
#[1] 0.08432606
#[1] 0.07044583

mean(Beta_Species$Beta_SIM_Bird);sd(Beta_Species$Beta_SIM_Bird)
#[1] 0.4455597
#[1] 0.1418905

mean(Beta_Species$Beta_SNE_Bird);sd(Beta_Species$Beta_SNE_Bird)
#[1] 0.1102871
#[1] 0.105618


fit1 <- lm(Beta_SOR_Plant~log(Distance_C)+log(Distance_area),data = Beta_Species)
fit2 <- lm(Beta_SIM_Plant~log(Distance_C)+log(Distance_area),data = Beta_Species)
fit3 <- lm(Beta_SNE_Plant~log(Distance_C)+log(Distance_area),data = Beta_Species)

summary(fit1)
summary(fit2)
summary(fit3)

summary(step(fit1))
summary(step(fit2))
summary(step(fit3))



plt.beta.plant.geo <- ggplot(Beta_Species,aes(x=log(Distance_C),y=Beta_SOR_Plant))+
    geom_point(size=2)+
    geom_point(aes(x=log(Distance_C),y=Beta_SIM_Plant),shape=0,size=2)+
    geom_point(aes(x=log(Distance_C),y=Beta_SNE_Plant),shape=2,size=2)+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SOR_Plant),method = "lm",col="black")+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SIM_Plant),
                method = "lm",col="blue",linetype="dashed")+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SNE_Plant),
                method = "lm",col="red",linetype="dashed")+
    theme_bw()+ theme(panel.grid=element_blank())+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    xlab("Log(geographic distance)")+ylab("Plant species beta diversity")+
    labs(tag="(A)")

plt.beta.plant.geo


Beta_Species_Plant <-melt(Beta_Species[,c(1:5)],id=c("i","j"))
names(Beta_Species_Plant)[3:4] <- c("Partition","Turnover")

Fig.B=ggplot(data = Beta_Species_Plant,
             aes(x=factor(Partition,levels = c("Beta_SOR_Plant","Beta_SIM_Plant","Beta_SNE_Plant")),
                 y=Turnover,fill=Partition))+
    geom_boxplot()+
    scale_color_manual(values = c("black","red","blue"))+
    scale_fill_manual(values = c("black","red","blue"))+
    xlab("")+ylab("Plant species beta diversity")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    labs(tag = "(B)")

Fig.B





fit1.bird <- lm(Beta_SOR_Bird~log(Distance_C)+log(Distance_area),data = Beta_Species)
fit2.bird <- lm(Beta_SIM_Bird~log(Distance_C)+log(Distance_area),data = Beta_Species)
fit3.bird <- lm(Beta_SNE_Bird~log(Distance_C)+log(Distance_area),data = Beta_Species)

summary(fit1.bird)
summary(fit2.bird)
summary(fit3.bird)

summary(step(fit1.bird))
summary(step(fit2.bird))
summary(step(fit3.bird))






plt.beta.bird.geo <- ggplot(Beta_Species,aes(x=log(Distance_C),y=Beta_SOR_Bird))+
    geom_point(size=2)+
    geom_point(aes(x=log(Distance_C),y=Beta_SIM_Bird),shape=0,size=2)+
    geom_point(aes(x=log(Distance_C),y=Beta_SNE_Bird),shape=2,size=2)+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SOR_Bird),
                method = "lm",col="black",linetype="dashed")+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SIM_Bird),
                method = "lm",col="blue",linetype="dashed")+
    geom_smooth(aes(x=log(Distance_C),y=Beta_SNE_Bird),
                method = "lm",col="red",linetype="dashed")+
    theme_bw()+ theme(panel.grid=element_blank())+
    xlab("Log(geographic distance)")+ylab("Bird species beta diversity")+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    labs(tag="(C)")
plt.beta.bird.geo


Beta_Species_Bird <-melt(Beta_Species[,c(1:2,6:8)],id=c("i","j"))
names(Beta_Species_Bird)[3:4] <- c("Partition","Turnover")

Fig.D=ggplot(data = Beta_Species_Bird,
             aes(x=factor(Partition,levels = c("Beta_SOR_Bird","Beta_SIM_Bird","Beta_SNE_Bird")),
                 y=Turnover,fill=Partition))+
    geom_boxplot()+
    scale_color_manual(values = c("black","red","blue"))+
    scale_fill_manual(values = c("black","red","blue"))+
    xlab("")+ylab("Bird species beta diversity")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    theme(axis.title.y=element_text(size = 12,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"))+
    labs(tag = "(D)")

Fig.D

plt.beta.plant.geo+Fig.B+plt.D

## 7. Multilayer network analysis
########################################
##########Multilayer network analysis



# this script shows how to calculate the {z, c}-values of species observed 
# interacting within a landscape (a set of habitats)


# load functions, code from Hackett et al. 2019, pls cite this.
source("Functions/sumbygroup.r")
source("Functions/mw_adjacencies.R")
source("Functions/network_role.R")
source("Functions/feature_scaling.R")

# load data
# interaction list for the Dujiangyan

df_network_two <- rbind(df_network2014,df_network2015)

Dujiangyan_ListObs <- sumbygroup(df_network_two,measurevar = "Animal_Number",groupvars = c("Patch","TreeSpecies","AnimalSpecies"))

write.csv(Dujiangyan_ListObs,"Dujiangyan_ListObs.csv",row.names = F)

HH_ListObs <- read.csv("Dujiangyan_ListObs.csv", stringsAsFactors = FALSE)

# normalise interaction frequencies
IntTypes <- unique(HH_ListObs$upper_guild)
HH_ListObs$freq_norm <- 0
for (it in IntTypes){
    freq_int <- HH_ListObs$freq[HH_ListObs$upper_guild == it]
    HH_ListObs$freq_norm[HH_ListObs$upper_guild == it] <- round(feature_scaling(freq_int, min(HH_ListObs$freq), max(HH_ListObs$freq)), 3) # feature scaling, 3 digits for fractional part
}

# species list for the mosaic of habitats in Hengistbury Head
# must be a data frame with at least one column with species names ('taxon'), it can be customised to your taste with more information
HH_SpList <- data.frame(taxon = unique(sort(c(HH_ListObs$lower_taxon, HH_ListObs$upper_taxon))))

# metaweb, all interaction types
a_obs <- mw_adjacencies(HH_ListObs, HH_SpList, FALSE, FALSE, TRUE) # 3D-array, dim = {Sp+Sa, Sp+Sa, H}, Sp = No. plant species, Sa = No. insect species, H = No. habitats
aw_obs <- mw_adjacencies(HH_ListObs, HH_SpList, TRUE, FALSE, TRUE) # weighted

# Among-habitat connectivity
C <- c_value(a_obs) # qualitative
Cw <- c_value(aw_obs) # quantitative

# Weighted mean of within-habitat degree/weight
Z <- z_value_landscape(a_obs) # qualitative
Zw <- z_value_landscape(aw_obs) # quantitative





## 8.Rarefiction analysis

sapply(Networks,FUN = sum)
metanetwork_df <- read.csv("metanetwork_dataframe.csv",header = T)
head(metanetwork_df)

metanetwork_df1 <-droplevels(metanetwork_df[which(metanetwork_df$patch!="XXHS"),])
metanetwork_df2 <-droplevels(metanetwork_df1[which(metanetwork_df1$patch!="XE5"),])


rarefied_metanetwork <- function(df,patch="patch",plant="Plant",
                                 animal="Animal",Value="feq",guild="guild",
                                 rows=20,nsim=100,on="on",under="under",
                                 shared="shared",quantitative=TRUE){
    patchname <- unique(df[,patch])
    full_webs <- list()
    drop_on <- list()
    drop_under <-list()
    drop_shared <-list()
    #test <-list()
    for (i in 1:nsim) {
        dflist <-list()
        for (j in 1:length(patchname)) {
            tep <- df[df[,patch] %in% patchname[j],]
            tep[sample(1:length(tep[,1]),rows),] ->dflist[[j]]
        }
        df_full <- do.call("rbind",dflist)
        df_drop_on <-droplevels(df_full[which(df_full[,guild]!=on),])
        df_drop_under <-droplevels(df_full[which(df_full[,guild]!=under),])
        df_drop_shared <-droplevels(df_full[which(df_full[,guild]!=shared),])
        #test[[i]] <-df_drop_shared
        full_webs[[i]] <- df2metanetwork(df=df_full,patch=patch,
                                         plant=plant,animal=animal,
                       Value=Value,quantitative=quantitative)
        drop_on[[i]] <- df2metanetwork(df=df_drop_on,patch=patch,
                                       plant=plant,animal=animal,
                                      Value=Value,quantitative=quantitative)
        drop_under[[i]] <- df2metanetwork(df=df_drop_under,patch=patch,
                                         plant=plant,animal=animal,
                                         Value=Value,quantitative=quantitative)
        drop_shared[[i]] <- df2metanetwork(df=df_drop_shared,patch=patch,
                                            plant=plant,animal=animal,
                                          Value=Value,quantitative=quantitative)
    }
    return(list(full_webs = full_webs,drop_on=drop_on,
                drop_under= drop_under,drop_shared=drop_shared))
    
}




rarefied_fixed_species_and_interactions <- function(df,patch="patch",plant="Plant",
                                 animal="Animal",Value="feq",guild="guild",
                                 rows=20,nsim=100,on="on",under="under",
                                 shared="shared",quantitative=TRUE){
    patchname <- unique(df[,patch])
    full_webs <- list()
    drop_on <- list()
    drop_under <-list()
    drop_shared <-list()
    #test <-list()
    for (i in 1:nsim) {
        dflist <-list()
        for (j in 1:length(patchname)) {
            tep <- df[df[,patch] %in% patchname[j],]
            plant_n <- length(unique(tep[,plant]))
            animal_n <- length(unique(tep[,animal]))
            success <- FALSE
            while (!success) {
                #sample process
                tep2 <- tep[sample(1:length(tep[,1]),rows),]
                plant_n2 <- length(unique(tep2[,plant]))
                animal_n2 <- length(unique(tep2[,animal]))
                
                # check for success
                success <- plant_n2==plant_n & animal_n==animal_n2
            }
            tep2->dflist[[j]]
        }
        df_full <- do.call("rbind",dflist)
        df_drop_on <-droplevels(df_full[which(df_full[,guild]!=on),])
        df_drop_under <-droplevels(df_full[which(df_full[,guild]!=under),])
        df_drop_shared <-droplevels(df_full[which(df_full[,guild]!=shared),])
        #test[[i]] <-df_drop_shared
        full_webs[[i]] <- df2metanetwork(df=df_full,patch=patch,
                                         plant=plant,animal=animal,
                                         Value=Value,quantitative=quantitative)
        drop_on[[i]] <- df2metanetwork(df=df_drop_on,patch=patch,
                                       plant=plant,animal=animal,
                                       Value=Value,quantitative=quantitative)
        drop_under[[i]] <- df2metanetwork(df=df_drop_under,patch=patch,
                                          plant=plant,animal=animal,
                                          Value=Value,quantitative=quantitative)
        drop_shared[[i]] <- df2metanetwork(df=df_drop_shared,patch=patch,
                                           plant=plant,animal=animal,
                                           Value=Value,quantitative=quantitative)
    }
    return(list(full_webs = full_webs,drop_on=drop_on,
                drop_under= drop_under,drop_shared=drop_shared))
    
}






rarefied_metanetwork1=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                 animal="animal",Value="feq",guild="guild",
                                 rows=29,nsim=1000,on="on",under="under",
                                 shared="shared",quantitative=TRUE)
    

rarefied_metanetwork2=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                           animal="animal",Value="feq",guild="guild",
                                           rows=91,nsim=1000,on="on",under="under",
                                           shared="shared",quantitative=TRUE)



rarefied_meta_fixed_species(metanetwork_df2,patch="patch",plant="plant",
                            animal="animal",Value="feq",guild="guild",
                            rows=91,nsim=1,on="on",under="under",
                            shared="shared",quantitative=TRUE)

M_rarefied1_all <- sapply(rarefied_metanetwork1$full_webs,function(x)computeModules(x, method="Beckett")@likelihood)
mean(M_rarefied1_all)
#0.5423321
sd(M_rarefied1_all)
#0.01876729

M_rarefied2_all <- sapply(rarefied_metanetwork2$full_webs,function(x)computeModules(x, method="Beckett")@likelihood)
mean(M_rarefied2_all)
#0.5060327
sd(M_rarefied2_all)
#0.01079944


multi_modular <- function(list,size=NULL){
    M<-sapply(list,function(x)computeModules(x, method="Beckett")@likelihood)
    res<-data.frame(mean(M),sd(M))
    names(res) <- c("mean","sd")
    res$samplesize <-size
    return(res)
}

###sample size gradients analysis
#####################################################################################
#####remove two patches
metanetwork2=df2metanetwork(metanetwork_df2,patch="patch",plant="plant",
                            animal="animal",Value="feq",quantitative=TRUE)

(computeModules(metanetwork2, method="Beckett")@likelihood)#0.4447518

rarefied_metanetwork2_20=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                           animal="animal",Value="feq",guild="guild",
                                           rows=20,nsim=1000,on="on",under="under",
                                           shared="shared",quantitative=TRUE)
rarefied_metanetwork2_30=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=30,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)

rarefied_metanetwork2_40=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=40,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)

rarefied_metanetwork2_50=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=50,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork2_60=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=60,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork2_70=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=70,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork2_80=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=80,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork2_90=rarefied_metanetwork(metanetwork_df2,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=90,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)


M20=multi_modular(rarefied_metanetwork2_20$full_webs,size = 20)
M30=multi_modular(rarefied_metanetwork2_30$full_webs,size = 30)
M40=multi_modular(rarefied_metanetwork2_40$full_webs,size = 40)
M50=multi_modular(rarefied_metanetwork2_50$full_webs,size = 50)
M60=multi_modular(rarefied_metanetwork2_60$full_webs,size = 60)
M70=multi_modular(rarefied_metanetwork2_70$full_webs,size = 70)
M80=multi_modular(rarefied_metanetwork2_80$full_webs,size = 80)
M90=multi_modular(rarefied_metanetwork2_90$full_webs,size = 90)

Res_resampling_modularity <- do.call("rbind",list(M20,M30,M40,M50,M60,M70,M80,M90,
                                                  data.frame(mean=0.45,sd=NaN,samplesize=NaN)))

write.csv(Res_resampling_modularity,"./R2/output/Res_resampling_modularity.csv")


Res_resampling_modularity$zscore <- (0.45-Res_resampling_modularity$mean)/Res_resampling_modularity$sd

ggplot(Res_resampling_modularity,aes(x=samplesize,y=mean))+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
    geom_smooth() +
    geom_point(size=4)+geom_hline(yintercept = 0.45)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.title.y=element_text(size = 14,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"),
          strip.text.x = element_text(size = 12,family = "serif",face = "bold"))



ggplot(Res_resampling_modularity,aes(x=samplesize,y=zscore))+
    geom_smooth() +
    geom_point(size=4)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.title.y=element_text(size = 14,family = "serif"),
          axis.text.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          axis.title.x=element_text(size = 12,family = "serif"),
          legend.text =element_text(size = 12,family = "serif"),
          legend.title = element_text(size = 12,family = "serif"),
          strip.text.x = element_text(size = 12,family = "serif",face = "bold"))

###################drop on the tree bird involved links

M20_drop_on=multi_modular(rarefied_metanetwork2_20$drop_on,size = 20)
M30_drop_on=multi_modular(rarefied_metanetwork2_30$drop_on,size = 30)
M40_drop_on=multi_modular(rarefied_metanetwork2_40$drop_on,size = 40)
M50_drop_on=multi_modular(rarefied_metanetwork2_50$drop_on,size = 50)
M60_drop_on=multi_modular(rarefied_metanetwork2_60$drop_on,size = 60)
M70_drop_on=multi_modular(rarefied_metanetwork2_70$drop_on,size = 70)
M80_drop_on=multi_modular(rarefied_metanetwork2_80$drop_on,size = 80)
M90_drop_on=multi_modular(rarefied_metanetwork2_90$drop_on,size = 90)

Res_resampling_modularity_drop_on <- do.call("rbind",list(M20_drop_on,M30_drop_on,M40_drop_on,
                                                  M50_drop_on,M60_drop_on,M70_drop_on,
                                                  M80_drop_on,M90_drop_on))
write.csv(Res_resampling_modularity_drop_on,"./R2/output/Res_resampling_modularity_drop_on.csv")


###################drop under the tree bird involved links
M20_drop_under=multi_modular(rarefied_metanetwork2_20$drop_under,size = 20)
M30_drop_under=multi_modular(rarefied_metanetwork2_30$drop_under,size = 30)
M40_drop_under=multi_modular(rarefied_metanetwork2_40$drop_under,size = 40)
M50_drop_under=multi_modular(rarefied_metanetwork2_50$drop_under,size = 50)
M60_drop_under=multi_modular(rarefied_metanetwork2_60$drop_under,size = 60)
M70_drop_under=multi_modular(rarefied_metanetwork2_70$drop_under,size = 70)
M80_drop_under=multi_modular(rarefied_metanetwork2_80$drop_under,size = 80)
M90_drop_under=multi_modular(rarefied_metanetwork2_90$drop_under,size = 90)

Res_resampling_modularity_drop_under <- do.call("rbind",list(M20_drop_under,M30_drop_under,M40_drop_under,
                                                  M50_drop_under,M60_drop_under,M70_drop_under,
                                                  M80_drop_under,M90_drop_under))
write.csv(Res_resampling_modularity_drop_under,"./R2/output/Res_resampling_modularity_drop_under.csv")


###################drop shared bird involved links
M20_drop_shared=multi_modular(rarefied_metanetwork2_20$drop_shared,size = 20)
M30_drop_shared=multi_modular(rarefied_metanetwork2_30$drop_shared,size = 30)
M40_drop_shared=multi_modular(rarefied_metanetwork2_40$drop_shared,size = 40)
M50_drop_shared=multi_modular(rarefied_metanetwork2_50$drop_shared,size = 50)
M60_drop_shared=multi_modular(rarefied_metanetwork2_60$drop_shared,size = 60)
M70_drop_shared=multi_modular(rarefied_metanetwork2_70$drop_shared,size = 70)
M80_drop_shared=multi_modular(rarefied_metanetwork2_80$drop_shared,size = 80)
M90_drop_shared=multi_modular(rarefied_metanetwork2_90$drop_shared,size = 90)

Res_resampling_modularity_drop_shared <- do.call("rbind",list(M20_drop_shared,M30_drop_shared,M40_drop_shared,
                                                  M50_drop_shared,M60_drop_shared,M70_drop_shared,
                                                  M80_drop_shared,M90_drop_shared))

write.csv(Res_resampling_modularity_drop_shared,"./R2/output/Res_resampling_modularity_drop_shared.csv")


############################################################################################################

###removed one patch only

metanetwork1=df2metanetwork(metanetwork_df1,patch="patch",plant="plant",
                            animal="animal",Value="feq",quantitative=TRUE)

(computeModules(metanetwork1, method="Beckett")@likelihood)#0.45

rarefied_metanetwork1_20=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=20,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork1_30=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=30,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)

rarefied_metanetwork1_40=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=40,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)

rarefied_metanetwork1_50=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=50,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork1_60=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=60,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork1_70=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=70,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork1_80=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=80,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)
rarefied_metanetwork1_90=rarefied_metanetwork(metanetwork_df1,patch="patch",plant="plant",
                                              animal="animal",Value="feq",guild="guild",
                                              rows=90,nsim=1000,on="on",under="under",
                                              shared="shared",quantitative=TRUE)


M20_1=multi_modular(rarefied_metanetwork1_20$full_webs,size = 20)
M30_1=multi_modular(rarefied_metanetwork1_30$full_webs,size = 30)
M40_1=multi_modular(rarefied_metanetwork1_40$full_webs,size = 40)
M50_1=multi_modular(rarefied_metanetwork1_50$full_webs,size = 50)
M60_1=multi_modular(rarefied_metanetwork1_60$full_webs,size = 60)
M70_1=multi_modular(rarefied_metanetwork1_70$full_webs,size = 70)
M80_1=multi_modular(rarefied_metanetwork1_80$full_webs,size = 80)
M90_1=multi_modular(rarefied_metanetwork1_90$full_webs,size = 90)

Res_resampling_modularity1 <- do.call("rbind",list(M20_1,M30_1,M40_1,M50_1,M60_1,M70_1,M80_1,M90_1))

write.csv(Res_resampling_modularity1,"./R2/output/Res_resampling_modularity1.csv")
ggplot(Res_resampling_modularity1,aes(x=samplesize,y=mean))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_smooth() +
  geom_point(size=4)+geom_hline(yintercept = 0.45)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y=element_text(size = 14,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x=element_text(size = 12,family = "serif"),
        legend.text =element_text(size = 12,family = "serif"),
        legend.title = element_text(size = 12,family = "serif"),
        strip.text.x = element_text(size = 12,family = "serif",face = "bold"))


###################drop on the tree bird involved links

M20_drop_on1=multi_modular(rarefied_metanetwork1_20$drop_on,size = 20)
M30_drop_on1=multi_modular(rarefied_metanetwork1_30$drop_on,size = 30)
M40_drop_on1=multi_modular(rarefied_metanetwork1_40$drop_on,size = 40)
M50_drop_on1=multi_modular(rarefied_metanetwork1_50$drop_on,size = 50)
M60_drop_on1=multi_modular(rarefied_metanetwork1_60$drop_on,size = 60)
M70_drop_on1=multi_modular(rarefied_metanetwork1_70$drop_on,size = 70)
M80_drop_on1=multi_modular(rarefied_metanetwork1_80$drop_on,size = 80)
M90_drop_on1=multi_modular(rarefied_metanetwork1_90$drop_on,size = 90)

Res_resampling_modularity_drop_on1 <- do.call("rbind",list(M20_drop_on1,M30_drop_on1,M40_drop_on1,
                                                          M50_drop_on1,M60_drop_on1,M70_drop_on1,
                                                          M80_drop_on1,M90_drop_on1))

write.csv(Res_resampling_modularity_drop_on1,"./R2/output/Res_resampling_modularity_drop_on1.csv")


###################drop under the tree bird involved links
M20_drop_under1=multi_modular(rarefied_metanetwork1_20$drop_under,size = 20)
M30_drop_under1=multi_modular(rarefied_metanetwork1_30$drop_under,size = 30)
M40_drop_under1=multi_modular(rarefied_metanetwork1_40$drop_under,size = 40)
M50_drop_under1=multi_modular(rarefied_metanetwork1_50$drop_under,size = 50)
M60_drop_under1=multi_modular(rarefied_metanetwork1_60$drop_under,size = 60)
M70_drop_under1=multi_modular(rarefied_metanetwork1_70$drop_under,size = 70)
M80_drop_under1=multi_modular(rarefied_metanetwork1_80$drop_under,size = 80)
M90_drop_under1=multi_modular(rarefied_metanetwork1_90$drop_under,size = 90)

Res_resampling_modularity_drop_under1 <- do.call("rbind",list(M20_drop_under1,M30_drop_under1,M40_drop_under1,
                                                             M50_drop_under1,M60_drop_under1,M70_drop_under1,
                                                             M80_drop_under1,M90_drop_under1))

write.csv(Res_resampling_modularity_drop_under1,"./R2/output/Res_resampling_modularity_drop_under1.csv")

###################drop shared bird involved links
M20_drop_shared1=multi_modular(rarefied_metanetwork1_20$drop_shared,size = 20)
M30_drop_shared1=multi_modular(rarefied_metanetwork1_30$drop_shared,size = 30)
M40_drop_shared1=multi_modular(rarefied_metanetwork1_40$drop_shared,size = 40)
M50_drop_shared1=multi_modular(rarefied_metanetwork1_50$drop_shared,size = 50)
M60_drop_shared1=multi_modular(rarefied_metanetwork1_60$drop_shared,size = 60)
M70_drop_shared1=multi_modular(rarefied_metanetwork1_70$drop_shared,size = 70)
M80_drop_shared1=multi_modular(rarefied_metanetwork1_80$drop_shared,size = 80)
M90_drop_shared1=multi_modular(rarefied_metanetwork1_90$drop_shared,size = 90)

Res_resampling_modularity_drop_shared1 <- do.call("rbind",list(M20_drop_shared1,M30_drop_shared1,M40_drop_shared1,
                                                              M50_drop_shared1,M60_drop_shared1,M70_drop_shared1,
                                                              M80_drop_shared1,M90_drop_shared1))
write.csv(Res_resampling_modularity_drop_shared1,"./R2/output/Res_resampling_modularity_drop_shared1.csv")


####################################################################################################################
######Nestedness analysis

multi_nestedness <- function(list,index="weighted NODF",size=NULL){
    Nestedness <-sapply(list,function(x)networklevel(x, index = index))
    res<-data.frame(mean(Nestedness),sd(Nestedness))
    names(res) <- c("mean","sd")
    res$samplesize <-size
    return(res)
}



##########################################################################
N20_1=multi_nestedness(rarefied_metanetwork1_20$full_webs,size = 20)
N30_1=multi_nestedness(rarefied_metanetwork1_30$full_webs,size = 30)
N40_1=multi_nestedness(rarefied_metanetwork1_40$full_webs,size = 40)
N50_1=multi_nestedness(rarefied_metanetwork1_50$full_webs,size = 50)
N60_1=multi_nestedness(rarefied_metanetwork1_60$full_webs,size = 60)
N70_1=multi_nestedness(rarefied_metanetwork1_70$full_webs,size = 70)
N80_1=multi_nestedness(rarefied_metanetwork1_80$full_webs,size = 80)
N90_1=multi_nestedness(rarefied_metanetwork1_90$full_webs,size = 90)

Res_resampling_nestedness1 <- do.call("rbind",list(N20_1,N30_1,N40_1,
                                                   N50_1,M60_1,N70_1,N80_1,N90_1))

write.csv(Res_resampling_nestedness1,"./R2/output/Res_resampling_nestedness1.csv")


###################drop on the tree bird involved links

N20_drop_on1=multi_nestedness(rarefied_metanetwork1_20$drop_on,size = 20)
N30_drop_on1=multi_nestedness(rarefied_metanetwork1_30$drop_on,size = 30)
N40_drop_on1=multi_nestedness(rarefied_metanetwork1_40$drop_on,size = 40)
N50_drop_on1=multi_nestedness(rarefied_metanetwork1_50$drop_on,size = 50)
N60_drop_on1=multi_nestedness(rarefied_metanetwork1_60$drop_on,size = 60)
N70_drop_on1=multi_nestedness(rarefied_metanetwork1_70$drop_on,size = 70)
N80_drop_on1=multi_nestedness(rarefied_metanetwork1_80$drop_on,size = 80)
N90_drop_on1=multi_nestedness(rarefied_metanetwork1_90$drop_on,size = 90)

Res_resampling_netstedness_drop_on1 <- do.call("rbind",list(N20_drop_on1,N30_drop_on1,N40_drop_on1,
                                                           N50_drop_on1,N60_drop_on1,N70_drop_on1,
                                                           N80_drop_on1,N90_drop_on1))
write.csv(Res_resampling_netstedness_drop_on1,"./R2/output/Res_resampling_nestedness_drop_on1.csv")


###################drop under the tree bird involved links
N20_drop_under1=multi_nestedness(rarefied_metanetwork1_20$drop_under,size = 20)
N30_drop_under1=multi_nestedness(rarefied_metanetwork1_30$drop_under,size = 30)
N40_drop_under1=multi_nestedness(rarefied_metanetwork1_40$drop_under,size = 40)
N50_drop_under1=multi_nestedness(rarefied_metanetwork1_50$drop_under,size = 50)
N60_drop_under1=multi_nestedness(rarefied_metanetwork1_60$drop_under,size = 60)
N70_drop_under1=multi_nestedness(rarefied_metanetwork1_70$drop_under,size = 70)
N80_drop_under1=multi_nestedness(rarefied_metanetwork1_80$drop_under,size = 80)
N90_drop_under1=multi_nestedness(rarefied_metanetwork1_90$drop_under,size = 90)

Res_resampling_nestedness_drop_under1 <- do.call("rbind",list(N20_drop_under1,N30_drop_under1,N40_drop_under1,
                                                              N50_drop_under1,N60_drop_under1,N70_drop_under1,
                                                              N80_drop_under1,N90_drop_under1))

write.csv(Res_resampling_nestedness_drop_under1,"./R2/output/Res_resampling_nestedness_drop_under1.csv")


###################drop shared bird involved links
N20_drop_shared1=multi_nestedness(rarefied_metanetwork1_20$drop_shared,size = 20)
N30_drop_shared1=multi_nestedness(rarefied_metanetwork1_30$drop_shared,size = 30)
N40_drop_shared1=multi_nestedness(rarefied_metanetwork1_40$drop_shared,size = 40)
N50_drop_shared1=multi_nestedness(rarefied_metanetwork1_50$drop_shared,size = 50)
N60_drop_shared1=multi_nestedness(rarefied_metanetwork1_60$drop_shared,size = 60)
N70_drop_shared1=multi_nestedness(rarefied_metanetwork1_70$drop_shared,size = 70)
N80_drop_shared1=multi_nestedness(rarefied_metanetwork1_80$drop_shared,size = 80)
N90_drop_shared1=multi_nestedness(rarefied_metanetwork1_90$drop_shared,size = 90)

Res_resampling_nestedness_drop_shared1 <- do.call("rbind",list(N20_drop_shared1,N30_drop_shared1,N40_drop_shared1,
                                                               N50_drop_shared1,N60_drop_shared1,N70_drop_shared1,
                                                               N80_drop_shared1,N90_drop_shared1))

write.csv(Res_resampling_nestedness_drop_shared1,"./R2/output/Res_resampling_nestedness_drop_shared1.csv")


Res_resampling_nestedness1$Network <- "Full"
Res_resampling_netstedness_drop_on1$Network <- "drop.on"
Res_resampling_nestedness_drop_under1$Network <- "drop.under"
Res_resampling_nestedness_drop_shared1$Network <- "drop.shared"

Res_resampling_nestedness1 <- rbind(Res_resampling_nestedness1,
                                    Res_resampling_netstedness_drop_on1,
                                    Res_resampling_nestedness_drop_under1,
                                    Res_resampling_nestedness_drop_shared1)



######################################################################################
######Fixed plant species richness
ls2metanetwork <- function(list,quantitative=T){#network list to a metanetwork
    patch<-names(list)
    n <- length(list)
    webs_df <- list()
    for(i in 1:n){
        tep <- (m2f(list[[i]]))
        tep$patch <- patch[i]
        webs_df[[i]] <- tep
    }
    plain_df <-do.call("rbind",webs_df)
    df <-sumbygroup(plain_df,measurevar = "Value",
                    groupvars = c("patch","Interaction_Unit"))
    
    xtabs(Value~Interaction_Unit + patch,data = df)-> web
    attr(web,"class") <- NULL
    attr(web,"call") <- NULL
    if (quantitative=="TRUE"){
        metanetwork <- web
    }else{
        web[web>0]=1
        metanetwork <- web
    }
    return(metanetwork)
    
}


sampler_web <- function(list,nplant=2,nanimal=2,quantitative=TURE,nsim=100){
    n <- length(list)
    webs <- list()
    for (i in 1:nsim) {
        list_tep <- list()
        for (j in 1:n) {
            tep1 <- list[[j]]
            success <- FALSE
            while (!success) {
                #sample process
                tep <- tep1[sample(1:nrow(tep1),nplant),sample(1:ncol(tep1),nanimal)]
                tep2 <- bipartite::empty(web=tep)
                
                plant_n2 <- nrow(tep2)
                animal_n2 <- ncol(tep2)
                
                # check for success
                success <- plant_n2==nplant & animal_n2==nanimal
            }
            tep2->list_tep[[j]]
            
        }
        webs[[i]] <-ls2metanetwork(list_tep)
    }
    
    return(webs)
}


sampler_web(Networks,nsim = 2)





rarefied_meta_fixed_species <- function(df,patch="patch",plant="Plant",
                                        animal="Animal",Value="feq",guild="guild",
                                        nplant=2,nsim=100,on="on",under="under",
                                        shared="shared",quantitative=TRUE){
    patchname <- unique(df[,patch])
    full_webs <- list()
    drop_on <- list()
    drop_under <-list()
    drop_shared <-list()
    #test <-list()
    for (i in 1:nsim) {
        dflist <-list()
        for (j in 1:length(patchname)) {
            tep <- df[df[,patch] %in% patchname[j],]
            plant_name <-unique(tep[,plant])
            #animal_name <- unique(tep[,animal])

            #success <- FALSE
            #while (!success) {
                #sample process
                tep2 <- tep[tep[,plant] %in% (sample(plant_name,nplant)),]
                #tep3 <- bipartite::empty(xtabs(feq~plant+animal,data = tep2))
                #attr(web,"class") <- NULL
                #attr(web,"call") <- NULL
                #plant_n2 <- nrow(tep3)
                #animal_n2 <- ncol(tep3)
                
                # check for success
                #success <- plant_n2==nplant
            #}
            tep2->dflist[[j]]
        }
        df_full <- do.call("rbind",dflist)
        df_drop_on <-droplevels(df_full[which(df_full[,guild]!=on),])
        df_drop_under <-droplevels(df_full[which(df_full[,guild]!=under),])
        df_drop_shared <-droplevels(df_full[which(df_full[,guild]!=shared),])
        #test[[i]] <-df_drop_shared
        full_webs[[i]] <- df2metanetwork(df=df_full,patch=patch,
                                         plant=plant,animal=animal,
                                         Value=Value,quantitative=quantitative)
        drop_on[[i]] <- df2metanetwork(df=df_drop_on,patch=patch,
                                       plant=plant,animal=animal,
                                       Value=Value,quantitative=quantitative)
        drop_under[[i]] <- df2metanetwork(df=df_drop_under,patch=patch,
                                          plant=plant,animal=animal,
                                          Value=Value,quantitative=quantitative)
        drop_shared[[i]] <- df2metanetwork(df=df_drop_shared,patch=patch,
                                           plant=plant,animal=animal,
                                           Value=Value,quantitative=quantitative)
    }
    return(list(full_webs = full_webs,drop_on=drop_on,
                drop_under= drop_under,drop_shared=drop_shared))
    
}


(rarefied_meta_fixed_species(metanetwork_df,patch="patch",plant="plant",
                            animal="animal",Value="feq",guild="guild",
                            nplant = 3,nsim=1,on="on",under="under",
                            shared="shared",quantitative=TRUE))# test


rarefied_fixed_q <- rarefied_meta_fixed_species(metanetwork_df,patch="patch",plant="plant",
                                                animal="animal",Value="feq",guild="guild",
                                                nplant = 3,nsim=1000,on="on",under="under",
                                                shared="shared",quantitative=TRUE)


rarefied_fixed_b <- rarefied_meta_fixed_species(metanetwork_df,patch="patch",plant="plant",
                                                animal="animal",Value="feq",guild="guild",
                                                nplant = 3,nsim=1000,on="on",under="under",
                                                shared="shared",quantitative=FALSE)

(rarefied_fixed_q_M <- multi_modular(rarefied_fixed_q$full_webs))
(rarefied_fixed_q_M_drop_on <- multi_modular(rarefied_fixed_q$drop_on))
(rarefied_fixed_q_M_drop_under <- multi_modular(rarefied_fixed_q$drop_under))
(rarefied_fixed_q_M_drop_shared <- multi_modular(rarefied_fixed_q$drop_shared))

(rarefied_fixed_q_N <- multi_nestedness(rarefied_fixed_q$full_webs))
(rarefied_fixed_q_N_drop_on <- multi_nestedness(rarefied_fixed_q$drop_on))
(rarefied_fixed_q_N_drop_under <- multi_nestedness(rarefied_fixed_q$drop_under))
(rarefied_fixed_q_N_drop_shared <- multi_nestedness(rarefied_fixed_q$drop_shared))

(rarefied_fixed_q_MN <- rbind(rarefied_fixed_q_M,rarefied_fixed_q_M_drop_on,
                             rarefied_fixed_q_M_drop_under,rarefied_fixed_q_M_drop_shared,
                             rarefied_fixed_q_N,rarefied_fixed_q_N_drop_on,
                             rarefied_fixed_q_N_drop_under,rarefied_fixed_q_N_drop_shared))
rarefied_fixed_q_MN$Index <- c(rep("Modularity",4),
                               rep("Nestedness",4))
rarefied_fixed_q_MN$Type <- "Quantitative"
rarefied_fixed_q_MN$Network <- c("Full","drop.on","drop.under","drop.shared",
                                 "Full","drop.on","drop.under","drop.shared")

write.csv(rarefied_fixed_q_MN,"./R2/output/rarefied_fixed_q_MN.csv")


    
#(rarefied_fixed_b_M <- multi_modular(rarefied_fixed_b))
#(rarefied_fixed_b_N <- multi_nestedness(rarefied_fixed_b,index = "NODF"))

(rarefied_fixed_b_M <- multi_modular(rarefied_fixed_b$full_webs))
(rarefied_fixed_b_M_drop_on <- multi_modular(rarefied_fixed_b$drop_on))
(rarefied_fixed_b_M_drop_under <- multi_modular(rarefied_fixed_b$drop_under))
(rarefied_fixed_b_M_drop_shared <- multi_modular(rarefied_fixed_b$drop_shared))

(rarefied_fixed_b_N <- multi_nestedness(rarefied_fixed_b$full_webs,index = "NODF"))
(rarefied_fixed_b_N_drop_on <- multi_nestedness(rarefied_fixed_b$drop_on,index = "NODF"))
(rarefied_fixed_b_N_drop_under <- multi_nestedness(rarefied_fixed_b$drop_under,index = "NODF"))
(rarefied_fixed_b_N_drop_shared <- multi_nestedness(rarefied_fixed_b$drop_shared,index = "NODF"))

(rarefied_fixed_b_MN <- rbind(rarefied_fixed_b_M,rarefied_fixed_b_M_drop_on,
                              rarefied_fixed_b_M_drop_under,rarefied_fixed_b_M_drop_shared,
                              rarefied_fixed_b_N,rarefied_fixed_b_N_drop_on,
                              rarefied_fixed_b_N_drop_under,rarefied_fixed_b_N_drop_shared))
rarefied_fixed_b_MN$Index <- c(rep("Modularity",4),
                               rep("Nestedness",4))
rarefied_fixed_b_MN$Type <- "Binary"
rarefied_fixed_b_MN$Network <- c("Full","drop.on","drop.under","drop.shared",
                                 "Full","drop.on","drop.under","drop.shared")

write.csv(rarefied_fixed_b_MN,"./R2/output/rarefied_fixed_b_MN.csv")
write.csv(rarefied_fixed_b_MN,"./R2/output/rarefied_fixed_b_MN.csv")


rarefied_fixed <-rbind(rarefied_fixed_q_MN,rarefied_fixed_b_MN)
write.csv(rarefied_fixed,"./R2/output/rarefied_fixed.csv")

rarefied_fixed <- read.csv("./R2/output/rarefied_fixed.csv")
pd <- position_dodge(0.2) # move them .05 to the left and right

fixed_species_modularity <- ggplot(subset(rarefied_fixed,Index=="Modularity"), 
       aes(x=factor(Network,levels = c("Full","drop.on","drop.under",
                                       "drop.shared")), 
           y=mean, colour=Type, group=Type, shape=Rarefication)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="grey80", 
                width=.0,size=.8, position=pd) +
  geom_point(position=pd,size=4)+xlab("")+
  ylab("Modularity")+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.y=element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x=element_text(size = 12,family = "serif"),
        legend.text =element_text(size = 12,family = "serif"),
        legend.title = element_text(size = 12,family = "serif"))

fixed_species_modularity 


fixed_species_nestedness <- ggplot(subset(rarefied_fixed,Index=="Nestedness"), 
       aes(x=factor(Network,levels = c("Full","drop.on","drop.under",
                                       "drop.shared")), 
           y=mean, colour=Type, group=Type, shape=Rarefication)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="grey80", 
                width=.0,size=.8, position=pd) +
  geom_point(position=pd,size=4)+xlab("")+
  ylab("Nestedness")+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.y=element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x=element_text(size = 12,family = "serif"),
        legend.text =element_text(size = 12,family = "serif"),
        legend.title = element_text(size = 12,family = "serif"))

fixed_species_nestedness
fixed_species_modularity+fixed_species_nestedness+plot_layout(ncol = 2)





save.image(file = "./R2/20200203analysis.Rdata")

load(file = "./R2/20200203analysis.Rdata")

save.image(file = "./R2/202002010analysis.Rdata")


for (i in 1:13) {
  #print(names(Networkdata2014)[i])
  print(dim(Networkdata2014[[i]]))
  
}

for (i in 1:13) {
  print(sum(Networkdata2014[[i]]))
  
}

for (i in 1:13) {
  print(link(Networkdata2014[[i]]))
  
}


for (i in 1:13) {
  #print(names(Networkdata2015)[i])
  print(dim(Networkdata2015[[i]]))
  
}



for (i in 1:13) {
  print(dim(Networks[[i]]))
  
}



specieslevel(q_metanetwork,index = "degree")


colSums(q_metanetwork)

for (i in 1:13) {
  print(sum(Networks[[i]]))
  
}


##################Sp. Species analysis

Plant_matrix <- read.csv("PlantSpe_matrix.csv", row.names = 1)

Bird_matrix <- read.csv("Bird_community.csv", row.names = 1)

op <- par(mfrow=c(2,1))
plotweb(sortweb(Plant_matrix, sort.order="dec"), method="normal",col.low = "#403E18",
        col.high = "#71C24C")

plotweb(sortweb(Bird_matrix, sort.order="dec"), method="normal",col.low = "#403E18",
        col.high = "#C2BA55")
par(op)

###Plant
set.seed(123)
nulls_plant <- (simulate(vegan::nullmodel(Plant_matrix, method="r1"),
                             nsim = 1000))

(con_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                            index = "connectance"))

(M_nulls_plant <- MFunwithNull_vegan(Plant_matrix,nulls_plant))


(nestedness_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                               index = "NODF")) 
(r1_plant <- rbind(con_nulls_plant,M_nulls_plant,nestedness_nulls_plant))

set.seed(123)
nulls_plant <- (simulate(vegan::nullmodel(Plant_matrix, method="r2"),
                         nsim = 1000))

(con_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                        index = "connectance"))

(M_nulls_plant <- MFunwithNull_vegan(Plant_matrix,nulls_plant))


(nestedness_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                               index = "NODF"))

(r2_plant <- rbind(con_nulls_plant,M_nulls_plant,nestedness_nulls_plant))


set.seed(123)
nulls_plant <- (simulate(vegan::nullmodel(Plant_matrix, method="quasiswap"),
                         nsim = 1000))

(con_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                        index = "connectance"))

(M_nulls_plant <- MFunwithNull_vegan(Plant_matrix,nulls_plant))


(nestedness_nulls_plant <- myFunwithNull_vegan(Plant_matrix,nulls_plant,
                                               index = "NODF"))

(quasiswap_plant <- rbind(con_nulls_plant,M_nulls_plant,nestedness_nulls_plant))


set.seed(123)
nulls_plant <- bipartite::nullmodel(Plant_matrix, method=5,
                         N = 1000)

(con_nulls_plant <- myFunwithNull(Plant_matrix,nulls_plant,
                                        index = "connectance"))

(M_nulls_plant <- MFunwithNull(Plant_matrix,nulls_plant))


(nestedness_nulls_plant <- myFunwithNull(Plant_matrix,nulls_plant,
                                               index = "NODF"))

(megen_plant <- rbind(con_nulls_plant,M_nulls_plant,nestedness_nulls_plant))

(plant_net_str <- rbind(r1_plant,r2_plant,quasiswap_plant,megen_plant))

write.csv(plant_net_str,"./R2/output/plant_net_str.csv")


###Bird
set.seed(123)
nulls_bird <- (simulate(vegan::nullmodel(Bird_matrix, method="r1"),
                         nsim = 1000))

(con_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                        index = "connectance"))

(M_nulls_bird <- MFunwithNull_vegan(Bird_matrix,nulls_bird))


(nestedness_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                               index = "NODF"))

(r1_bird <- rbind(con_nulls_bird,M_nulls_bird,nestedness_nulls_bird))

set.seed(123)
nulls_bird <- (simulate(vegan::nullmodel(Bird_matrix, method="r2"),
                         nsim = 1000))

(con_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                        index = "connectance"))

(M_nulls_bird <- MFunwithNull_vegan(Bird_matrix,nulls_bird))


(nestedness_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                               index = "NODF"))

(r2_bird <- rbind(con_nulls_bird,M_nulls_bird,nestedness_nulls_bird))


set.seed(123)
nulls_bird <- (simulate(vegan::nullmodel(Bird_matrix, method="quasiswap"),
                         nsim = 1000))

(con_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                        index = "connectance"))

(M_nulls_bird <- MFunwithNull_vegan(Bird_matrix,nulls_bird))


(nestedness_nulls_bird <- myFunwithNull_vegan(Bird_matrix,nulls_bird,
                                               index = "NODF"))

(quasiswap_bird <- rbind(con_nulls_bird,M_nulls_bird,nestedness_nulls_bird))


set.seed(123)
nulls_bird <- bipartite::nullmodel(Bird_matrix, method=5,
                                    N = 1000)

(con_nulls_bird <- myFunwithNull(Bird_matrix,nulls_bird,
                                  index = "connectance"))

(M_nulls_bird <- MFunwithNull(Bird_matrix,nulls_bird))


(nestedness_nulls_bird <- myFunwithNull(Bird_matrix,nulls_bird,
                                         index = "NODF"))

(megen_bird <- rbind(con_nulls_bird,M_nulls_bird,nestedness_nulls_bird))

(bird_net_str <- rbind(r1_bird,r2_bird,quasiswap_bird,megen_bird))

write.csv(bird_net_str,"./R2/output/bird_net_str.csv")


plant_patch.degree <- specieslevel(Plant_matrix,index = "degree")[[2]]
set.seed(123)
mod.b <- computeModules(Plant_matrix)
cz.b <- czvalues(mod.b,weighted=F,level="lower")
c.plant_patch_b <- cz.b$c
z.plant_patch_b <- cz.b$z

set.seed(123)
nestcontribution <- nestedcontribution(Plant_matrix, nsimul = 999) -> n #same as quantitative
nestcontribution_plant_patch <- nestcontribution[[2]]


##perform a PCA analysis
contribution2metanet_plant_patch.b <- data.frame(plant_patch.degree,c.plant_patch_b,
                                     z.plant_patch_b,nestcontribution_plant_patch)

str(contribution2metanet_plant_patch.b)

contribution2metanet_plant_patch.b.pca <- prcomp(as.matrix(contribution2metanet_plant_patch.b), 
                                     center = TRUE,scale. = TRUE)

summary(contribution2metanet_plant_patch.b.pca)

contribution2metanet_plant_patch.b.pca$x

write.csv(contribution2metanet.b.pca$x,"./R2/output/PCA_b.csv")

