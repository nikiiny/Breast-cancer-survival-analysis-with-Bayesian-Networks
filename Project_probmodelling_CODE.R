#################################################################
############# PROBABILSTIC MODELLING PROJECT ####################
#################################################################


##########################################################################
###################### DATA PREPARATION ####################################

# read the dataset
breast_canc = read.csv('Breast Cancer METABRIC.csv')
features = colnames(breast_canc)
dim(breast_canc) #2509, 34

summary(breast_canc['Overall.Survival..Months.']) 
# there are 528 NA values

# subtitute the empty observations "" with the NA value, so that they
#can be removed easily aftwerwards
breast_canc[breast_canc==''] <- NA


# start by removing all the observations having NA value in
#correspondence to the target variable
breast_canc_data = breast_canc[!is.na(breast_canc['Overall.Survival..Months.']),]
summary(breast_canc_data['Overall.Survival..Months.'])

table(breast_canc_data['Sex']) 
#all the patients are female so we will eliminate this column afterwards

# remove all the patients which have died from causes other than breast cancer
table(breast_canc_data['Patient.s.Vital.Status']) 
breast_canc_data = breast_canc_data[breast_canc_data['Patient.s.Vital.Status']!='Died of Other Causes',]
# check
table(breast_canc_data['Patient.s.Vital.Status']) 

summary(breast_canc_data['Overall.Survival..Months.'])
# there is a single NA which we delete
breast_canc_data = breast_canc_data[complete.cases(breast_canc_data['Overall.Survival..Months.']),]

# check the distribution of the target variable
# the median value is 114 months survival (9.5 years) and the first quantile is 
#55 months survival (around 4.6 years). we create a binary variable with 1 if the
#patient has survived >120 months and 0 otherwise
summary(breast_canc_data$Overall.Survival..Months.)
breast_canc_data$Overall.Survival..Months.[breast_canc_data$Overall.Survival..Months. <60] = 0
breast_canc_data$Overall.Survival..Months.[breast_canc_data$Overall.Survival..Months. >=60] = 1
# change the name of the target
names(breast_canc_data)[24] = 'Overall.Survival.Months.TARGET'
table(breast_canc_data$Overall.Survival.Months.TARGET)

# remove:
# - the patient ID [1] 
# - the cancer type [4]
# - the choort [9]
# - the ER status measured by IHC [10]
# - the HER2 Status [14]
# - the survival status [25]
# - relapse free status in months [28]
# - the sex [30]
# - the patient's vital status [34]

breast_canc_data = breast_canc_data[,-c(1,4,9,10,14,25,28,30,34)]
dim(breast_canc_data) #1483 25


# the variable age [1] is discretized according to the quantiles in:
# - <=50
# - 50-60
# - 60-70
# - > 70
summary(breast_canc_data$Age.at.Diagnosis)
breaks = c(0,50,60,70,120)
labels = c('<=50', '50-60', '60-70', '>70')
breast_canc_data$Age.at.Diagnosis = cut(breast_canc_data$Age.at.Diagnosis, breaks, labels, right = FALSE) 

table(breast_canc_data$Age.at.Diagnosis)


# the nottingham prognostic index is transformed into an integer
#by rounding it to the unit
breast_canc_data$Nottingham.prognostic.index = round(breast_canc_data$Nottingham.prognostic.index)
unique(breast_canc_data$Nottingham.prognostic.index)


# the size of the cancer is discretised into intervals of the same
#length (10 cm) and an upper rounding from 100 cm on
breast_canc_data$Tumor.Size = round(breast_canc_data$Tumor.Size)
breaks = c(0,10,20,30,40,50,60,70,80,90,100,200)
labels = c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100','>100')
breast_canc_data$Tumor.Size = cut(unlist(breast_canc_data$Tumor.Size), breaks, labels, right = FALSE)
table(breast_canc_data$Tumor.Size)


# for the variable of mutation counts, the values from 1 to 15
#are kept as they are since there are many patients for each category,
#while the values greater than 15 are discretised in the same class
breast_canc_data$Mutation.Count[breast_canc_data$Mutation.Count > 15 & !is.na(breast_canc_data$Mutation.Count)] = '>15' 
table(breast_canc_data$Mutation.Count)


#the variable lymph nodes examined positive is kept for values up to 10 and then 
#discretized as follows:
# - 10-20
# - >20
table(breast_canc_data$Lymph.nodes.examined.positive)
breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,50)
labels = c('1','2','3','4','5','6','7','8','9','10','10-20','>20')
breast_canc_data$Lymph.nodes.examined.positive = cut(as.integer(breast_canc_data$Lymph.nodes.examined.positive), breaks, labels, right = TRUE) 
table(breast_canc_data$Lymph.nodes.examined.positive)


# check dimension and summary of variables after data cleansing
dim(breast_canc_data) #1483, 25
summary(breast_canc_data)



##########################################################################
############## 1. BAYESIAN NETWORK WITH COMPLETE DATA ######################
#############################################################################

# remove all the missing values
breast_canc_complete = breast_canc_data[complete.cases(breast_canc_data), ]

#check dimension and if there are not anymore NA values
dim(breast_canc_complete) #420, 25
breast_canc_complete[is.na(breast_canc_complete)] #no NA

# transform into factors all the categorical variables
features = colnames(breast_canc_complete)
for(i in 1:length(features)) {
  breast_canc_complete[,i] = factor(breast_canc_complete[,i])
}
#check
summary(breast_canc_complete)


##############

# download the packages
library('gRbase')
library('gRapHD')
library('qgraph')
library("Rgraphviz") 
library("RBGL")
library("gRain")
library("igraph")
library('ggm')
library('bnlearn')
library('pcalg')

library('plyr')
library('hash')

features = colnames(breast_canc_complete)


# the categories for each variable are mapped to integer numbers starting from 0,
#since some algorithms require the dataset to be in this format.
# a dictionary is built in order to map the original values to the new values and
#retrieve them later

# build a dictionary by using hash package
h <- hash() 

for(i in 1:dim(breast_canc_complete)[2]) {
  col_hash = hash()
  # for each variable retrieve its categories
  lev = levels(breast_canc_complete[,i])
  # map each category to an integer value from 0 to n-1
  legend = plyr::mapvalues(lev, from = lev, to = 0:(length(lev)-1))
  # create a mapping between of old categories and new categories
  col_hash[['levels']] = lev
  col_hash[['mapping']] = legend
  # assign the number of the column to the mapping pair in the dictionary
  h[[as.character(i)]] = col_hash
  # substitute the factors in the dataset with the new values
  levels(breast_canc_complete[,i]) = plyr::mapvalues(lev, from = lev, to = 0:(length(lev)-1))
}

#check
summary(breast_canc_complete)
h[['9']] #example: retrieve HER2 status 


########################### STRUCTURE LEARNING ####################################

# define a function to plot the network graph. the target variable is colored
#in orange
my_plot = function(mygraph, color_nodes='lightblue', n_target_var=19){
  par(mar=c(0,0,0,0))
  E(mygraph)$color='black'
  V(mygraph)$color = color_nodes
  V(mygraph)$size = 15
  V(mygraph)$color[n_target_var] = 'orange'
  par(cex=0.9)  
  plot.igraph(mygraph, edge.arrow.size=0.50)
}


############################## 1) PC-ALGORITHM
library('pcalg')
library('parallel')


# create a function to build the matrix for the white list
build_matrix = function(dim_matrix, relationships, symmetric=FALSE) {
  # build an empty matrix
  m = matrix(rep(0,dim_matrix*dim_matrix), nrow=dim_matrix)
  c=1:2
  # iterate over pairs in the vector relationships
  for(i in 1:(length(relationships)/2)) {
    # retrieve each pair of relations (from,to)
    pair = relationships[c]
    #set to 1 the corresponding cell in the matrix to indicate 
    #there is a connection between the 2
    m[pair[1],pair[2]] = 1
    if(symmetric) {
      # the pc algorithm requires a symmetric matrix
      m[pair[2],pair[1]] = 1
    }
    # go to the next pair
    c = c+2
  }
  return(m)
}

# set the connected nodes
rel_wl = c(22,19, # relapse free -> survival
           25,19, # stage -> survival
           8,25, # neoplasm histological grade -> stage
           24,25, # size -> stage
           15,25, # lymph nodes -> stage
           20,19, # PR status -> survival
           7,19, # ER status -> survival
           9,19, # HER2 status -> survival
           6,19, # Pam50...Claudin.low.subtype -> survival
           3,18, # cancer type detailed -> onco tree code
           11,22, # Hormone Therapy -> relapse free
           21,22, # Radio Therapy -> relapse free
           5,22, # Chemotherapy -> relapse free
           7,11, # ER status -> Hormone Therapy
           20,11) # PR status -> Hormone Therapy

# build the WHITE LIST for the pc algorithm
wl1 = build_matrix(dim_matrix=25,relationships=rel_wl, symmetric=TRUE)

# set the variables to build the BLACK LIST 
phisiological_char = c(1,3,4,7,8,9,10,12,13,14,15,16,20,23,24) 
phisiological_len = length(phisiological_char)

bl1 = c()
for(j in c(2,5,6,17,18,21,25)) {
  for(i in 1:length(phisiological_char)) {
    bl1 = c(bl1, j, phisiological_char[i])
  }
}

# build the black list
bl1 = build_matrix(dim_matrix=25,relationships=bl1, symmetric=TRUE)


#### MODEL 1: data-driven structure

# transform the data into a matrix to use them with the PC algorithm
M = data.matrix(breast_canc_complete)
# all the values of the categories in the matrix are increased by 1 unit, so
#we subtract that unit
M2 = M-1

# set the sufficient statistics
counter = c()
for(i in 1:dim(breast_canc_complete)[2]) {
  counter[i] = length(levels(breast_canc_complete[,i]))}
  
suffStat <- list(dm = M2, nlev = counter, adaptDF = FALSE)

# create a function to build the models by using the PC algorithm
pc_ = function(suffStat=suffStat, fixedEdges = NULL, fixedGaps = NULL, 
               return_bn_model=FALSE, arcs_to_drop=NULL) {
  # apply PC algorithm
  pc <- pcalg::pc(suffStat=suffStat, indepTest = disCItest, alpha = 0.01, 
                  labels = colnames(breast_canc_complete), verbose = TRUE,
                  fixedEdges = fixedEdges, fixedGaps = fixedGaps, u2pd='retry') 
  # transform into bn object for other operations
  g = as.bn(pc)
  # the network resulting from the PC can contain bidirected arcs, so
  #they should be transformed into directed arcs in order to fit the parameters
  #later
  if (!is.null(arcs_to_drop)) {
    for(i in 1:nrow(arcs_to_drop)) {
    g = drop.arc(g, from=arcs_to_drop[i,1], to=arcs_to_drop[i,2], debug = FALSE)
    g = set.arc(g, from=arcs_to_drop[i,2], to=arcs_to_drop[i,1], debug = FALSE)
    }
  }
  # return the bn object which can be fitted
  if (return_bn_model) { 
    return(g) 
  } 
  # return the igraph object which can be plotted
  else if (return_bn_model == FALSE) {
    g = as(amat(g), 'igraph')
    return(g)
  }
  
}

set.seed(12)
pc.1 = pc_(suffStat=suffStat)
x11()
my_plot(pc.1, color_nodes = 'green') 

#### MODEL 2: set connections a priori

set.seed(22)
pc.2 = pc_(suffStat=suffStat,fixedEdges = wl1, fixedGaps = bl1)
x11()
my_plot(pc.2, color_nodes = 'green')


############################## 2) HILL-CLIMBING
library('bnlearn')
library('gRim')

# create a function to apply the HC algorithm
hc_ = function(whitelist=NULL,blacklist=NULL,return_bn_model=FALSE,dataset=breast_canc_complete) {
  # build the structure
  hc <- hc(dataset, whitelist=whitelist, blacklist=blacklist,
           restart=10,perturb = 1)
  # transform into an igraph object
  g = as(amat(hc), 'igraph')
  # return the bn object to be fitted
  if (return_bn_model) { 
    return(hc) 
  } 
  #return the igraph object to be plotted
  else {
    return(g)
  }
}

# create the WHITE LIST for HC and MMHC
wl2 = data.frame(from = c("Relapse.Free.Status","Tumor.Stage","Neoplasm.Histologic.Grade","Tumor.Size","Lymph.nodes.examined.positive","PR.Status","ER.Status","HER2.status.measured.by.SNP6","Pam50...Claudin.low.subtype","Cancer.Type.Detailed","Hormone.Therapy","Chemotherapy","Radio.Therapy","PR.Status","ER.Status"), 
                to = c("Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Tumor.Stage","Tumor.Stage","Tumor.Stage","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Oncotree.Code","Relapse.Free.Status","Relapse.Free.Status","Relapse.Free.Status","Hormone.Therapy","Hormone.Therapy"))

# set the parameters to build the BLACK LIST
phisiological_char = c("Lymph.nodes.examined.positive","Age.at.Diagnosis","Cancer.Type.Detailed","Cellularity","ER.Status","Neoplasm.Histologic.Grade",
                       "HER2.status.measured.by.SNP6","Tumor.Other.Histologic.Subtype","Inferred.Menopausal.State","Integrative.Cluster","Primary.Tumor.Laterality","Mutation.Count",
                       "PR.Status","Tumor.Size")
phisiological_len = length(phisiological_char)
# build the black list
bl2 = data.frame(from = c( rep(c("Nottingham.prognostic.index","Oncotree.Code","X3.Gene.classifier.subtype","Tumor.Stage","Pam50...Claudin.low.subtype","Type.of.Breast.Surgery", "Chemotherapy","Hormone.Therapy","Radio.Therapy"), 
                             c(rep(phisiological_len,phisiological_len,9))), phisiological_char[-2]), 
                to = c( rep(phisiological_char,9), rep("Age.at.Diagnosis", length(phisiological_char[-2]))))


#### MODEL 1: data-driven structure
hc.1 <- hc_()
x11()
my_plot(hc.1)

#### MODEL 2: set connections a priori

hc.2 = hc_(wl2,bl2)
x11()
my_plot(hc.2) 


############################## 3) MAX-MIN HILL CLIMBING

# build a function to apply the MMHC algorithm
mmhc_ = function(whitelist=NULL,blacklist=NULL,return_bn_model=FALSE, dataset=breast_canc_complete) {
  # build the structure
  mmhc <- mmhc(dataset, whitelist=whitelist, blacklist=blacklist)
  # transform into an igraph object
  g = as(amat(mmhc), 'igraph')
  # return the BN object to be fitted
  if (return_bn_model) { 
    return(mmhc) 
  } 
  # return the igraph object to be plotted
  else {
    return(g)
  }
}


#### MODEL 1: data-driven structure
mmhc.1 <- mmhc_()
x11()
my_plot(mmhc.1, color_nodes='pink')

#### MODEL 2: set connections a priori
mmhc.2 = mmhc_(wl2,bl2)
x11()
my_plot(mmhc.2, color_nodes='pink')



##################### PARAMETERS LEARNING and VALIDATION ########################

# direct the bi-directed edges of the models built by the PC algorithm
# then return for each model the BN object for parameters' fitting
set.seed(22)
no2 = data.frame('from'=c("Inferred.Menopausal.State"),
                'to'=c("Age.at.Diagnosis"))
pc.2 = pc_(suffStat=suffStat,fixedEdges = wl1, fixedGaps = bl1,arcs_to_drop = no2,
           return_bn_model = TRUE) #51.75, 38.9

no1 = data.frame('from'=c( "Type.of.Breast.Surgery","Chemotherapy","Chemotherapy","Inferred.Menopausal.State","Overall.Survival.Months.TARGET","Hormone.Therapy"), 
                 'to'=c("Radio.Therapy","Hormone.Therapy","ER.Status","Age.at.Diagnosis","Relapse.Free.Status","ER.Status")) 
set.seed(1)
pc.1 = pc_(suffStat=suffStat,arcs_to_drop = no1,return_bn_model = TRUE) #50

hc.1 = hc_(return_bn_model=TRUE) 
hc.2 = hc_(wl2,bl2,return_bn_model=TRUE) 
mmhc.1 = mmhc_(return_bn_model=TRUE) 
mmhc.2 = mmhc_(wl2,bl2,return_bn_model=TRUE) 


# split dataset into training+validation set and test set
n = dim(breast_canc_complete)[1]
breast_canc_complete.train_val = breast_canc_complete[1:(round(n*0.8,0)-1),]
breast_canc_complete.test = breast_canc_complete[round(n*0.8,0):n,]

# Build the FUNCTIONS FOR VALIDATING THE MODELS
library(caret)
library('e1071')

# Build 2 functions to fit the parameters and make predictions

fit_predict = function(BN, train_set, val_set) {
  # fit the parameters of the model on the training set
  fitted = bn.fit(BN,train_set,method='bayes')
  # predict the target variable of the validaiton set
  y_pred=bnlearn:::predict.bn.fit(fitted,node="Overall.Survival.Months.TARGET",val_set)
  return(y_pred)
}

fit_impute_predict = function(BN, train_set, val_set) {
  # fit the parameters of the model on the training set
  fitted = bn.fit(BN,train_set,method='bayes')
  # impute the missing values in the validation set
  imputed.test = bnlearn:::impute(fitted, val_set)
  # predict the target variable of the validaiton set
  y_pred=bnlearn:::predict.bn.fit(fitted, node="Overall.Survival.Months.TARGET", imputed.test)
  return(y_pred)
}

# Build a function to return the most important statistics, which are the
#balanced accuracy, precision, recall, F1 score
avg_statistics = function(y_pred, y_true) {
  # build the confusion matrix by using the predicted and true values of the target
  conf_matr = confusionMatrix(y_pred,y_true)
  # retrieve the metrics of interest
  avg_acc = conf_matr$byClass['Balanced Accuracy']
  recall = conf_matr$byClass['Recall']
  precision = conf_matr$byClass['Precision']
  F1 = conf_matr$byClass['F1']
  #return the vector containing all of them
  stats = c(avg_acc,precision, recall, F1)
  names(stats) = c('balanced accuracy','precision', 'recall','F1')
  return(stats)
}

# Build a function to perform a CROSS-VALIDATION
kfold_CV = function(BN, train_val.set, n_folds=10, missing_data=FALSE) {
  m = dim(train_val.set)[1]
  # create intervals to divide the dataset
  v = trunc(seq(1,m,length.out=n_folds+1))
  # build a vector to store the metrics
  stats = c()
  
  for(i in 1:(length(v)-1)){
    # for each iteration leave out one of the k folds as the validation
    #set and use the rest as the training set.
    val = v[i]:v[i+1]
    val.set = train_val.set[val,]
    train = setdiff(1:m,val)
    train.set = train_val.set[train,]
    # if the dataset contains missing values, impute the data of the validation
    #set for predicting the target variable after the fitting.
    if(missing_data==TRUE) {
      y_pred = fit_impute_predict(BN,train.set, val.set)
    }
    # fit the data and predict the classification of the validation set
    else {
      y_pred = fit_predict(BN,train.set, val.set)
    }
    # return the average statistics and store the into the vector stats
    stat = avg_statistics(y_pred,val.set$Overall.Survival.Months.TARGET)
    stats = c(stats, stat)
  }
  # sum each different metric and divide by the number of iterations
  #of the cross-validation
  avg_balaccuracy = sum(stats[seq(1,length(stats),by=4)])/n_folds
  avg_precision = sum(stats[seq(2,length(stats),by=4)])/n_folds
  avg_recall = sum(stats[seq(3,length(stats),by=4)])/n_folds
  avg_F1 = sum(stats[seq(4,length(stats),by=4)])/n_folds
  # return the vector containing the average metrics for the CV
  avg_stats = c(avg_balaccuracy,avg_precision,avg_recall,avg_F1)
  names(avg_stats) = c('avg_balaccuracy','avg_precision','avg_recall','avg_F1')
  
  return(avg_stats)
}


# build a dataframe to store the result of inner CV for all the models and 
#compare their performance
set.seed(100)
cv_results = data.frame()
cv_results = rbind(cv_results, kfold_CV(pc.1, breast_canc_complete.train_val))
cv_results = rbind(cv_results, kfold_CV(pc.2, breast_canc_complete.train_val))
cv_results = rbind(cv_results, kfold_CV(hc.1, breast_canc_complete.train_val))
cv_results = rbind(cv_results, kfold_CV(hc.2, breast_canc_complete.train_val))
cv_results = rbind(cv_results, kfold_CV(mmhc.1, breast_canc_complete.train_val))
cv_results = rbind(cv_results, kfold_CV(mmhc.2, breast_canc_complete.train_val))
colnames(cv_results) = c('avg_balaccuracy','avg_precision','avg_recall','avg_F1')
rownames(cv_results) = c('pc.1', 'pc.2', 'hc.1', 'hc.2', 'mmhc.1', 'mmhc.2')
cv_results

# the best model according to the CV, excluding pc.1 which contains too few
#variables connected to the target, is mmhc.2.

#####################

# Proceed to evaluate the effective performance on the test set through an outer
#CV. First, the variables corresponding to nodes not connected to the main network
#are eliminated

new_data_complete = breast_canc_complete[,-c(3,4,10,13,14,16,18,23)] 
# build new white list
new_wl2 = data.frame(from = c("Relapse.Free.Status" ,"Tumor.Stage","Neoplasm.Histologic.Grade","Tumor.Size","Lymph.nodes.examined.positive","PR.Status","ER.Status","HER2.status.measured.by.SNP6","Pam50...Claudin.low.subtype","Hormone.Therapy","Chemotherapy","Radio.Therapy","PR.Status","ER.Status"), 
                 to = c("Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Tumor.Stage","Tumor.Stage","Tumor.Stage","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Overall.Survival.Months.TARGET","Relapse.Free.Status","Relapse.Free.Status","Relapse.Free.Status","Hormone.Therapy","Hormone.Therapy"))
phisiological_char_2 = c("Lymph.nodes.examined.positive","Age.at.Diagnosis","ER.Status","Neoplasm.Histologic.Grade",
                       "HER2.status.measured.by.SNP6","Inferred.Menopausal.State","PR.Status","Tumor.Size")
phisiological_len_2 = length(phisiological_char_2)
# build new black list
new_bl2 = data.frame(from = c( rep(c("Nottingham.prognostic.index","Tumor.Stage","Pam50...Claudin.low.subtype","Type.of.Breast.Surgery", "Chemotherapy","Hormone.Therapy","Radio.Therapy"), 
                               c(rep(phisiological_len_2,phisiological_len_2,7))), phisiological_char_2[-2]), 
                 to = c( rep(phisiological_char_2,7), rep("Age.at.Diagnosis", length(phisiological_char_2[-2]))))

# refit the MMHC model and plot the new graph
mmhc.2 = mmhc_(new_wl2,new_bl2,dataset=new_data_complete)
x11()
my_plot(mmhc.2, color_nodes='pink',n_target_var=12)

# return the BN model for the MMHC
mmhc.2 = mmhc_(new_wl2,new_bl2,dataset=new_data_complete, return_bn_model = TRUE)
# perform an outer cross-validation by using the training+validation set and
#the test set to evaluate the effective performance
set.seed(458)
performance = kfold_CV(mmhc.2, new_data_complete) 
performance

################

# Perform again the CV for different classifications threshold in order to see
#whether the metrics improve for specific thresholds

# Build a function to make classifications based on the threshold value
change_threshold_prediction = function(BN, train.set, test.set, threshold,missing.values=FALSE){
  # fit the model's parameters
  fitted = bn.fit(BN,train.set,method='bayes')
  # if there are missing values impute the missing values in the test set
  if(missing.values) {
    test.set = bnlearn:::impute(fitted, test.set)
  }
  # make predictions and return the probabilities of the target variable to
  #belong to each class
  y_pred = bnlearn:::predict.bn.fit(fitted, node="Overall.Survival.Months.TARGET",test.set, prob=TRUE)
  # retrieve the probabilities
  new_ypred = attributes(y_pred)$prob[2,]
  # classify an observation according to the value of the threshold
  new_ypred[new_ypred > threshold] =1
  new_ypred[new_ypred <= threshold] =0
  new_ypred = factor(new_ypred)
  # return the statistics of interest 
  return(avg_statistics(new_ypred, test.set$Overall.Survival.Months.TARGET))
}

# Build a new cross-validation function to obtain average metrics for different
#threshold values
kfold_CV_wthreshold = function(BN, train_val.set, n_folds=10, missing.values=FALSE) {
  m = dim(train_val.set)[1]
  # create intervals to divide the dataset
  v = trunc(seq(1,m,length.out=n_folds+1))
  # build an empty matrix to store the values of the metrics
  stats = matrix(rep(0,12*4), ncol=4)
  
  for(i in 1:(length(v)-1)){
    # use as test set one of the k folds and the rest as training set for
    #each iteration
    val = v[i]:v[i+1]
    val.set = train_val.set[val,]
    train = setdiff(1:m,val)
    train.set = train_val.set[train,]
    
    # build a dataframe to store the metrics for each threshold level
    thresh = data.frame()
    for(i in seq(35,90,by=5)/100) {
      thresh = rbind(thresh,change_threshold_prediction(BN,train.set, val.set,
                                                        threshold = i,missing.values=missing.values))
    }
    # sum the metrics values to the previous values
    stats = stats + as.matrix(thresh)
  }
  stats = as.data.frame(stats)
  rownames(stats) = seq(35,90,by=5)/100
  colnames(stats) = c('balanced accuracy','precision','recall','F1')
  
  # divide the metrics by the number of iterations of the CV to obtain
  #the average value
  return(stats/n_folds)
}

set.seed(1)
kfold_CV_wthreshold(mmhc.2, new_data_complete, n_folds=10, missing.values=FALSE)

#####################

library(pROC)
# Build a function to refit the model and then make predictions on the test set,
#then plot the ROC curve
plot_roc = function(BN, train_test.set, missing_variables=FALSE) {
  # create training and test set
  n = dim(train_test.set)[1]
  train_val = train_test.set[1:(round(n*0.8,0)-1),]
  test = train_test.set[round(n*0.8,0):n,]
  # if there are missing values impute the test set before making predictions
  if(missing_variables){
    y_pred=fit_impute_predict(BN,train_val,test)
  }
  # if there are no missing values just fit and make predictions
  else{
  y_pred=fit_predict(BN,train_val,test)
  }
  # plot the ROC curve and its AUC
  pROC_obj <- roc(as.numeric(y_pred), as.numeric(test$Overall.Survival.Months.TARGET),
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)
    # plot the confidence bands for the ROC
    sens.ci <- ci.se(pROC_obj)
    plot(sens.ci, type="shape", col="lightblue")
    plot(sens.ci, type="bars")
}


# plot the ROC and AUC for the complete version of the dataset
set.seed(2)
x11()
plot_roc(mmhc.2, new_data_complete)

# plot the ROC and AUC for the imputed version of the dataset
# NB: before doing it run the preprocessing for the dataset containing missing
#values in the next section in order to have the correct format of the variables 
set.seed(15)
x11()
plot_roc(mmhc.2, breast_canc_data[,-c(3,4,10,13,14,16,18,23)], missing_variables = TRUE)



##########################################################################
############## 2. BAYESIAN NETWORK WITH MISSING VARIABLES #################
#############################################################################

###############
# check the number of downloads of the package bnstruct
library(httr)
# downloads in the last 3 years
url = 'http://cranlogs.r-pkg.org/downloads/total/2018-01-01:2021-03-28/bnstruct'
x=GET(url) 
d=content(x)
d
# downloads in the last year
url = 'http://cranlogs.r-pkg.org/downloads/total/2020-01-01:2021-03-28/bnstruct'
x=GET(url) 
d=content(x)
d
# in the last 3 years there has been over 38661 downloads and in the last year 23503
# this means that probably this package is reliable, so it will be used
#in order to perform inference on the bayesian model by using missing values

###############

library('bnstruct')
library('hash')

# transform into factors all the categorical variables
features = colnames(breast_canc_data)
for(i in 1:length(features)) {
  breast_canc_data[,i] = factor(breast_canc_data[,i])
}

#############

# Plot the histograms of the variables divided by survival
library('ggplot2')
library('gridExtra')

plot_hist = function(variable, lab_x){
  p = ggplot(breast_canc_data, aes(x=variable,
                                   group=Overall.Survival.Months.TARGET,
                                   fill=Overall.Survival.Months.TARGET))+
    geom_bar(width = 0.5,position='fill')+ 
    scale_fill_discrete(name = "Survival (>5 years)")+
    labs(x=lab_x, title=' ')
  
  return(p)
}

p1 = plot_hist(breast_canc_data$Tumor.Stage,'Tumor Stage')
p2 = plot_hist(breast_canc_data$Relapse.Free.Status,'Relapse Free Status')
p3 = plot_hist(breast_canc_data$Pam50...Claudin.low.subtype,'Pam 50-Claudin-low subtype')
p4 = plot_hist(breast_canc_data$ER.Status,'ER Status')
p5 = plot_hist(breast_canc_data$PR.Status,'PR Status')
p6 = plot_hist(breast_canc_data$HER2.status.measured.by.SNP6,'HER2 Status')
x11()
grid.arrange(p1,p2,p3,p4,p5,p6,
             nrow=2)

###############

# map the categories of each variable to integer numbers starting
#from 0 since the algorithms need a matrix as input
h2 <- hash() 

for(i in 1:dim(breast_canc_data)[2]) {
  col_hash = hash()
  lev = levels(breast_canc_data[,i])
  legend = plyr::mapvalues(lev, from = lev, to = 0:(length(lev)-1))
  col_hash[['levels']] = lev
  col_hash[['mapping']] = legend
  h2[[as.character(i)]] = col_hash
  levels(breast_canc_data[,i]) = plyr::mapvalues(lev, from = lev, to = 0:(length(lev)-1))
}
# check
summary(breast_canc_data)
h2[['1']] #example: age at diagnosis

dim(breast_canc_data) #1483 25

########################### STRUCTURE LEARNING ####################################

# Create a function in order to transform the dataframe into the
#proper format
build_data_matrix_forBN = function(dataset) {
  # create a matrix from the dataframe
  Matrix = data.matrix(dataset)-1
  # retrieve the number of categories for each variable
  counter = c()
  for(i in 1:dim(dataset)[2]) {
    counter[i] = length(levels(dataset[,i]))
  }
  # transform the matrix into the BNDataset format used for the bnstruct
  #package
  data = BNDataset(data=Matrix, discreteness=rep('c', dim(dataset)[2]),
          variables = colnames(dataset),
          num.variables=1)
  # set the number of categories for each level
  data@node.sizes = counter
  return(data)
}

data = build_data_matrix_forBN(breast_canc_data)

# create the WHITE LIST
rel_wl = c(22,19, # relapse free -> survival
           25,19, # stage -> survival
           8,25, # neoplasm histological grade -> stage
           24,25, # size -> stage
           15,25, # lymph nodes -> stage
           20,19, # PR status -> survival
           7,19, # ER status -> survival
           9,19, # HER2 status -> survival
           6,19, # Pam50...Claudin.low.subtype -> survival
           3,18, # cancer type detailed -> onco tree code
           11,22, # Hormone Therapy -> relapse free
           21,22, # Radio Therapy -> relapse free
           5,22, # Chemotherapy -> relapse free
           7,11, # ER status -> Hormone Therapy
           20,11) # PR status -> Hormone Therapy

# in this case the matrix is not symmetric
wl3 = build_matrix(dim_matrix=25,relationships=rel_wl, symmetric=FALSE)

# create the matrix for the BLACK LIST.
# the layers are:
# 1: indexes or therapies
# 2: physiological states
# 3: age at diagnosis
# 4: target variable and relapse free status
# allowed connections are: 2->1, 3->1, 3->2
layers = c(3,1,2,2,
           1,1,2,2,
           2,2,1,2,
           2,2,2,2,
           1,1,4,2,
           1,4,1,2,1)
order = matrix(rep(0,4*4), nrow = 4)
order[2,1]=1
order[3,1]=1
order[1,4]=1
order[2,4]=1
order[3,4]=1
order[4,4]=1
order[3,2]=1
order[2,2]=1
order[1,1]=1


############################## 1) MMHC

#### MODEL 1: data-driven structure

# Create a function to apply the different algorithms
missing_bn = function(data,algo, scoring.func, mandatory.edges=NULL,
                      layering=NULL,layer.struct=NULL,return_bn_model=FALSE) {
  # learn the network structure
  net = learn.network(data, algo = algo, scoring.func=scoring.func,
                      mandatory.edges=mandatory.edges,layering=layering,
                      layer.struct=layer.struct)
  colnames(net@dag) = net@variables
  # extract the adjacency matrix
  bnet = net@dag
  # convert into igraph object to be plotted
  g = graph_from_adjacency_matrix(bnet)
  if (return_bn_model) { 
    # convert into BN object to be fitted
    return( as.bn(as(graph_from_adjacency_matrix(net@dag),'graphNEL')) ) 
  } 
  else {
    return(g)
  }
}

mmhc_AIC1 = missing_bn(data,algo='mmhc', scoring.func = 'AIC')
x11()
my_plot(mmhc_AIC1, color_nodes='pink')

mmhc_BIC1 = missing_bn(data,algo='mmhc', scoring.func = 'BIC')
x11()
my_plot(mmhc_BIC1, color_nodes='pink')


#### MODEL 2: set connections a priori
mmhc_BIC2 = missing_bn(data,algo='mmhc', scoring.func = 'BIC', 
                      mandatory.edges=wl3, layering=layers,
                      layer.struct = order)
x11()
my_plot(mmhc_BIC2, color_nodes='pink')

mmhc_AIC2 = missing_bn(data,algo='mmhc', scoring.func = 'AIC', 
                       mandatory.edges=wl3, layering=layers,
                       layer.struct = order)
x11()
my_plot(mmhc_AIC2, color_nodes='pink')


############################## 2) HC

#### MODEL 1: data-driven structure
hc_AIC1 = missing_bn(data,algo='hc', scoring.func = 'AIC') 
x11()
my_plot(hc_AIC1)

hc_BIC1 = missing_bn(data,algo='hc', scoring.func = 'BIC') 
x11()
my_plot(hc_BIC1)

#### MODEL 2: set connections a priori
hc_AIC2 = missing_bn(data,algo='hc', scoring.func = 'AIC', 
                      mandatory.edges=wl3, layering=layers,
                      layer.struct = order)
x11()
my_plot(hc_AIC2)

hc_BIC2 = missing_bn(data,algo='hc', scoring.func = 'BIC', 
                    mandatory.edges=wl3, layering=layers,
                    layer.struct = order)
x11()
my_plot(hc_BIC2)

##################### PARAMETERS LEARNING and VALIDATION ########################

# return the models as BN objects for parameters' fitting
mmhc_AIC1 = missing_bn(data,algo='mmhc', scoring.func = 'AIC',return_bn_model = TRUE) 
mmhc_BIC1 = missing_bn(data,algo='mmhc', scoring.func = 'BIC',return_bn_model = TRUE) 
mmhc_BIC2 = missing_bn(data,algo='mmhc', scoring.func = 'BIC', 
                       mandatory.edges=wl3, layering=layers,
                       layer.struct = order,return_bn_model = TRUE)
mmhc_AIC2 = missing_bn(data,algo='mmhc', scoring.func = 'AIC', 
                       mandatory.edges=wl3, layering=layers,
                       layer.struct = order,return_bn_model = TRUE) 
hc_AIC1 = missing_bn(data,algo='hc', scoring.func = 'AIC',return_bn_model = TRUE) 
hc_BIC1 = missing_bn(data,algo='hc', scoring.func = 'BIC',return_bn_model = TRUE) 
hc_AIC2 = missing_bn(data,algo='hc', scoring.func = 'AIC', 
                     mandatory.edges=wl3, layering=layers,
                     layer.struct = order,return_bn_model = TRUE) 
hc_BIC2 = missing_bn(data,algo='hc', scoring.func = 'BIC', 
                     mandatory.edges=wl3, layering=layers,
                     layer.struct = order,return_bn_model = TRUE) 

# divide the data into training+validation set and test set
n = dim(breast_canc_data)[1]
breast_canc_mis.train_val = breast_canc_data[1:(round(n*0.8,0)-1),]
breast_canc_mis.test = breast_canc_data[round(n*0.8,0):n,]

# perform inner CV for all the models to learn the parameters and make predictions
#on the validation set, and then compare their performance
set.seed(190)
cv_results = data.frame()
cv_results = rbind(cv_results, kfold_CV(hc_AIC1, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(hc_BIC1, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(hc_AIC2, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(hc_BIC2, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(mmhc_AIC1, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(mmhc_BIC1, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(mmhc_AIC2, breast_canc_mis.train_val,missing_data=TRUE))
cv_results = rbind(cv_results, kfold_CV(mmhc_BIC2, breast_canc_mis.train_val,missing_data=TRUE))
colnames(cv_results) = c('avg_balaccuracy','avg_precision','avg_recall','avg_F1')
rownames(cv_results) = c('hc_AIC1', 'hc_BIC1', 'hc_AIC2', 'hc_BIC2', 
                         'mmhc_AIC1', 'mmhc_BIC1','mmhc_AIC2','mmhc_BIC2')
cv_results

# the chosen model is hc_AIC2

# the effective performance is tested by refitting the parameters on the 
#training+validation set and making prediction on the test set.
# The variables not included in the network are eliminated from the dataset
new_data_missing = breast_canc_data[,-c(1,2,4,10,12)] 

# create the new WHITE LIST
rel_wl = c(17,14, #"Relapse.Free.Status" -> target
           20,14, #"Tumor.Stage"-> target
           5,20, #"Neoplasm.Histologic.Grade"->"Tumor.Stage"
           19,20, #"Tumor.Size"->"Tumor.Stage"
           10,20, #"Lymph.nodes.examined.positive"->"Tumor.Stage"
           15,14, #"PR.Status"->target
           4,14, #"ER.Status"->target
           6,14, #"HER2.Status"->target
           3,14, #"Pam50...Claudin.low.subtype"->target
           1,13, #"Cancer.Type.Detailed" ->"Oncotree.Code"
           7,17, #"Hormone.Therapy" -> "Relapse.Free.Status"
           16,17, #"Radio.Therapy"-> "Relapse.Free.Status"
           2,17, # "Chemotherapy"-> "Relapse.Free.Status"
           4,7, #"ER.Status" -> "Hormone.Therapy"
           15,7) #"PR.Status" -> "Hormone.Therapy"

new_wl3 = build_matrix(dim_matrix=20,relationships=rel_wl, symmetric=FALSE)

# create the new BLACK LIST
layers = c(2,1,1, 
           2,2,2,1, 
           2,2,2,2,
           1,1,3,2,
           1,3,1,2,1)
order = matrix(rep(0,3*3), nrow = 3)
order[2,1]=1
order[1,3]=1
order[2,3]=1
order[2,2]=1
order[1,1]=1

# prepare the new data
data = build_data_matrix_forBN(new_data_missing)

# fit the model on the new data and plot it
set.seed(15)
hc_AIC2 = missing_bn(data,algo='hc', scoring.func = 'AIC', 
                     mandatory.edges=new_wl3, layering=layers,
                     layer.struct = order) 
x11()
my_plot(hc_AIC2,n_target_var=14)


# fit the model and perform an outer cross validation to assess its effective 
#performance on the imputed dataset
set.seed(15)
hc_AIC2 = missing_bn(data,algo='hc', scoring.func = 'AIC', 
                       mandatory.edges=new_wl3, layering=layers,
                       layer.struct = order,return_bn_model = TRUE) 

set.seed(11)
performance = kfold_CV(hc_AIC2, new_data_missing, missing_data=TRUE) 
performance 

# retrieve the performance of the model also for the complete dataset
set.seed(13)
performance = kfold_CV(hc_AIC2, breast_canc_complete[,-c(2,4,10,12,1)]) 
performance 

###############

# plot the ROC and AUC curve for the imputed dataset
library(pROC)
set.seed(10) 
x11()
plot_roc(hc_AIC2, new_data_missing,missing_variables = TRUE)

# plot the ROC and AUC curve for the complete dataset
set.seed(26) 
x11()
plot_roc(hc_AIC2, breast_canc_complete[,-c(2,4,10,12,1)])


# perform a CV for all the different values of thresholds and assess
#their average performance
# on the complete dataset
set.seed(2)
kfold_CV_wthreshold(hc_AIC2, breast_canc_complete[,-c(2,4,10,12,1)], n_folds=10, missing.values=FALSE)

# on the imputed dataset
set.seed(3)
kfold_CV_wthreshold(hc_AIC2, new_data_missing, n_folds=10, missing.values=TRUE)


################################# INFERENCE #######################################

# create a function to print the result of queries after having 
#set the evidence
print_evidence = function(fitted_model, nodes, states) {
  # transform the object into a gRain object for making queries
  fitted_grain = as.grain(fitted_model)
  print('Marginal prob TARGET:')
  # print the marginal probabilities of the target variable
  print(querygrain(fitted_grain, nodes=c("Overall.Survival.Months.TARGET"),type='marginal'))
  # set evidence
  ev.1 = setFinding(fitted_grain, nodes=nodes, 
                    states=states)
  print('Conditional prob TARGET:')
  # print the new conditional probability of the target variable given
  #the evidence
  print(querygrain(ev.1, nodes=c("Overall.Survival.Months.TARGET"), type='conditional'))
  
}

# Create a function to return a dataframe containing the conditional
#probabilities of the target given all the possible values assumed by
#another variable
df_evidence = function(fitted_model, nodes, n_node) {
  # transform into a gRain object
  fitted_grain = as.grain(fitted_model)
  cond.prob = data.frame()
  # set the evidence for all the possible values of the considered variable
  for(i in levels(breast_canc_complete[,nodes])){
    ev.1 = setFinding(fitted_grain, nodes=nodes, 
                      states=i)
    # retrieve the new conditional probability of the target
    cp=querygrain(ev.1, nodes="Overall.Survival.Months.TARGET", type='conditional')
    # store the result in the dataset
    cond.prob = rbind(cond.prob,cp)}
  colnames(cond.prob) = c('survival<=5','survival>5')
  rownames(cond.prob)=h[[n_node]]$levels
  # return the dataframe
  return(cond.prob)
}


###### MODEL 1:
fitted_mmhc.2 = bn.fit(mmhc.2,new_data_complete,method='bayes')
fitted_grain = as.grain(fitted_mmhc.2)
# return the marginal probability of the target variable
querygrain(fitted_grain,nodes="Overall.Survival.Months.TARGET", type='marginal')

# HER2 status
df_evidence(fitted_mmhc.2,nodes="HER2.status.measured.by.SNP6",n_node='9')

# Tumour stage
df_evidence(fitted_mmhc.2,nodes="Tumor.Stage",n_node='25')

# relapse free status not recurred
df_evidence(fitted_mmhc.2,nodes="Relapse.Free.Status",n_node='22')

# ER status neg
df_evidence(fitted_mmhc.2,nodes="ER.Status",n_node='7')

# PR status neg
df_evidence(fitted_mmhc.2,nodes="PR.Status",n_node='20')

# PAM 50 
df_evidence(fitted_mmhc.2,nodes="Pam50...Claudin.low.subtype",n_node='6')


###### MODEL 2:
fitted_hc_AIC2 = bn.fit(hc_AIC2,new_data_missing,method='bayes')
fitted_grain = as.grain(fitted_hc_AIC2)
# return the marginal probability of the target variable
querygrain(fitted_grain,nodes="Overall.Survival.Months.TARGET", type='marginal')

# HER2 status
df_evidence(fitted_hc_AIC2,nodes="HER2.status.measured.by.SNP6",n_node='9')

# Tumour stage
df_evidence(fitted_hc_AIC2,nodes="Tumor.Stage",n_node='25')

# relapse free status not recurred
df_evidence(fitted_hc_AIC2,nodes="Relapse.Free.Status",n_node='22')

# ER status neg
df_evidence(fitted_hc_AIC2,nodes="ER.Status",n_node='7')

# PR status neg
df_evidence(fitted_hc_AIC2,nodes="PR.Status",n_node='20')

# PAM 50 
df_evidence(fitted_hc_AIC2,nodes="Pam50...Claudin.low.subtype",n_node='6')
