#载入包以及读入数据
library(caret)
data(mdrr)
dim(mdrrDescr)
#删除方差为0的变量
zerovar=nearZeroVar(mdrrDescr)
newdata1=mdrrDescr[,-zerovar]
#可以发现现在只剩下297个变量了
dim(newdata1)
#首先删除强相关的变量
descrCorr = cor(newdata1)
highCorr = findCorrelation(descrCorr, 0.90)
newdata2 = newdata1[, -highCorr]

#随后解决多重共线性，本例中不存在多重共线性问题
comboInfo = findLinearCombos(newdata2)
newdata2=newdata2[, -length(comboInfo$remove)]
#我们还需要将数据进行标准化并补足缺失值，这时可以用preProcess命令
Process = preProcess(newdata2)
newdata3 = predict(Process, newdata2)

ctrl= rfeControl(functions = rfFuncs, method = "repeatedcv",verbose = FALSE, returnResamp = "final")
#functions是确定用什么样的模型进行自变量排序，本例选择的模型是随机森林即rfFuncs，可以选择的
#还有lmFuncs（线性回归），nbFuncs（朴素贝叶斯），treebagFuncs（装袋决策树），caretFuncs
#(自定义的训练模型）。svmRadial
#method是确定用什么样的抽样方法，本例使用cv即交叉检验, 还有提升boot以及留一交叉检验LOOCV

#最后使用rfe命令进行特征选择
Profile = rfe(newdata3, mdrrClass,  rfeControl = ctrl)
print(Profile)
plot(Profile)

#预处理完成后，我们就需要使用train函数选择合适的机器学习算法进行训练以及使用predict函数对结果进行预测
#首先按照比例划分训练集与测试集
newdata4=newdata3[,Profile$optVariables]
inTrain = createDataPartition(mdrrClass, p = 3/4, list = FALSE)
trainx = newdata4[inTrain,]
testx = newdata4[-inTrain,]
trainy = mdrrClass[inTrain]
testy = mdrrClass[-inTrain]

#作图查看前6个变量的分布情况
featurePlot(trainx[,1:6],trainy,plot='box')

#######以下importance
#在正式训练前，首先需要使用trainControl函数定义模型训练参数，method确定多次交叉检验的抽样方法，number确定了划分的重数， repeats确定了反复次数。
#method 重抽样方法：“boot”, “boot632”, “optimism_boot”, “boot_all”, “cv”, “repeatedcv”, “LOOCV”, “LGOCV” (for repeated training/test splits), “none” (only fits one model to the entire training set), “oob” (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), timeslice, “adaptive_cv”, “adaptive_boot” or “adaptive_LGOCV”
#number folds的数量或重抽样的迭代次数
#repeats 仅作用于k折交叉验证：代表要计算的完整折叠集的数量
#p 仅作用于分组交叉验证：代表训练集的百分比
#search Either “grid” or “random”，表示如何确定调整参数网格
fitControl = trainControl(method = "repeatedcv", number = 10, repeats = 3,returnResamp = "all")

#使用train训练模型，本例中使用的时gbm算法，我们可以对一些参数进行手动调优，包括interaction.depth,n.trees,shrinkage，n.minobsinnode等参数，也可以使用默认参数
## rf, gbm, C5.0是决策树模型中的算法, nnet前馈反向传播神经网络算法,svmLinear,svmRadial SVM的非线性核预测
## knn, rpart, lm, 
gbmFit = train(trainx,trainy,method = "svmRadial",trControl = fitControl,verbose = FALSE,preProc = c("center", "scale"),weights=0.5)

#随后预测模型准确度
cell_labels.class <- as.character(predict.train(gbmFit, 
                                                newdata = testx, 
                                                type = "raw"))
cell_labels.prob <- predict.train(gbmFit, 
                                  newdata = testx, 
                                  type = "prob")



predictions<-predict(gbmFit, newdata = testx)
pred<-confusionMatrix(predictions,testy)
### 画ROC曲线
library(pROC)
roc=roc(testy, factor (predictions, ordered = T))
roc
Specificity=roc$specificities
Sensitivity=roc$sensitivities
library(ggplot2)
p <- ggplot (data = NULL, mapping= aes (x= 1-Specificity, y= Sensitivity))
p+geom_line(colour = 'red', size = 1) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))+ 
  geom_abline(intercept = 0,slope = 1)+
  annotate('text', x = 0.5, y = 0.25, label=paste('AUC=', round (roc$auc, 2)))+ 
  labs (x = '1-Specificity',y = 'Sensitivity', title = 'ROC curve') +
  theme (plot.title = element_text (hjust = 0.5, face = 'bold', colour = 'brown'))


rffit <- train(x = train_sce[,good.marker], 
               y = factor(train_sce$gating),
               method = "rf", 
               ntree = 500,
               tuneLength = 10,
               trControl = fitControl,
               allowParallel = TRUE)
cell_labels.class <- as.character(predict.train(rffit, 
                                                newdata = unlabel_data[,good.marker], 
                                                type = "raw"))
cell_labels.prob <- predict.train(rffit, 
                                  newdata = unlabel_data[,good.marker], 
                                  type = "prob")


