from tkinter import TRUE
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score, recall_score, f1_score

true_lable = [0, 0, 0, 0, 1, 1, 1, 2, 2]
prediction = [0, 0, 1, 2, 1, 1, 2, 1, 2]


measure_result = classification_report(true_lable, prediction)
print('measure_result = \n', measure_result)

print("----------------------------- precision（精确率）-----------------------------")
precision_score_average_None = precision_score(true_lable, prediction, average=None)
precision_score_average_micro = precision_score(true_lable, prediction, average='micro')
precision_score_average_macro = precision_score(true_lable, prediction, average='macro')
precision_score_average_weighted = precision_score(true_lable, prediction, average='weighted')
print('precision_score_average_None = ', precision_score_average_None)
print('precision_score_average_micro = ', precision_score_average_micro)
print('precision_score_average_macro = ', precision_score_average_macro)
print('precision_score_average_weighted = ', precision_score_average_weighted)

print("\n\n----------------------------- recall（召回率）-----------------------------")
recall_score_average_None = recall_score(true_lable, prediction, average=None)
recall_score_average_micro = recall_score(true_lable, prediction, average='micro')
recall_score_average_macro = recall_score(true_lable, prediction, average='macro')
recall_score_average_weighted = recall_score(true_lable, prediction, average='weighted')
print('recall_score_average_None = ', recall_score_average_None)
print('recall_score_average_micro = ', recall_score_average_micro)
print('recall_score_average_macro = ', recall_score_average_macro)
print('recall_score_average_weighted = ', recall_score_average_weighted)

print("\n\n----------------------------- F1-value-----------------------------")
f1_score_average_None = f1_score(true_lable, prediction, average=None)
f1_score_average_micro = f1_score(true_lable, prediction, average='micro')
f1_score_average_macro = f1_score(true_lable, prediction, average='macro')
f1_score_average_weighted = f1_score(true_lable, prediction, average='weighted')
print('f1_score_average_None = ', f1_score_average_None)
print('f1_score_average_micro = ', f1_score_average_micro)
print('f1_score_average_macro = ', f1_score_average_macro)
print('f1_score_average_weighted = ', f1_score_average_weighted)


# 写循环
import pandas as pd
import numpy as np
from tkinter import TRUE
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import balanced_accuracy_score

data = pd.read_csv('/mnt/samba/temp/xjj/drug/drug_result/HDACi_chemo618/10fold-CV/SVM_tenfold_30.txt',header='infer',sep="\t")
type(data)
data[:2] ##提取第二行的数据
data['text_label'] #特定的列
data[['text_label','predict_label']] #data.loc[:,['text_label','predict_label']] 这辆是等价的
data[data.columns[2:4]] #data[['text_label','predict_label']]
#true_lable = data.iloc[i-1, 2].split(",")
#prediction = data.iloc[i-1, 3].split(",")

result = pd.DataFrame(columns=('idx',"precision_None",'precision_micro','precision_macro','precision_weighted',"recall_None","recall_micro","recall_macro","recall_weighted","f1_None","f1_micro","f1_macro","f1_weighted"))
for i in range(1,len(data)):
    #true_lable = data.iloc[i-1, 6:9].tolist()
    #prediction = data.iloc[i-1, 9:12].tolist()
    true_lable = data.iloc[i-1, 2].split(",")
    prediction = data.iloc[i-1, 3].split(",")
    #measure_result = classification_report(true_lable, prediction)
    #print('measure_result = \n', measure_result)
    precision_score_average_None = precision_score(true_lable, prediction, average=None)
    precision_score_average_micro = precision_score(true_lable, prediction, average='micro')
    precision_score_average_macro = precision_score(true_lable, prediction, average='macro')
    precision_score_average_weighted = precision_score(true_lable, prediction, average='weighted')

    recall_score_average_None = recall_score(true_lable, prediction, average=None)
    recall_score_average_micro = recall_score(true_lable, prediction, average='micro')
    recall_score_average_macro = recall_score(true_lable, prediction, average='macro')
    recall_score_average_weighted = recall_score(true_lable, prediction, average='weighted')

    f1_score_average_None = f1_score(true_lable, prediction, average=None)
    f1_score_average_micro = f1_score(true_lable, prediction, average='micro')
    f1_score_average_macro = f1_score(true_lable, prediction, average='macro')
    f1_score_average_weighted = f1_score(true_lable, prediction, average='weighted')

    a = pd.DataFrame({'idx':[i],"precision_None":[precision_score_average_None],'precision_micro':[precision_score_average_micro],'precision_macro':[precision_score_average_macro],'precision_weighted':[precision_score_average_weighted],
                       "recall_None":[recall_score_average_None],"recall_micro":[recall_score_average_micro],"recall_macro":[recall_score_average_macro],"recall_weighted":[recall_score_average_weighted],
                       "f1_None":[f1_score_average_None],"f1_micro":[f1_score_average_micro],"f1_macro":[f1_score_average_macro],"f1_weighted":[f1_score_average_weighted]})
    result = result.append(a,ignore_index=True)



result.to_csv('/mnt/samba/temp/xjj/drug/drug_result/HDACi_chemo618/classifier/try.txt', sep='\t', index=0)

