import tkinter.filedialog
import numpy as np
import copy

'''
    This function is to calculate the frequent item 
    satisfying support with of length 1
'''
def aprioriAlgo_fq1(data, sup, total):
    fqItemSet =[]
    countSet = []

    #Calculating the Count for the 'Up' and 'Down' Genes
    for i in range (0,data.shape[1]-1):
        countUp=0
        countDown =0
        for j in range(0,data.shape[0]):
            if (data[j][i] == 'Up'):
                countUp+=1
            elif(data[j][i] == 'Down'):
                countDown+=1
        if ((countUp/total) *100 >= sup):
            countSet.append(countUp)
            fqItemSet.append([labelItemSetUp[i]])
        if ((countDown/total) *100 >= sup):
            fqItemSet.append([labelItemSetDown[i]])
            countSet.append(countDown)
    #Checking whether the last column diseases are more than the support
    diseases=set(data[:data.shape[0], data.shape[1]-1])
    for d in diseases:
        countD =0
        for j in range(0, data.shape[0]):
            if (d == data[j][data.shape[1]-1]):
                countD += 1
        if ((countD/total)*100) >= sup:
            fqItemSet.append([str('G101_'+d)])
            countSet.append(countD)
    return fqItemSet, countSet

'''
    This function is to calculate the frequent item 
    satisfying support with of length 2
'''
def aprioriAlgo_fq2(itemset, data, sup, total):
    res = []

    for i in range(len(itemset) - 1):
        a=itemset[i][0].split('_')
        for j in range(i + 1, len(itemset)):
            b=itemset[j][0].split('_')
            if a[0] != b[0]:
                res.append([itemset[i][0], itemset[j][0]])

    fqItemSet = []
    countSet = []
    # Calculating the Count for the 'Up' and 'Down' Genes
    for item in res:
        item1 = item[0].split('_')
        item2 = item[1].split('_')
        count = 0
        for i in range(0, data.shape[0]):
            if (data[i][int(item1[0][1:]) -1] == str(item1[1]) and
                    data[i][int(item2[0][1:]) -1] == str(item2[1])):
                count += 1
        if ((count/total)*100) >= sup:
            fqItemSet.append(item)
            countSet.append(count)
    return fqItemSet, countSet

'''
    This function is to calculate the frequent item 
    satisfying support with of length >2 using Apriori Algorithm
'''
def aprioriAlgo_fqAll(itemset, data, sup, total):
    fqItemset =[]
    countSet = []
    y = len(itemset[0])
    for i in range(0, len(itemset) - 1):
        for j in range(i + 1, len(itemset)):
            if itemset[i][0:y - 1] == itemset[j][0:y - 1] and itemset[i][y - 1] != itemset[j][y - 1]:
                val = itemset[i][0:y - 1] + [itemset[i][y - 1]] + [itemset[j][y - 1]]
                #val = list(set(itemset[i] + itemset[j]))
                # Calculating the Count for the 'Up' and 'Down' Genes
                itemList=[]
                for item in val:
                    itemList.append(item.split('_'))
                count = 0
                for k in range(0, data.shape[0]):
                    flag = True
                    for element in itemList:
                        if data[k][int(element[0][1:])-1] != str(element[1]):
                            flag = False
                            break
                    if flag:
                        count += 1
                if ((count*100)/total) >= sup:
                    fqItemset.append(val)
                    countSet.append(count)
    return fqItemset, countSet

## Reading the data from the file
fileName = tkinter.filedialog.askopenfilename()
fileData = []
with open(fileName,'r') as f:
    for line in f:
        line = line.strip().split("\t")
        fileData.append(line)


row=len(fileData)
col=len(fileData[0])-1

#Creating numpy Matrix from the data set
data = np.zeros(shape=(row,col))
data = np.array(fileData)

labelItemSetUp = []
labelItemSetDown = []

# Label for features of the Genes
for i in range(0, data.shape[1]):
    labelItemSetUp.append('G' + str(i + 1) + '_Up')
    labelItemSetDown.append('G' + str(i + 1) + '_Down')

sup=[30,40,50,60,70]

for i in range (0, len(sup)):
    count =1;
    totalItem = 0
    print("Support is {}%".format(sup[i]))

    # calculating itemset for Length 1
    itemset, countSet = aprioriAlgo_fq1(data, sup[i], data.shape[0])
    print("Number of length-{} frequent Itemset: {}".format(count,len(itemset)))
    totalItem+=len(itemset)
    count += 1
    # calculating Itemset for length 2
    itemset2, countSet = aprioriAlgo_fq2(itemset, data, sup[i], data.shape[0])
    if len(itemset2) > 0:
        print("Number of length-{} frequent Itemset: {}".format(count, len(itemset2)))
        count+=1
        totalItem += len(itemset2)
    while(len(itemset2) > 1):
        itemset2, countSet = aprioriAlgo_fqAll(itemset2, data, sup[i], data.shape[0])
        if len(itemset2) > 0:
            print("Number of length-{} frequent Itemset: {}".format(count, len(itemset2)))
            count+=1
            totalItem += len(itemset2)
    print("Total Number of all ItemSet:{}".format(totalItem))
