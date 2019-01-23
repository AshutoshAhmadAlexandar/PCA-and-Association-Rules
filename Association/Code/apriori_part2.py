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
'''
#############################################################################
#    Association Rule
#############################################################################
'''

'''
    Function to calculate the Cobinations for the
    RULE Generation using itemset
'''
def combinations(target, data):
    for i in range(len(data)):
        new_target = copy.copy(target)
        new_data = copy.copy(data)
        new_target.append(data[i])
        new_data = data[i + 1:]
        newTarget.append(new_target)
        combinations(new_target,new_data)

'''
    Function to Calculate the Confidence 
    of a Rule and send 'True' or 'False'
    if the itemset is frequent or not frequent 
'''
def confidenceCalc(body, item, confFactor):
    body_string = ','.join(body)
    item_string = ','.join(item)
    if float((int(res_dict[item_string])/int(res_dict[body_string]))*100) >= confFactor:
        return True
    else:
        return False


count = 1
sup = 50
conf = 70

print("Support used-{} and Confidence used-{}".format(sup,conf))
feqItemListSupported=[]
countOfItemSupported=[]
itemset, countSet = aprioriAlgo_fq1(data, sup, data.shape[0])
print("Number of length-{} frequent Itemset: {}".format(count,len(itemset)))
count += 1
feqItemListSupported.append(itemset)
countOfItemSupported.append(countSet)
# calculating Itemset for length 2
itemset2, countSet = aprioriAlgo_fq2(itemset, data, sup, data.shape[0])
if len(itemset2) > 0:
    print("Number of length-{} frequent Itemset: {}".format(count, len(itemset2)))
    count+=1
    feqItemListSupported.append(itemset2)
    countOfItemSupported.append(countSet)
while(len(itemset2) > 1):
    itemset2, countSet = aprioriAlgo_fqAll(itemset2, data, sup, data.shape[0])
    if len(itemset2) > 0:
        print("Number of length-{} frequent Itemset: {}".format(count, len(itemset2)))
        count+=1
        feqItemListSupported.append(itemset2)
        countOfItemSupported.append(countSet)

#print(feqItemListSupported)
res_dict = {}
for i in range (len(feqItemListSupported)):
    for j in range(len(feqItemListSupported[i])):
        s = ','.join(feqItemListSupported[i][j])
        res_dict.setdefault(s, 0)
        res_dict[s] = countOfItemSupported[i][j]
#print (res_dict)

#print(len(feqItemListSupported))
combination=[]
headList=[]
bodyList=[]
for i in range(1,len(feqItemListSupported)):
    for j in range(0, len(feqItemListSupported[i])):
        newTarget =[]
        combinations(combination, feqItemListSupported[i][j])
        head=[]
        body=[]
        for ele in newTarget:
            val = [item for item in feqItemListSupported[i][j] if item not in ele]
            if val:
                ret = confidenceCalc(val, feqItemListSupported[i][j], conf)
                if ret:
                    body.append(ele)
                    head.append(val)
        headList.append(head)
        bodyList.append(body)


'''
    Template 1 function 
'''
def template_1(type1, type2, list):
    countConf = 0
    rule = []
    if (type1.upper() == 'RULE'):
        for i in range (0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                count=0
                for p in range (0, len(list)):
                    if list[p] in bodyList[i][j] or list[p] in headList[i][j]:
                        count+=1
                if (str(type2.upper()) == '1'):
                    if count == 1:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j],headList[i][j]))
                elif (str(type2.upper()) == 'NONE'):
                    if count == 0:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
                elif str(type2.upper()) == 'ANY':
                    if count >= 1:
                        countConf += 1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    elif (type1.upper() == 'BODY'):
        for i in range (0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                count=0
                for p in range (0, len(list)):
                    if list[p] in bodyList[i][j]:
                        count+=1
                if (str(type2.upper()) == '1'):
                    if count == 1:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j],headList[i][j]))
                elif (str(type2.upper()) == 'NONE'):
                    if count == 0:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
                elif str(type2.upper()) == 'ANY':
                    if count >= 1:
                        countConf += 1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    elif (type1.upper() == 'HEAD'):
        for i in range (0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                count=0
                for p in range (0, len(list)):
                    if list[p] in headList[i][j]:
                        count+=1
                if (str(type2.upper()) == '1'):
                    if count == 1:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j],headList[i][j]))
                elif (str(type2.upper()) == 'NONE'):
                    if count == 0:
                        countConf+=1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
                elif str(type2.upper()) == 'ANY':
                    if count >= 1:
                        countConf += 1
                        rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    return rule, countConf
'''
    Template 2 Function
'''
def template_2(type1, type2):
    countConf = 0
    rule = []
    if (type1.upper() == 'RULE'):
        for i in range(0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                if (len(bodyList[i][j]) + len(headList[i][j])) > int(type2)-1:
                    countConf += 1
                    rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    elif (type1.upper() == 'BODY'):
        for i in range(0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                if (len(bodyList[i][j])) > int(type2)-1:
                    countConf += 1
                    rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    elif (type1.upper() == 'HEAD'):
        for i in range(0, len(bodyList)):
            for j in range(0, len(bodyList[i])):
                if (len(headList[i][j])) > int(type2)-1:
                    countConf += 1
                    rule.append("{} => {}".format(bodyList[i][j], headList[i][j]))
    return rule, countConf




'''
####################################################################
#                   Main loop for the User Interface   
####################################################################
'''

run_program = True

while run_program:
    print('Press 1 for Template 1:')
    print('Press 2 for Template 2:')
    print('Press 3 for Template 3:')
    print('Press 4 to exit')
    user_input = input()
    if int(user_input) == 1:
        inner_program = True
        print('Entering template 1:')
        while inner_program:
            print('Press 1 to use existing templates')
            print("Press 2 to give new template")
            print("Press 3 to go back")
            user_input2 = input()
            if int(user_input2) == 1:
                rule, confCount = template_1("RULE", "ANY", ["G59_Up"])
                print("RULE Has ANY of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("RULE", "NONE", ["G59_Up"])
                print("RULE Has NONE of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("RULE", "1", ["G59_Up", "G10_Down"])
                print("RULE Has 1 of [G59_Up, G10_Down]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("HEAD", "ANY", ["G59_Up"])
                print("HEAD Has ANY of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("HEAD", "NONE", ["G59_Up"])
                print("HEAD Has NONE of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("HEAD", "1", ["G59_Up", "G10_Down"])
                print("HEAD Has 1 of [G59_Up, G10_Down]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("BODY", "ANY", ["G59_Up"])
                print("BODY Has ANY of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("BODY", "NONE", ["G59_Up"])
                print("BODY Has NONE of [G59_Up]")
                print(rule)
                print("The Count of the rules are ", confCount)

                rule, confCount = template_1("BODY", "1", ["G59_Up", "G10_Down"])
                print("BODY Has 1 of [G59_Up, G10_Down]")
                print(rule)
                print("The Count of the rules are ", confCount)

            elif int(user_input2) == 2:
                inner_check1 = True
                while inner_check1:
                    first_part = input('Enter [RULE/HEAD/BODY]')
                    if str(first_part).upper() not in ['RULE', 'BODY', 'HEAD']:
                        inner_check = False
                    else:
                        inner_check2 = True
                        while inner_check2:
                            second_part = input('Enter [ANY/NONE/1]')
                            if str(second_part).upper() not in ['ANY','NONE','1']:
                                inner_check2 = False
                            else:
                                third_part = input('Enter elements separate by comma(do not put space):')
                                third_part = list(third_part.split(','))
                                res = []
                                for i in third_part:
                                    a, b = i.split('_')
                                    c = a.title() + '_' + b.title()
                                    res.append(c)
                                print('New Template is: {} - {} - {}'.format(first_part,second_part,third_part))
                                inner_check1 = False
                                inner_check2 = False
                                rule, confCount = template_1(str(first_part), str(second_part), res)
                                print(rule)
                                print("The Count of the rules are ", confCount)
            elif int(user_input2) == 3:
                inner_program = False
    elif int(user_input) == 2:
        inner_program = True
        print('Entering template 2:')
        while inner_program:
            print('Press 1 to use existing templates')
            print("Press 2 to give new template")
            print("Press 3 to go back")
            user_input2 = input()
            if int(user_input2) == 1:
                rule, confCount = template_2("RULE", "3")
                print("RULE, 3")
                print(rule)
                print("The Count of the rules are ", confCount)
                rule, confCount = template_2("HEAD", "2")
                print("HEAD, 2")
                print(rule)
                print("The Count of the rules are ", confCount)
                rule, confCount = template_2("BODY", "1")
                print("BODY, 1")
                print(rule)
                print("The Count of the rules are ", confCount)
            elif int(user_input2) == 2:
                inner_check1 = True
                while inner_check1:
                    first_part = input('Enter [RULE/HEAD/BODY]')
                    if str(first_part).upper() not in ['RULE', 'BODY', 'HEAD']:
                        inner_check = False
                    else:
                        inner_check2 = True
                        while inner_check2:
                            second_part = input('Enter [Size]')
                            if not second_part.isdigit():
                                inner_check2 = False
                            else:
                                print('New Template is: {} - {}'.format(first_part,second_part))
                                inner_check1 = False
                                inner_check2 = False
                                rule,confCount = template_2(first_part,second_part)
                                print(rule)
                                print("The Count of the rules are ", confCount)
            elif int(user_input2) == 3:
                inner_program = False
    elif int(user_input) == 3:
        inner_program = True
        print('Entering template 3:')
        while inner_program:
            print('Press 1 to use existing templates')
            print("Press 2 to give new template")
            print("Press 3 to go back")
            user_input2 = input()
            if int(user_input2) == 1:
                rule1, confCount1 = template_1("HEAD", "ANY", ["G10_Down"])
                rule2, confCount2 = template_1("BODY", "1", ["G59_Up"])
                rule3 = set(rule1).union(set(rule2))
                print("1or1 HEAD has ANY of [G10_Down] or BODY has 1 of [G59_Up]")
                print(rule3)
                print("The Count of the rules are ", len(rule3))

                rule1, confCount1 = template_1("HEAD", "ANY", ["G10_Down"])
                rule2, confCount2 = template_1("BODY", "1", ["G59_Up"])
                rule3 = set(rule1).intersection(set(rule2))
                print("1and1 HEAD has ANY of [G10_Down] and BODY has 1 of [G59_Up]")
                print(rule3)
                print("The Count of the rules are ", len(rule3))

                rule1, confCount1 = template_1("HEAD", "ANY", ["G10_Down"])
                rule2, confCount2 = template_2("BODY", "2")
                rule3 = set(rule1).union(set(rule2))
                print("1or2 HEAD has ANY of [G10_Down] or BODY 2")
                print(rule3)
                print("The Count of the rules are ", len(rule3))

                rule1, confCount1 = template_1("HEAD", "ANY", ["G10_Down"])
                rule2, confCount2 = template_2("BODY", "2")
                rule3 = set(rule1).intersection(set(rule2))
                print("1and2 HEAD has ANY of [G10_Down] and BODY 2")
                print(rule3)
                print("The Count of the rules are ", len(rule3))

                rule1, confCount1 = template_2("HEAD", "1")
                rule2, confCount2 = template_2("BODY", "2")
                rule3 = set(rule1).union(set(rule2))
                print("2or2 HEAD 1 or BODY 2")
                print(rule3)
                print("The Count of the rules are ", len(rule3))

                rule1, confCount1 = template_2("HEAD", "1")
                rule2, confCount2 = template_2("BODY", "2")
                rule3 = set(rule1).intersection(set(rule2))
                print("2and2 HEAD 1 and BODY 2")
                print(rule3)
                print("The Count of the rules are ", len(rule3))
            elif int(user_input2) == 2:
                inner_check1 = True
                while inner_check1:
                    first_input = input('Enter combination: [1or1, 1and1, 1or2, 1and2, 2or2, 2and2]')
                    if str(first_input) not in ['1or1', '1and1', '1or2', '1and2', '2or2', '2and2']:
                        inner_check = False
                    else:
                        inner_check2 = True
                        a, b, c = first_input[0], first_input[-1], first_input[1:len(first_input)-1]
                        while inner_check2:
                            if a == '1':
                                first_part = input('Enter [RULE/HEAD/BODY]-[ANY/NONE/1]-elem2,elem2,etc')
                                first_list = first_part.split('-')
                                t = list(first_list[2].split(','))
                                rule1, confCount1 = template_1(first_list[0], first_list[1],t)
                                if int(b) == 1:
                                    second_part = input('Enter [RULE/HEAD/BODY]-[ANY/NONE/1]-elem2,elem2,etc')
                                    second_list = second_part.split('-')
                                    t1 = list(second_list[2].split(','))
                                    rule2, confCount2 = template_1(second_list[0], second_list[1], t1)
                                    print('Template is: {} - {} - {} - {} - {} - {} - {}'.format(first_input,
                                                                                                 first_list[0],
                                                                                                 first_list[1],
                                                                                                 first_list[2],
                                                                                                 second_list[0],
                                                                                                 second_list[1],
                                                                                                 second_list[2]))
                                    if c == 'or':
                                        rule3 = set(rule1).union(set(rule2))
                                    elif c == 'and':
                                        rule3 = set(rule1).intersection(set(rule2))
                                    print(rule3)
                                    print("The Count of the rules are ", len(rule3))

                                    inner_check2 = False
                                    inner_check1 = False
                                elif b == '2':
                                    second_part = input('Enter [RULE/HEAD/BODY]-Size')
                                    second_list = second_part.split('-')
                                    rule2, confCount2 = template_2(second_list[0], second_list[-1])
                                    print('Template is: {} - {} - {} - {} - {} - {}'.format(first_input, first_list[0],
                                                                                            first_list[1],
                                                                                            first_list[2],
                                                                                            second_list[0],
                                                                                            second_list[-1]))
                                    if c == 'or':
                                        rule3 = set(rule1).union(set(rule2))
                                    elif c == 'and':
                                        rule3 = set(rule1).intersection(set(rule2))
                                    print(rule3)
                                    print("The Count of the rules are ", len(rule3))
                                    inner_check2 = False
                                    inner_check1 = False

                            elif a == '2':
                                first_part = input('Enter [RULE/HEAD/BODY]-Size')
                                first_list = first_part.split('-')
                                rule1, confCount1 = template_2(first_list[0], first_list[-1])
                                if b == '1':
                                    second_part = input('Enter [RULE/HEAD/BODY]-[ANY/NONE/1]-elem2,elem2,etc')
                                    second_list = second_part.split('-')
                                    t = list(second_list[2].split(','))
                                    rule2, confCount2 = template_1(second_part[0], second_part[1], t)
                                    print('Template is: {} - {} - {} - {} - {} - {}'.format(first_input, first_list[0],
                                                                                            first_list[-1],
                                                                                            second_list[0],
                                                                                            second_list[1],
                                                                                            second_list[2]))
                                    if c == 'or':
                                        rule3 = set(rule1).union(set(rule2))
                                    elif c == 'and':
                                        rule3 = set(rule1).intersection(set(rule2))
                                    print(rule3)
                                    print("The Count of the rules are ", len(rule3))

                                    inner_check2 = False
                                    inner_check1 = False
                                elif b == '2':
                                    second_part = input('Enter [RULE/HEAD/BODY]-Size')
                                    second_list = second_part.split('-')
                                    rule2, confCount2 = template_2(second_list[0], second_list[-1])

                                    print('Template is: {} - {} - {} - {} - {}'.format(first_input, first_list[0],
                                                                                            first_list[-1],
                                                                                            second_list[0],
                                                                                            second_list[-1]))
                                    if c == 'or':
                                        rule3 = set(rule1).union(set(rule2))
                                    elif c == 'and':
                                        rule3 = set(rule1).intersection(set(rule2))
                                    print(rule3)
                                    print("The Count of the rules are ", len(rule3))

                                    inner_check2 = False
                                    inner_check1 = False
            elif int(user_input2) == 3:
                inner_program = False
    elif int(user_input) == 4:
        run_program = False