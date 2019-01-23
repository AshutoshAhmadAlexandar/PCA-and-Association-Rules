CSE-601 

Read Me - Association Analysis Part 2
=============================================
-Run the file apriori_part2.py using python:
	python apriori_part2.py
-You will be prompted to select the file:
	select filename
	e.g. choose file: associationruletestdata.txt
-The itemsets count will be displayed with support and confidence factors as follows:
	Support used-50 and Confidence used-70
	Number of length-1 frequent Itemset: 109
	Number of length-2 frequent Itemset: 63
	Number of length-3 frequent Itemset: 2
-A prompt to select the template number will appear as follows:
	Press 1 for Template 1:
	Press 2 for Template 2:
	Press 3 for Template 3:
	Press 4 to exit
-A prompt to select between (preloaded rules, customize rules) will appear as follows:
	Press 1 to use existing templates
	Press 2 to give new template
	Press 3 to go back
-For Existing templates preloaded rules results will be displayed.
-For Customized rules you will be prompted to enter the parameters based on the template
-Template 1 has 3 parameters:	
				First parameter: <RULE/BODY/HEAD>
				Second parameter: <ANY/NONE/1>
				Third parameter: <List of comma separated items>
-Template 2 has 2 parameters:	
				First parameter: <RULE/BODY/HEAD>
				Second parameter: <Number>
-Template 3 is a combination of Template 1 and Template 2 using AND or OR.
				First parameter: <Combination of template 1 and 2>
				e.g. 1or2, 1and2, 1and1, 2or2 etc.
-The rules which satisfy the entered template number and the parameters will be printed along with the count.