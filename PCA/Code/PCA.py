import tkinter.filedialog
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD as SVD
from sklearn.manifold import TSNE as TSNE


########################################
#    PCA Algorithm
########################################

def algoPCA(data):
    # calculating the mean for each dimensions
    meanVect = np.mean(data, axis=0)
    dataCov = np.subtract(data, meanVect)
    covMat = np.cov(dataCov.T)
    eigenVal, eigenVect = np.linalg.eig(covMat)
    maxInd = 0
    secondMaxInd = 0;

    # Getting the two indexes with higher eigen values
    # and copying eigen Vectors across those indexes

    for i in range(1, len(eigenVal)):
        if (eigenVal[maxInd] > eigenVal[i]):
            if maxInd == secondMaxInd:
                secondMaxInd = i
            if eigenVal[i] > eigenVal[secondMaxInd]:
                secondMaxInd = i
        else:
            secondMaxInd = maxInd
            maxInd = i

    newEVect = np.zeros(shape=(eigenVect.shape[1], 2))
    newEVect[:, 0] = np.array(eigenVect[:, maxInd])
    newEVect[:, 1] = np.array(eigenVect[:, secondMaxInd])

    reduData = (np.dot(data, newEVect))
    return reduData


#############################################
# Utility function to plot the data
#############################################
def plotGraph(reduData, typeOfPlot, labels):
    setlabel = set(labels)
    f, ax = plt.subplots(figsize=(10, 10))
    for name in setlabel:
        x = reduData[labels[:] == name, 0]
        y = reduData[labels[:] == name, 1]
        ax.scatter(x, y, marker='o', label=name)

    plt.title(typeOfPlot + " plot for : " + str(fileName))
    plt.legend(ncol=1, fontsize=12)
    plt.show()


## Reading the data from the file
fileName = tkinter.filedialog.askopenfilename()
fileData = []
with open(fileName, 'r') as f:
    for line in f:
        line = line.strip().split("\t")
        fileData.append(line)

row = len(fileData)
col = len(fileData[0]) - 1

# Creating numpy Matrix from the data set
data = np.zeros(shape=(row, col))
data = np.array(fileData)
labels = (data[:, col])

data = np.delete(data, col, axis=1)
data = data.astype(float)

## Calling the PCA Algorithm
reducedData = algoPCA(data)

## Scatter Plot for the new Data
plotGraph(reducedData, "PCA", labels)

## Using package TruncatedSVD
svd = SVD(n_components=2)
svd_data = svd.fit_transform(data)
plotGraph(svd_data, "SVD", labels)

tsne = TSNE(n_components=2, init="pca", n_iter=500)
tsne_data = tsne.fit_transform(data)
plotGraph(tsne_data, "TNSE", labels)