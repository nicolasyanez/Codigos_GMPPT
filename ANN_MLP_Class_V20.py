import numpy as np
from matplotlib import pyplot as plt  # now lets make some plots
plt.style.use("ggplot")


# spc_mat (0-2), MaxInd (3), Vmax (4), Pmax(5), localMaxV (6-8), localMaxP (9-11)

dataset = np.loadtxt("dataset.csv", delimiter=",")


clas = 1
index = [[0,3],[1,4],[2,5]]


ind = dataset[:,3]==float(clas)
Xmax = dataset[ind,4:]
Xmax = Xmax[:,[0,1]]
ymax = dataset[ind,3]

Xlocal = dataset[~ind,6:]

ylocal = dataset[~ind,3]
Xlocal = Xlocal[:,index[clas]] 
ind0 = Xlocal[:,0] == 0.        #Indices donde el voltaje sea 0.0
X = np.delete(Xlocal,ind0,0)    #Se eliminan

y = np.delete(ylocal,ind0)


ind2 = y == 2.
y = np.delete(y,ind2)
X = np.delete(X,ind2,0) 


X = np.vstack((X,Xmax))
y = np.concatenate((y, ymax))

## Multi Layer Perceptron Clasiffier

from sklearn.neural_network import MLPClassifier
clf = MLPClassifier(solver='sgd', alpha=1e-5,hidden_layer_sizes=(3,3), random_state=1, activation='logistic', max_iter=2000)
clf = clf.fit(X, y)

print('Coeficiente:', clf.coefs_)
print('Bias:', clf.intercepts_)
print('N° de capas:', clf.n_layers_)
print('N° de entradas:', clf.n_features_in_)
print('N° de salidas:', clf.n_outputs_)