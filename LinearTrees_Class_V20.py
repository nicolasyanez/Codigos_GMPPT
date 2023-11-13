# -*- coding: utf-8 -*-
"""
Created on Mon May 16 19:04:06 2022

@author: cristianguarnizo
"""
import numpy as np
from matplotlib import pyplot as plt  # now lets make some plots
plt.style.use("ggplot")


# spc_mat (0-2), MaxInd (3), Vmax (4), Pmax(5), localMaxV (6-8), localMaxP (9-11)

dataset = np.loadtxt("dataset.csv", delimiter=",")


clas = 1
index = [[0,3],[1,4],[2,5]]

plt.figure()
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


plt.plot(Xmax[:,0],Xmax[:,1], '+r')
plt.plot(X[:,0],X[:,1], 'ob',mfc = 'none')
X = np.vstack((X,Xmax))
y = np.concatenate((y, ymax))

## Decision trees based on linear classifiers

from sklearn.linear_model import RidgeClassifier
from lineartree import LinearTreeClassifier
from sklearn.datasets import make_classification
clf2 = LinearTreeClassifier(base_estimator=RidgeClassifier())
clf2.fit(X, y)

def plot_mesh(model, y, X, title):

    n_classes = 3
    plot_colors = "ryb"
    plot_step = 0.02

    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                          np.arange(y_min, y_max, plot_step))
    plt.tight_layout(h_pad=0.5, w_pad=0.5, pad=2.5)

    Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    cs = plt.contourf(xx, yy, Z, cmap=plt.cm.RdYlBu)

    plt.xlabel('features 0'); plt.ylabel('features 1')

    # Plot the training points
    for i, color in zip(range(n_classes), plot_colors):
        idx = np.where(y == i)
        plt.scatter(X[idx, 0], X[idx, 1], c=color,
                    cmap=plt.cm.RdYlBu, edgecolor='black', s=15)
    plt.title(title)
    
plt.figure(figsize=(8,6))
#plt.ylim([30,50])
plot_mesh(clf2, y, X, title="Class 20V")

plt.figure()
clf2.plot_model()

leaves = clf2.summary(only_leaves=True)

for m in leaves:
    if hasattr(leaves[m]['models'],'coef_'):
        print(f'Modelo {m}: ',leaves[m]['models'])
        print(f'Coeficientes: ',leaves[m]['models'].coef_)
        print('Intercepto: ',leaves[m]['models'].intercept_)
    else:
        print(f'Modelo {m}: ',leaves[m]['models'])

# #imprimir el arbol en codigo de C
# from micromlgen import port
# print(port(clf2))