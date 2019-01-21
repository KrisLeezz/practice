# -coding= utf-8 
import csv
import numpy
import math
from matplotlib import pyplot
from scipy.signal import argrelextrema

csv_file=open('G:/Git/Fourier/try/try/imfs.csv','r')
csv_read=csv.reader(csv_file)
imfs=[]
for row in csv_read:
	imfs.append(row)
imfs=numpy.array(imfs)
imfs=imfs.astype(float)
print imfs.shape
#significanceIMF
nDof=len(imfs[0])#number of data points=degree of freedom，数据的个数,列

#1.计算每个IMF的能量谱
columns=len(imfs)#imf的个数，行
logep=[0]*columns
for i in range (0,columns):
	logep[i]=0
	for j in range (0,nDof):
		logep[i]=logep[i]+imfs[i][j]*imfs[i][j]#第i个IMF中所有数的平方的和
	logep[i]=logep[i]/nDof
#2.将第一个IMF设置为白噪声，
sfactor =logep[0]
for i in range (0,columns):
	logep[i]=0.5636*logep[i]/sfactor
#3.计算每一个IMF的mean period,根据局部的极大值,所有的都是局部极大值
logep2=[0]*columns
for i in range(0,columns):
	x=imfs[i]
	max_position=argrelextrema(x, numpy.greater)#https://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array
	temp=numpy.size(max_position)
	logep2[i]=nDof/temp
print logep2
#4.乘以转换系数，1.4427=log(e)/log(2)
logep=numpy.array(logep)
logep.astype(float)
logep2=numpy.array(logep2)
logep2.astype(float)
for i in range (0,columns):
	logep[i]=1.4427*math.log(logep[i])
	logep2[i]=1.4427*math.log(logep2[i])

pyplot.figure()
pyplot.scatter(logep2,logep)
pyplot.show()
