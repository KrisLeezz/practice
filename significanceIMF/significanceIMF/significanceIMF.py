# -coding= utf-8 
from __future__ import division
import csv
import numpy
import math
from matplotlib import pyplot
from scipy.signal import argrelextrema

#----------------------读取文件------------------------------------------
csv_file=open('G:/Git/Fourier/try/try/imfs.csv','r')
csv_read=csv.reader(csv_file)
imfs=[]
for row in csv_read:
	imfs.append(row)
imfs=numpy.array(imfs)
imfs=imfs.astype(float)
print imfs.shape
#----------------------完毕---------------------------------------------
#significanceIMF
nDof=len(imfs[0])#number of data points=degree of freedom，数据的个数,列
#1.计算每个IMF的能量谱
columns=len(imfs)#imf的个数，行
#------------------------计算能量谱En-------------------------------------
logep=[0]*columns#记录能量谱
for i in range (0,columns):#公式2.2
	logep[i]=0
	for j in range (0,nDof):
		logep[i]=logep[i]+imfs[i][j]*imfs[i][j]#第i个IMF中所有数的平方的和
	logep[i]=logep[i]/nDof
#-----------------------计算能量谱En完毕----------------------------------
#2.将第一个IMF设置为白噪声，进行标准化,参考K.Coughlin and K. K. Tung, 2004.
#高频信息往往包含的是噪声，the power in this mode can then be used to calibrate the moise distribution
sfactor =logep[0]
for i in range (0,columns):
	logep[i]=0.5636*logep[i]/sfactor#信噪比S/N，0.5636是什么玩意？？？？可能是一个平均值之类的
#3.----计算每一个IMF的mean period,根据相邻的极大值确定平均周期,所有的都是局部极大值-----
#  原代码IMF4之后采用了另外一种计算方式
logep2=[0]*columns
for i in range(0,columns):
	x=imfs[i]
	max_position=argrelextrema(x, numpy.greater)#https://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array
	temp=numpy.size(max_position)
	logep2[i]=nDof/temp
print logep2
#4.乘以转换系数，1.4427=log(e)/log(2)，自然对数？？？
logep=numpy.array(logep)
logep.astype(float)
logep2=numpy.array(logep2)
logep2.astype(float)
for i in range (0,columns):
	logep[i]=1.4427*math.log(logep[i])
	logep2[i]=1.4427*math.log(logep2[i])
#5.计算置信曲线
#计算PDF
def pdf_dist_value(yPos,yBar,nDof):
	#计算PDF
	ylen= len(yPos)
	#公式3.4
	eBar=numpy.exp(yBar)#Ebar
	evalue=numpy.exp(yPos)#E

	#计算每个y对应的pdf
	tmp3=[0]*ylen
	for iii in range (0,ylen):#公式3.4
		tmp1 = evalue[iii]/eBar-yPos[iii]#calculate---tmp1=E/(E-bar-y)
		tmp2 = -tmp1*nDof*eBar/2#calculate--------tmp2=(-1/2)*E-bar*Ndof*tmp1
		tmp3[iii] = 0.5*nDof*eBar*math.log(nDof) + tmp2
	#保证计算的收敛性，除以rscale，并没有改变PDF,PDF只是一个相对值，
	rscale=max(tmp3)
	tmp4=tmp3-rscale
	PDF=numpy.exp(tmp4)
	return PDF#每一段x中对应的5000个

def confidenceLine(percenta,Npt):
	#参数初始化
	percenta=percenta#置信度
	pdMax=round(math.log(nDof))+1 #x轴的最大值，IMF一定小于点的个数，ln对应lnX
	#x-ln(T)
	pdIntv=numpy.linspace(1,pdMax,100)#x轴分为100份
	#y-ln(E),满足y=-x，公式2.8
	yBar=-pdIntv
	#y的值服从卡方分布，对每个x来说求可能的y的，存在[yUpper,yLower]
	yUpper=[0]*100
	yLower=[0]*100
	for i in range(0,100):
		yUpper[i]=0#卡方从0开始
		yLower[i]=-3-pdIntv[i]*pdIntv[i]#chi-square decayed to zero,Wu use -3-xx as lower bound
	#计算每个x置信极限的位置
	sigline1=[0]*100#记录x-ln(T)
	sigline2=[0]*100#记录y-ln(E)
	for i in range (0,100):
		sigline1[i]=pdIntv[i]
		#将[yUpper,yLowe]分为5000段来研究白噪声的概率密度函数
		yPos=numpy.linspace(yUpper[i],yLower[i],5000)
		dyPos=yPos[1]-yPos[2]#记录区间间隔
		yPDF=pdf_dist_value(yPos,yBar[i],nDof)#调用pdf_dist_value计算pdf
		sum=0
		for jj in range (0,5000):
			sum=sum+yPDF[jj]#积累的PDF值
		#插值[yUpper,yLower],当给出一个精确的percenta，求对应的值
		jj1=0
		jj2=1
		psum1=0#置信上限
		psum2=yPDF[0]#这里为什么要选择第一个呢，是因为这是置信下限，挨个往下找
		pratio1=psum1/sum
		pratio2=psum2/sum
		#插值
		while pratio2<percenta:
			jj1=jj1+1
			jj2=jj2+1
			psum1=psum1+yPDF[jj1]
			psum2=psum2+yPDF[jj2]
			pratio1=psum1/sum
			pratio2=psum2/sum
			yref=yPos[jj1]
		sigline2[i]=yref + dyPos*(pratio2-percenta)/(pratio2-pratio1)#插值结果
		sigline2[i]= sigline2[i] + 0.066*pdIntv[i] + 0.12#原代码的意思‘confidence values assign to the line’
	#乘以转换系数1.4427
	sigline1=1.4427*numpy.asarray(sigline1)#ln(T)
	sigline2=1.4427*numpy.asarray(sigline2)#显著线的ln(E)
	return sigline2,sigline1

sigline95,pdIntv=confidenceLine(0.05,nDof)#95%
sigline90,pdIntv=confidenceLine(0.1,nDof)#90%
#画图
pyplot.figure()
p1=pyplot.scatter(logep2,logep,label='ln(E)')
p2=pyplot.plot(pdIntv,sigline95,'r.',label='sigline95')
p3=pyplot.plot(pdIntv,sigline90,'b-',label='sigline90')
pyplot.xlabel('ln(T)')
pyplot.ylabel('ln(E)')
pyplot.title('Significance IMF: ln(T)-ln(E)',fontsize='large')
pyplot.legend()
pyplot.show()
#为什么要1.4427？？？？