#-*- coding: utf-8 -*-
import netCDF4
import datetime
from netCDF4 import Dataset
from netCDF4 import num2date
from matplotlib import pyplot
import numpy 
from mpl_toolkits.basemap import Basemap
import csv
ncdata=Dataset("mrro_Lmon_bcc-csm1-1_historical_r3i1p1_185001-201212.nc")
print ncdata
variable=ncdata.variables
runoff=ncdata.variables['mrro'][:]
time=ncdata.variables['time'][:]#366
lons=ncdata.variables['lon'][:]
lats=ncdata.variables['lat'][:]
runoff_unit=ncdata.variables['mrro'].units

dates=num2date(ncdata.variables['time'][:],ncdata.variables['time'].units)
print dates[-1]

lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(lat_0=lat_0, lon_0=lon_0)
lon, lat =numpy.meshgrid(lons, lats)#meshgrid()格网化
xi, yi = m(lon, lat)

for i in range(0,1955):
	cs = m.pcolor(xi, yi, numpy.squeeze(runoff[i]))#这里rhum_0表示的是最后一天6：00的rhum数据///pcolor()根据传入的data，绘制伪彩色图像//squeeze去除维度为1的维
	# Add Grid Lines
	# 绘制经纬线
	m.drawparallels(numpy.arange(-90., 91., 20.), labels=[1,0,0,0], fontsize=10)
	m.drawmeridians(numpy.arange(-180., 181., 40.), labels=[0,0,0,1], fontsize=10)

	# Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()

	# Add Colorbar
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	cbar.set_label(runoff_unit)
	# Add Title
	date=str(dates[i])
	pyplot.title('%s runoff'%date)
	#pyplot.savefig('%d_test.png'%i)
	pyplot.show()
	print i
	
print ('finish')
ncdata.close()
