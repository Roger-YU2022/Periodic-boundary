#!/usr/bin/env python
# _*_coding:UTF-8 _*_
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

executeOnCaeStartup()
from interaction import *


def Get_Cube_dimension(Model0914tension, Part1):  # 该函数返回rve的三维维尺寸范围
	node = mdb.models[Model0914tension].rootAssembly.instances[Part1].nodes
	Xmin = 0.0
	Xmax = 0.0
	Ymin = 0.0
	Ymax = 0.0
	Zmin = 0.0
	Zmax = 0.0
	for i in range(len(node)):
		x = node[i].coordinates[0]
		y = node[i].coordinates[1]
		z = node[i].coordinates[2]
		if (Xmin > x):
			Xmin = x
		elif (Xmax < x):
			Xmax = x
		if (Ymin > y):
			Ymin = y
		elif (Ymax < y):
			Ymax = y
		if (Zmin > z):
			Zmin = z
		elif (Zmax < z):
			Zmax = z
	return (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)


def Reference2D(Model0914tension):  # 创建参考点的函数
	from part import TWO_D_PLANAR, DEFORMABLE_BODY
	# Create reference parts and assemble
	NameRef1 = 'Ref-1'
	NameRef2 = 'Ref-2'
	p = mdb.models[Model0914tension].Part(dimensionality=TWO_D_PLANAR, name=NameRef1, type=
	DEFORMABLE_BODY)
	p.ReferencePoint(point=(5.0, -6.0, 0.0))
	p1 = mdb.models[Model0914tension].Part(dimensionality=TWO_D_PLANAR, name=NameRef2, type=
	DEFORMABLE_BODY)
	p1.ReferencePoint(point=(-5.0, -6.0, 0.0))
	mdb.models[Model0914tension].rootAssembly.Instance(dependent=ON, name=NameRef1,
	                                            part=mdb.models[Model0914tension].parts[NameRef1])
	mdb.models[Model0914tension].rootAssembly.Instance(dependent=ON, name=NameRef2,
	                                            part=mdb.models[Model0914tension].parts[NameRef2])
	# Create set of reference points
	mdb.models[Model0914tension].rootAssembly.Set(name=NameRef1, referencePoints=(
		mdb.models[Model0914tension].rootAssembly.instances[NameRef1].referencePoints[1],))
	mdb.models[Model0914tension].rootAssembly.Set(name=NameRef2, referencePoints=(
		mdb.models[Model0914tension].rootAssembly.instances[NameRef2].referencePoints[1],))


def fun_ne(Model0914tension, Part1, dimension_list):
	p = mdb.models[Model0914tension].rootAssembly
	# 找到边界点
	node = p.instances[Part1].nodes
	ne1 = []
	ne2 = []
	ne3 = []
	ne4 = []
	ne5 = []
	ne6 = []
	xmin = dimension_list[0]
	xmax = dimension_list[1]
	ymin = dimension_list[2]
	ymax = dimension_list[3]
	zmin = dimension_list[4]
	zmax = dimension_list[5]
	for i in range(len(node)):
		x = node[i].coordinates[0]
		y = node[i].coordinates[1]
		z = node[i].coordinates[2]
		if (abs(x - xmin) < tor):
			ne1.append(i)
		if (abs(x - xmax) < tor):
			ne2.append(i)
		if (abs(y - ymin) < tor):
			ne3.append(i)
		if (abs(y - ymax) < tor):
			ne4.append(i)
		if (abs(z - zmin) < tor):
			ne5.append(i)
		if (abs(z - zmax) < tor):
			ne6.append(i)
	print(len(ne1), len(ne2), len(ne3), len(ne4))
	return (ne1, ne2, ne3, ne4, ne5, ne6)

def periodic_2Dfun(Model0914tension, Part1, dimension_list, ne_all, num, ltf, rbn):
	p = mdb.models[Model0914tension]
	r = p.rootAssembly
	# 找到边界点
	xmin = dimension_list[0]
	xmax = dimension_list[1]
	ymin = dimension_list[2]
	ymax = dimension_list[3]
	zmin = dimension_list[4]
	zmax = dimension_list[5]
	if num == 1:
		lst1 = ne_all[0]
		lst2 = ne_all[1]
	elif num == 2:
		lst1 = ne_all[2]
		lst2 = ne_all[3]
	node = r.instances[Part1].nodes
	for n in range(min(len(lst1), len(lst2))):
		x0 = node[lst1[n]].coordinates[0]
		y0 = node[lst1[n]].coordinates[1]
		r.Set(nodes=node[lst1[n]:lst1[n] + 1], name='Set-' + str(ltf) + '-' + str(n + 1))
		min_distance = np.sqrt((xmin - xmax) ** 2 + (ymin - ymax) ** 2) * 2
		for j in range(len(lst2)):
			x1 = node[lst2[j]].coordinates[0]
			y1 = node[lst2[j]].coordinates[1]
			distance = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
			if (distance < min_distance):
				min_distance = distance
				index = j
		r.Set(nodes=node[lst2[index]:lst2[index] + 1], name='Set-' + str(rbn) + '-' + str(n + 1))
		if num == 2:
			if x0 == 0:  # 避免完全耦合，少约束一次，以施加边界条件
				continue
			else:
				p.Equation(name='Eq-' + str(ltf) + str(rbn) + '-X' + str(n + 1),
				           terms=((1, 'Set-' + str(ltf) + '-' + str(n + 1), 1),
				                  (-1, 'Set-' + str(rbn) + '-' + str(n + 1), 1)))
				p.Equation(name='Eq-' + str(ltf) + str(rbn) + '-Y' + str(n + 1),
				           terms=((1, 'Set-' + str(ltf) + '-' + str(n + 1), 2),
				                  (-1, 'Set-' + str(rbn) + '-' + str(n + 1), 2)))
		else:
			p.Equation(name='Eq-' + str(ltf) + str(rbn) + '-X' + str(n + 1),
			           terms=((1, 'Set-' + str(ltf) + '-' + str(n + 1), 1),
			                  (-1, 'Set-' + str(rbn) + '-' + str(n + 1), 1)))
			p.Equation(name='Eq-' + str(ltf) + str(rbn) + '-Y' + str(n + 1),
			           terms=((1, 'Set-' + str(ltf) + '-' + str(n + 1), 2),
			                  (-1, 'Set-' + str(rbn) + '-' + str(n + 1), 2)))


# --------------------Main protram--------------
# fields = (
# 	('x min:', '0.0'), ('x max:', '1.0'), ('y min:', '0.0'), ('y max:', '1.0'))
# xmin, xmax, ymin, ymax = getInputs(fields=fields, label="Input the parameters")
# xmin = float(xmin)
# xmax = float(xmax)
# ymin = float(ymin)
# ymax = float(ymax)

tor = 0.00001

Model0914tension, Part1 = getInputs(
	fields=(('Model Name:', 'Model-1'), ('Instance Name:', 'Part-1-mesh-1-1')),
	label='Enter model Name and instance Name',
	dialogTitle='Enter the names.')
Dimension = Get_Cube_dimension(Model0914tension, Part1)
print('Xmin,Xmax,Ymin,Ymax,Zmin,Zmax=', Dimension)
# 创建并建立参考点的set
Reference2D(Model0914tension)
ne_all = fun_ne(Model0914tension, Part1, Dimension)
# 定义Left/Right 边界节点的PBC
periodic_2Dfun(Model0914tension, Part1, Dimension, ne_all, 1, 'L', 'R')
# 定义Top/Bottom 边界节点的PBC
periodic_2Dfun(Model0914tension, Part1, Dimension, ne_all, 2, 'T', 'B')
