import numpy as np
from LinReg import LinReg
import matplotlib.pyplot as plt
import os

true_mass = 71.1296-70.8639
conc = np.asarray([1,5,10,25,50, 100, 200, 1,1,1,5,5,5,10,10,10,25,25,25,50,50,50,100, 100, 200,200,200]) * true_mass * 1000/256.0
standard = np.asarray([0.004021097446, 0.01683960625, 0.03588608995, 0.09062817743, 0.27408, 0.601964, 1.10209, 0.00499760334, 0.005498753541, 0.004176975915, 0.01795402176, 0.01828315729, 0.0167947616, 0.03091238593, 0.03678892959, 0.03365368369, 0.08035973519, 0.08409251011, 0.08309703597, 0.175323413, 0.1868157196, 0.2012053372, 0.3985708183, 0.4225342248, 0.6136564896, 0.8452847497, 0.8014108039])
pinacolada = 0.00808965
peach = 0.578065
rainbow = 0.0874794
SB = 0.08952753
strawberry = 0.02141209

model = LinReg(standard, conc)
model.compute_model()

print "******************\nUnits = ug/mL\n*******************\n"

print "LOD/LOQ:"
model.calculate_LOD()
print "\n\n"

print 'Pina Colada: '
model.interpolate(pinacolada)

print 'Peaches and Cream: '
model.interpolate(peach)

print 'Rainbow: '
model.interpolate(rainbow)

print 'Strawberry Banana: '
model.interpolate(SB)

print 'Strawberry: '
model.interpolate(strawberry)

root_direc = "/Users/nicolasquach/Documents/stanford/senior_yr/spring/CHEM134/final_project/plots"
save_path = os.path.join(root_direc, 'regression.png')

model.plot_reverse(save_path)
model.compute_reverse_model()