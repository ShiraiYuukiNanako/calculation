import formula as f
import math

l1 = [3.81,0.1,3.98,0.95]
l2 = [1.7*math.pow(10,-6),0.546,3.89*math.pow(10,-5),8.72*math.pow(10,-3)]
C = 1*math.pow(10,-3)

print(f.free_metal(C,l1,l2))