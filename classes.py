'''
Created on Jun 23, 2014

@author: giuliano
'''

class Step:
    def __init__(self, step=0, variable=None, value=None, energy = None):
        self.step = step
        self.variable = variable
        self.value = value
    def __str__(self):
        return " {0.step}        {0.variable}          {0.value}".format(self)
    
class Point:
    def __init__(self, step=None, coordinate=None, energy = None):
        self.step = step
        self.coordinate = coordinate
        self.energy = energy 
    def __str__(self):
        return " {0.step}                {0.coordinate}              {0.energy}".format(self)
    
class UB3LYP:
    def __init__(self, n = 0, E = None):
        self.E = E
        self.n = n
    def __str__(self):
            
            return " {0.n}            {0.E}".format(self)

class Coordinates:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    def __str__(self):
        return "{0.x}   {0.y}   {0.z}".format(self)
    
class Electron:
    def __init__(self, center_number, atomic_number, atomic_type, coordinates):
        self.center_number = center_number
        self.atomic_number = atomic_number
        self.atomic_type = atomic_type
        self.coordinates = coordinates
    def __str__(self):
        return "{0.center_number}          {0.atomic_number}       {0.atomic_type}   {0.coordinates}".format(self)
class Angle:
    def __init__(self, combinations='', angle=0):
        self.combinations = combinations
        self.angle = angle
    def __str__(self):
        return "{0.combinations}            {0.angle}".format(self)

class UBHFLYP:
    def __init__(self, n = 0, E = 0):
        self.E = E
        self.n = n
    def __str__(self):
        return " {0.n}            {0.E}".format(self)