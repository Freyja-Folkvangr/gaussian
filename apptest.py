# code's author: Giuliano Tognarelli Buono-core
import logging
import sys
from tkinter import *
import time


global verbose
verbose = True

global file
file = "Select a file to begin"


old_stdout = sys.stdout

class StdoutRedirector(object):

    def __init__(self, text_area):
        self.text_area = text_area

    def write(self, str):
        self.text_area.insert(END, str)
        self.text_area.see(END)

def getline(file, line_num):
    if line_num < 1: return ''
    current_line = 0
    for line in open(file):
        current_line += 1
        if current_line == line_num: return line
    return ''

def checkfile(afile):
    import os.path
    if os.path.exists(afile):
        x, *y = afile.split('.')
        if y[0] != "log":
                load_file()
        global file
        file = afile
        return True
    else:
        load_file()

def load_file():
        global file
        from tkinter.filedialog import askopenfilename
        from tkinter.messagebox import showerror
        try:
                fname = askopenfilename(filetypes=(("Gaussian log/output files", "*.log"),
                                           ("All files", "*.*")))
        except (RuntimeError, IOError) as inst:
                print ("There was an error while oppening a {1}".format(file))
                print(type(inst))
                print(inst)
        if fname:
            import os.path
            if os.path.exists(fname):
                        file = fname
                        global textBox1_v
                        textBox1_v.set(file)
            else:
                return False

def op():
        class Coordinates:
            def __init__(self, x=0, y=0, z=0):
                self.x = x
                self.y = y
                self.z = z
                    
            def __str__(self):
                return "<{0.x}, {0.y}, {0.z}>".format(self)

        class Electron:
                def __init__(self, center_number, atomic_number, atomic_type, coordinates):
                        self.center_number = center_number
                        self.atomic_number = atomic_number
                        self.atomic_type = atomic_type
                        self.coordinates = coordinates
                def __str__(self):
                        return "{0.center_number}          {0.atomic_number}       {0.atomic_type}   {0.coordinates}".format(self)
        def report_Standard_orientation(Standard_orientation, Matrix_number):
            if verbose == True: print("================STANDARD ORIENTATION {}================".format(Matrix_number))
            else: print("================STANDARD ORIENTATION================")
            print("NUMBER    Atom    Type     Coords")
            for item in Standard_orientation:
                print(item)
        def log_Standard_orientation(Standard_orientation, Matrix_number):
            if verbose == True:
                    log.write("================STANDARD ORIENTATION {}================\n".format(Matrix_number))
                    log.write("NUMBER    Atom    Type     Coords\n")
                    for item in Standard_orientation:
                            log.write("{}\n".format(item))
        def report_orbitals(aocc, bocc, avirt, bvirt):
                try:
                    global n
                    n = n + 1
                    print("=======================RESULTS======================")
                    if verbose == True: print("Itineration {}:".format(n))
                    if verbose == True:
                        log.write("\nAO:\n")
                        log.write(str(aocc))
                        log.write("\n\nAV:\n")
                        log.write(str(avirt))
                        log.write("\n\nBO:\n")
                        log.write(str(bocc))
                        log.write("\n\nBV:\n")
                        log.write(str(bvirt))
                        log.write("\n=============================EOI=============================\n")
                    if bocc != []: print ("Alpha results")
                    print("HOMO: {}".format(aocc[len(aocc) - 1]))
                    print("LUMO: {}".format(avirt[0]))
                    print("μ= {0}   Ηη= {1}".format(((avirt[0] + aocc[len(aocc) - 1]) * 0.5), ((avirt[0] - aocc[len(aocc) - 1]) * 0.5)))
                    print()
                    if bocc != []:
                        print ("Beta results")
                        print("HOMO: {}".format(bocc[len(bocc) - 1]))
                        print("LUMO: {}".format(bvirt[0]))
                        print("μ= {0}   Ηη= {1}".format(((bvirt[0] + bocc[len(bocc) - 1]) * 0.5), ((bvirt[0] - bocc[len(bocc) - 1]) * 0.5)))
                except (RuntimeError, TypeError, ValueError) as inst:
                        print ("There was an error while reporting the results")
                        print (type(inst))
                        print (inst)
        try:
            global n
            n = 0
            Line_number = 0
            global file
            global verbose
            if verbose == True:
                global log
                log = open("gauss09sAWK.log", "w")
                from datetime import datetime
                log.write("=============================Gauss09 sAWK=============================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
                log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
                log.close()
                log = open("gauss09sAWK.log", "a")
                print("NOTE: Logs are turned on")
                print("NOTE 2: Saving Gauss09 sAWK logs in 'gauss09sAWK.log'")    
            checkfile(file)
            with open(file, "r") as f:
                    Energy = (None, None, 0) #[E, type, found multiple HF values?] Types: 1=Hartree-Fock
                    aocc = []
                    bocc = []
                    avirt = []
                    bvirt = []
                    Standard_orientation = []
                    Matrix_number = 0
                    for line in f:
                        Line_number += 1
                        if "Alpha  occ. eigenvalues" in line:
                            x, *y = line.split(' --  ')
                            args = y[0].split('  ')
                            if verbose == True: log.write("AO arguments: {}\n".format(args))  # AO means alpha occ eigenvalues
                            for item in args:
                                    if item != '':
                                            try:
                                                    aocc.append(float(item))
                                            except (ValueError) as err:
                                                    if verbose == True:
                                                            log.write("\nWarning: There was an expected format-related-problem with data treatment\n")
                                                    a, *b = item.split(' ')
                                                    aocc.append(float(a))
                                                    aocc.append(float(b[0]))
                                            except:
                                                   print ("Unknown Error!")
                                                   print (type(err))
                                                   print (err)
                        elif "Alpha virt. eigenvalues" in line:
                            x, *y = line.split(' --  ')
                            args = y[0].split('  ')
                            if verbose == True: log.write("AV arguments: {}\n".format(args))  # AV means alpha virt eigenvalues
                            for item in args:
                                    if item != '':
                                            try:
                                                    avirt.append(float(item))
                                            except (ValueError) as err:
                                                    if verbose == True:
                                                            log.write("\nWarning: There was an expected format-related-problem with data treatment\n")
                                                    a, *b = item.split(' ')
                                                    aocc.append(float(a))
                                                    aocc.append(float(b[0]))
                                            except:
                                                   print ("Unknown Error!")
                                                   print (type(err))
                                                   print (err)
                        elif "Beta  occ. eigenvalues" in line:
                            x, *y = line.split(' --  ')
                            args = y[0].split('  ')
                            if verbose == True: log.write("BO arguments: {}\n".format(args))  # BO means beta occ eigenvalues
                            for item in args:
                                if item != '':
                                        try:
                                                bocc.append(float(item))
                                        except(ValueError) as err:
                                                if verbose == True:
                                                        log.write("\nWarning: there was an expected error format-related-problem with data treatment\n")
                                                a, *b = item.split(' ')
                                                aocc.append(float(a))
                                                aocc.append(float(b[0]))
                                        except:
                                                print ("Unknown error!")
                                                print (type(err))
                                                print (err)

                        elif "Beta virt. eigenvalues" in line:
                            x, *y = line.split(' --  ')
                            args = y[0].split('  ')
                            if verbose == True: log.write("BV arguments: {}\n".format(args))  # BV means beta virt eigenvalues
                            for item in args:
                                    if item != '':
                                            try:
                                                    bvirt.append(float(item))
                                            except (ValueError) as err:
                                                    if verbose == True:
                                                            log.write("\nWarning: There was an expected format-related-problem with data treatment\n")
                                                            log.write("{1}\n{2}\nDone some automatic fixes..\n".format(type(err), err))
                                                    a, *b = item.split(' ')
                                                    aocc.append(float(a))
                                                    aocc.append(float(b[0]))
                                            except:
                                                   print ("Unknown Error!")
                                                   print (type(err))
                                                   print (err)
                                                   
                        elif "Standard orientation" in line:
                            import linecache
                            j = 5
                            Matrix_number += 1
                            while "-----" not in linecache.getline(file, Line_number + j):
                                q, *r = linecache.getline(file, Line_number + j).split(' ')
                                q = []
                                for item in r:
                                        if item != '':
                                                q.append(item)
                                coordinates = Coordinates(float(q[3]), float(q[4]), float(q[5]))
                                electron = Electron(int(q[0]), int(q[1]), int(q[2]), coordinates)
                                Standard_orientation.append(electron)
                                j = j + 1
                            
                        elif "HF=" in line:
                                x, *y = line.split('HF=')
                                x, *y = y[0].split('\\')
                                try:
                                    tmp=float(x)
                                except (ValueError):
                                    p, *q = x.split(",")
                                    x = p
                                try:
                                    if Energy[0] != None and Energy[1] != None and Energy[2] == 0:
                                        Energy = (float(x), 1, 1)
                                        if verbose == True:
                                            log.write("HF: {}\n".format(Energy[0]))
                                    else:
                                        Energy = (float(x), 1, 0)
                                        if verbose == True:
                                            log.write("HF: {}\n".format(Energy[0]))
                                except (ValueError) as err:
                                    print("{} in HF (E)".format(type(err)))
                                    if verbose == True:
                                        log.write("{1} in HF (E)\nError: {2}\n Tuple: {3}\n".format(type(err), err, Energy))
                                except (TypeError) as err:
                                    print("{} in HF (E)".format(type(err)))
                                    if verbose == True:
                                        log.write("{1} in HF (E)\nError: {2}\n Tuple: {3}\n".format(type(err), err, Energy))
                                
                        elif "GradGradGradGradGradGradGrad" in line or "Initial guess <" in line:
                                if verbose == True and (aocc != [] and avirt != []):
                                        report_orbitals(aocc, bocc, avirt, bvirt)
                                if Standard_orientation != [] and verbose == True:
                                        report_Standard_orientation(Standard_orientation, Matrix_number)
                                        log_Standard_orientation(Standard_orientation, Matrix_number)
                                Standard_orientation = []
                                aocc = []
                                bocc = []
                                avirt = []
                                bvirt = []
                                                
                        elif "Normal termination of Gaussian" in line:
                                if aocc != [] and avirt != []:
                                        report_orbitals(aocc, bocc, avirt, bvirt)
                                if Standard_orientation != []: report_Standard_orientation(Standard_orientation, Matrix_number)
                                if verbose == True: log_Standard_orientation(Standard_orientation, Matrix_number)
                                Standard_orientation = []
                                aocc = []
                                bocc = []
                                avirt = []
                                bvirt = []
                                
                                
                        else: pass
                    log.close()
                    if Energy[1] == 1: print("Hartree-Fock= {}".format(Energy[0]))
                    if Energy[2] == 1: print("-Many HF found, see details in log file")
                    print("Finished")
        except(RuntimeError, TypeError) as inst:
                print ("There was an error while reading the optimization values")
                print (type(inst))
                print (inst)

def irc():
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
        
        class Coordinates:
            def __init__(self, x=0, y=0, z=0):
                self.x = x
                self.y = y
                self.z = z
                
            def __str__(self):
                return "<{0.x}, {0.y}, {0.z}>".format(self)

        class Electron:
                def __init__(self, center_number, atomic_number, atomic_type, coordinates):
                        self.center_number = center_number
                        self.atomic_number = atomic_number
                        self.atomic_type = atomic_type
                        self.coordinates = coordinates
                def __str__(self):
                        return "{0.center_number}          {0.atomic_number}       {0.atomic_type}      {0.coordinates}".format(self)
        def report_zMatrix(zMatrix, Matrix_number):
                if verbose == True: print("======================Z-MATRIX {}======================".format(Matrix_number))
                else: print("======================Z-MATRIX======================")
                print("NUMBER    Atom    Type      Coords")
                for item in zMatrix:
                        print(item)
        def log_zMatrix(zMatrix, Matrix_number):
            if verbose == True:
                    log.write("======================Z-MATRIX {}======================\n".format(Matrix_number))
                    log.write("NUMBER    Atom    Type      Coords\n")
                    for item in zMatrix:
                            log.write("{}\n".format(item))
            else: pass
        def report_angles(internal_angles, Matrix_number):
            if verbose == True: print("==================INTERNAL-ANGLES {}===================".format(Matrix_number))
            else: print("==================INTERNAL-ANGLES===================")
            print("ANGLE                 VALUE")
            for item in internal_angles:
                print(item)
        def log_angles(internal_angles, Matrix_number):
            log.write("==================INTERNAL-ANGLES {}===================".format(Matrix_number))
            log.write("ANGLE              VALUE\n")
            for item in internal_angles:
                log.write("{}\n".format(item))
        try:
                global n
                n = 0
                Matrix_number = 0
                global file
                global verbose
                if verbose == True:
                        global log
                        log = open("gauss09sAWK.log", "w")
                        from datetime import datetime
                        log.write("======================Gaussian09 simple AWK======================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
                        log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
                        log.close()
                        log = open("gauss09sAWK.log", "a")
                        print("NOTE: Logs are turned on")
                        print("NOTE 2: Saving Gauss09 sAWK logs in 'gauss09sAWK.log'")    
                checkfile(file)
                with open(file, "r") as f:
                    internal_angles = []
                    zMatrix = []
                    ubhflyp=[int(0)]
                    Energy = (None, None, 0) #[E, type, found multiple HF values?] Types: 1=Hartree-Fock
                    for line in f:
                            n = n + 1
                            if "Z-Matrix orientation" in line:
                                    Matrix_number = Matrix_number + 1
                                    import linecache
                                    j = 6
                                    while "-----" not in linecache.getline(file, n + j):
                                            q, *r = linecache.getline(file, n + j).split(' ')
                                            q = []
                                            for item in r:
                                                    if item != '':
                                                            q.append(item)
                                            coordinates = Coordinates(float(q[3]), float(q[4]), float(q[5]))
                                            electron = Electron(int(q[0]), int(q[1]), int(q[2]), coordinates)
                                            zMatrix.append(electron)
                                            j = j + 1
                            elif "Interatomic angles" in line:
                                    import linecache
                                    j = n + 1
                                    while "     " in linecache.getline(file, j):
                                            x, *y = linecache.getline(file, j).split('  ')
                                            for item in y:
                                                if item != '':
                                                    w, *z = item.split('=')
                                                    if ' ' in w:
                                                        a, *b = w.split(' ')
                                                        w=b[0]
                                                    if (z != [] and z != '' and z != '\n') and ' ' in z[0]:
                                                        a, *b = z[0].split(' ')
                                                        if a != '' and a != '\n': z[0] = a
                                                        else: z[0] = b[0]
                                                    if z == []: pass
                                                    else: internal_angles.append(Angle(w, float(z[0])))
                                            j += 1
                                            
                            elif "SCF Done:" in line and "E(UB+HF-LYP)" in line:
                                x, *y = line.split("=")
                                #print(y)
                                x, *y = y[0].split(" ")
                                for item in y:
                                    if item == "": y.remove(item)
                                e=UBHFLYP(ubhflyp[0] +1, float(y[0]))
                                ubhflyp.append(e)
                                ubhflyp[0] += 1
                                    
                            elif "IRC-IRC-IRC-IRC-IRC" in line or "Initial guess <" in line:
                                    if verbose == True:
                                        if zMatrix != []: log_zMatrix(zMatrix, Matrix_number)
                                        if internal_angles != []: log_angles(internal_angles, Matrix_number)
                                        log.write("\n\n")
                                        if internal_angles != []: report_angles(internal_angles, Matrix_number)
                                        if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)
                                        print()

                                    zMatrix = []
                                    internal_angles = []
                                    
                            elif "HF=" in line:
                                x, *y = line.split("HF=")
                                x, *y = y[0].split("\\")
                                try:
                                    tmp=float(x)
                                except(ValueError):
                                    p, *q = x.split(",")
                                    x = p
                                try:
                                    if Energy[0] != None and Energy[1] != None and Energy[2] == 0:
                                        Energy = (float(x), 1, 1)
                                        if verbose == True:
                                            log.write("HF: {}\n".format(Energy[0]))
                                    else:
                                        Energy = (float(x), 1, 0)
                                        if verbose == True:
                                            log.write("HF: {}\n".format(Energy[0]))
                                except(ValueError) as err:
                                    
                                    print("{1} in HF (E)\n{2}\n{3}".format(type(err), err.args, Energy))
                                    if verbose ==  True:
                                        log.write("{1} in HF (E)\nError: {2}\n Tuple: {3}\n".format(type(err), err.args, Energy))
                                except(TypeError) as err:
                                    print("{1} in HF (E)\n{2}\n{3}".format(type(err), err.args, Energy))
                                    if verbose == True:
                                        log.write("{1} in HF (E)\nError: {2}\n Tuple: {3}\n".format(type(err), err.args, Energy))

                                    

                            elif "Normal termination of Gaussian" in line:
                                if verbose == True:
                                    if zMatrix != []: log_zMatrix(zMatrix, Matrix_number)
                                    if internal_angles != []: log_angles(internal_angles, Matrix_number)
                                    log.write("\n\n\n")
                                if internal_angles != []: report_angles(internal_angles, Matrix_number)
                                if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)

                                if verbose == True:
                                    log.write("\n=====================E(UB+HF+LYP)=====================")
                                    print("\n=====================E(UB+HF+LYP)=====================")
                                    log.write("Step                 E")
                                    for item in ubhflyp:
                                        if isinstance(item, UBHFLYP):
                                            log.write("{}\n".format(item))
                                            print("{}".format(item))
                                        else: pass
                                else:
                                    print("\n=====================E(UB+HF+LYP)=====================")
                                    print("{}".format(ubhflyp[len(ubhflyp)-1]))

                                if Energy[1] == 1: print("Hartree-Fock= {}".format(Energy[0]))
                                if Energy[2] == 1: print("-Many HF found, see details in log file")


                if verbose == True: log.close()
                print("Finished")
        except (RuntimeError, TypeError, ValueError, OSError) as inst:
                print ("There was an error on IRC")
                print (type(inst))
                print (inst.args)
def scan():
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

        
    global file
    global verbose
    if verbose == True:
        global log
        log = open("gauss09sAWK.log", "w")
        from datetime import datetime
        log.write("======================Gaussian09 simple AWK======================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
        log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
        log.close()
        log = open("gauss09sAWK.log", "a")
        print("NOTE: Logs are turned on")
        print("NOTE 2: Saving Gauss09 sAWK logs in 'gauss09sAWK.log'") 
    checkfile(file)
    with open(file, "r") as f:
        n = 0
        steps=[]
        points=[]
        ub3lyp=[int(0)]
        Energy = (None, None, 0) #[E, type, found multiple HF values?] Types: 1=Hartree-Fock
        for line in f:
            n += 1
            if "Variable" in line and "Step" in line and "Value" in line and "No. Steps" not in line and "Step-Size" not in line:
                import linecache
                j = 2
                while "--------------------------------------------------" not in linecache.getline(file, n + j) and "A total of" not in linecache.getline(file, n + j):
                    if line != " ":
                        q, *r = linecache.getline(file, n + j).split(" ")
                    k = 0
                    while k < len(r):
                        if r[k] == '': r.remove(r[k])
                        else: k += 1
                    p=Step(int(r[1]),int(r[0]),float(r[2]))
                    steps.append(p)
                    j += 1
            elif "Summary" in line and "potential" in line and "surface" in line and "scan"in line:
                import linecache
                j = 3
                while "----" not in linecache.getline(file, n + j):
                    x, *y = linecache.getline(file, n + j).split(" ")
                    j += 1
                    k = 0
                    while k < len(y):
                        if y[k] == "": y.remove(y[k])
                        else: k += 1
                    p=Point(int(y[0]), float(y[1]), None)
                    points.append(p)
            elif "HF=" in line:
                x, *y = line.split("HF=")
                x, *y = y[0].split("\\")
                try:
                    tmp=float(x)
                except(ValueError):
                    p, *q = x.split(",")
                    x = p
                try:
                    if Energy[0] != None and Energy[1] != None and Energy[2] == 0:
                        Energy = (float(x), 1, 1)
                        if verbose == True:
                            log.write("HF= {}\n".format(Energy[0]))
                    else:
                        Energy = (float(x), 1, 0)
                        if verbose == True:
                            log.write("HF= {}\n".format(Energy[0]))
                            
                except(ValueError, TypeError) as err:
                    print("{1} in HF (E)\nError: {2}\nTuple: {3}".format(type(err), err.args, Energy))
                    if verbose == True:
                        log.write("{1} in HF (E)\nError: {2}\nTuple: {3}".format(type(err), err.args, Energy))
            elif "SCF Done:" in line and "E(UB3LYP)" in line:
                x, *y = line.split("=")
                x, *y = y[0].split(" ")
                for item in y:
                    if item == "": y.remove(item)
                e=UB3LYP(ub3lyp[0] + 1, float(y[0]))
                ub3lyp.append(e)
                ub3lyp[0] += + 1
                
            elif "Normal termination of Gaussian" in line:
                if points != []:
                    for item in ub3lyp:
                        if isinstance(item, UB3LYP):
                            points[item.n-1].energy = item.E
                    if verbose == True: log.write("\n======================POINTS=======================\nPoint            Coordinate                Energy\n")
                    print("\n======================POINTS=======================")
                    print("Point            Coordinate                Energy\n")
                    for item in points:
                        print("{}".format(item))
                        if verbose == True: log.write("{}\n".format(item))
                if verbose == True:
                    log.write("\n=====================E(UB3LYP)=====================\n")
                    log.write("Point              Value\n")
                    print("\n=====================E(UB3LYP)=====================")
                    print("Point              Value")
                    for item in ub3lyp:
                        if isinstance(item, UB3LYP):
                            print("{}".format(item))
                            log.write("{}\n".format(item))
                        else: pass
                        
                if steps != []:
                    if verbose == True:
                        log.write("=======================STEPS=======================\n")
                        log.write("Step     Var           Value\n")
                        print("=======================STEPS=======================")
                        print("Step     Var           Value")
                        for item in steps:
                            print("{}".format(item))
                            if verbose == True: log.write("{}\n".format(item))
                if Energy[1] == 1: print("\nHartree-Fock= {}".format(Energy[0]))
                if Energy[2] == 1: print("-Many HF found, see details in log file")


    print("Finished")
    if verbose == True: log.close()
                
    


def main():
        global file
        
        def go():
                textBox1.delete(1.0, END) 
                print("========================START=======================")
                global verbose
                if int(checkBox1_v.get()) == 1: verbose = True
                else: verbose = False
                
                if int(radioButton_v.get()) == 0: op()
                elif (radioButton_v.get()) == 1: irc()
                else: scan()
        #main window
        root = Tk()
        root.title("Gaussian09 simple AWK")
        root.geometry("430x588")
        root.resizable(0, 0)

        #file browser
        global textBox1_v
        textBox1_v = StringVar()
        textBox1 = Entry(root, textvariable=textBox1_v, width=45, state="disabled").grid(padx=1, pady=5, sticky=NW, columnspan=40)
        textBox1_v.set(file)

        #buttons
        button1 = Button(root, text="...", command=lambda:load_file()).grid(padx=380, pady=2, sticky=NW, row=0, column=0)
        button2 = Button(root, text="Go", command=lambda: go()).grid(padx=380, pady=0, sticky=NW, column=0, row=1, columnspan=4)

        #check boxes
        checkBox1_v = IntVar()
        checkBox1 = Checkbutton(root, text="Verbose", variable=checkBox1_v, onvalue=1, offvalue=0).grid(padx=300, pady=3, sticky=NW, row=1, column=0, columnspan=3)
        
        radioButton_v = IntVar()

        #radio buttons
        radioButton1 = Radiobutton(root, text="SCAN", variable=radioButton_v, value=2).grid(padx=104, pady=3, sticky=NW, row=1, column=0, columnspan=2)
        radioButton1 = Radiobutton(root, text="OPT", variable=radioButton_v, value=0).grid(padx=0, pady=3, sticky=NW, row=1, column=0, columnspan=2)
        radioButton2 = Radiobutton(root, text="IRC", variable=radioButton_v, value=1).grid(padx=52, pady=3, sticky=NW, row=1, column=0, columnspan=2)

        #text boxes
        textBox1 = Text(root, height=34, width=60, wrap='word')
        textBox1.grid(padx=0, pady=0, row=2, column=0, sticky=NW, columnspan=20, rowspan=70)
        sys.stdout = StdoutRedirector(textBox1)
        
        root.mainloop()


if __name__ == '__main__':
    main()
sys.stdout = old_stdout
