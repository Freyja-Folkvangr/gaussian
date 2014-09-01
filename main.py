__doc__ = info = '''
Created on Mar 26, 2014
@author: giuliano

Tabbed interface script
www.sunjay-varma.com

See the changelog and download updates at:
http://www.github.com/evergreen2/gaussian

Please report any bug you may find in this program,
verify your results and report any issue.
That's much appreciated!
'''

import sys
from classes import *
from tab_constructor import *
from tkinter import *
from elements import *


old_stdout = sys.stdout

global t0
global file
file = "None"

global dualStr
dualStr = 'first'

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

def checkfile(file):
    import os.path
    if os.path.exists(file):
        x, *y = file.split('.')
        if y[len(y)-1] != "log" and y[len(y)-1] != "txt":
                load_file()
        else:
            return True
    else:
        load_file()

def load_file():
    import os.path
    global file
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
    try:
        global textBox2
        if 'ABORTED' in textBox2.get(1.0,END): textBox2.delete(1.0,END)
        if int(checkBox8_v.get()) == 1:
            print("\n   Dual mode counts the data between two files so that:")
            print("-> Select the {} Density file (*.txt)".format(dualStr))
            fname = askopenfilename(filetypes=(("Density Files", "*.txt"),
                                           ("All files", "*.*")))
            file=fname
            if os.path.exists(file):
                textBox1_v.set(file)
            else:
                print("=======================ABORTED======================")
                return False
        else:
            print("\n-> Open Gaussian log/output files (*.log or *.txt)")
            fname = askopenfilename(filetypes=(("*.log output files", "*.log"),
                                               ("*.txt output files", "*.txt"),
                                               ("All files", "*.*")))
            file=fname
            if os.path.exists(file):
                textBox1_v.set(file)
            else:
                print("=======================ABORTED======================")
                return False

    except (RuntimeError, IOError) as inst:
            print ("There was an error while oppening the file.\n {}".format(file))
            print(type(inst))
            print(inst)

def report_zMatrix(zMatrix, Matrix_number):
    if verbose == True:
        print("======================Z-MATRIX {}======================".format(Matrix_number))
        log.write("======================Z-MATRIX {}======================\n".format(Matrix_number))
        log.write("NUMBER    Atom    Type      Coords\n")
    else: print("======================Z-MATRIX======================")
    print("NUMBER    Atom    Type      Coords")
    for item in zMatrix:
        print(item)
        if verbose == True: log.write("{}\n".format(item))

def report_angles(internal_angles, Matrix_number):
    if verbose == True:
        print("==================INTERNAL-ANGLES {}===================".format(Matrix_number))
        log.write("==================INTERNAL-ANGLES {}===================\n".format(Matrix_number))
        log.write("ANGLE                 VALUE\n")
    else: print("==================INTERNAL-ANGLES===================")
    print("ANGLE                 VALUE")
    for item in internal_angles:
        print(item)
        if verbose == True: log.write("{}\n".format(item))



def report_Standard_orientation(Standard_orientation, Matrix_number):
    if verbose == True:
        print("================STANDARD ORIENTATION {}================".format(Matrix_number))
        log.write("================STANDARD ORIENTATION {}================\n".format(Matrix_number))
        log.write("NUMBER    Atom    Type     Coords\n")
    else: print("================STANDARD ORIENTATION================")
    print("NUMBER    Atom    Type     Coords")
    for item in Standard_orientation:
        print(item)
        if verbose == True: log.write("{}\n".format(item))
def report_orbitals(aocc, bocc, avirt, bvirt):
    global orbitals
    global itineration
    itineration += 1
    print("======================HOMO-LUMO=====================")
    if verbose == True:
        log.write("\n======================HOMO-LUMO {}=====================\n".format(itineration))
        print("Itineration {}:".format(itineration))
    if verbose == True:
        log.write("AO:\n")
        if aocc == []: log.write("None")
        else: log.write(str(aocc))
        log.write("\n\nAV:\n")
        if avirt == []: log.write("None")
        else: log.write(str(avirt))
        log.write("\n\nBO:\n")
        if bocc == []: log.write("None")
        else: log.write(str(bocc))
        log.write("\n\nBV:\n")
        if bvirt == []: log.write("None\n\n")
        else: log.write(str(bvirt))
        if bocc != []:
            orbitals.append({"Type":"ALPHA", "HOMO":float(aocc[len(aocc) - 1]), "LUMO":float(avirt[0]), "U":float(((avirt[0] + aocc[len(aocc) - 1]) * 0.5)), "HN":float(((avirt[0] - aocc[len(aocc) - 1]))), "Itineration":int(itineration)})
            orbitals.append({"Type":"BETA", "HOMO":float(bocc[len(bocc) - 1]), "LUMO":float(bvirt[0]), "U":float(((bvirt[0] + bocc[len(bocc) - 1]) * 0.5)), "HN":float(((bvirt[0] - bocc[len(bocc) - 1]))), "Itineration":int(itineration)})
        else:
            orbitals.append({"Type":"ALPHA", "HOMO":float(aocc[len(aocc) - 1]), "LUMO":float(avirt[0]), "U":float(((avirt[0] + aocc[len(aocc) - 1]) * 0.5)), "HN":float(((avirt[0] - aocc[len(aocc) - 1]))), "Itineration":int(itineration)})
    if bocc != []: print ("    |Alpha results\n____|________________")
    print("HOMO|  {}".format(aocc[len(aocc) - 1]))
    print("LUMO|  {}\n____|________________".format(avirt[0]))
    print("μ= {0}  &  Ηη= {1}".format(((avirt[0] + aocc[len(aocc) - 1]) * 0.5), ((avirt[0] - aocc[len(aocc) - 1]))))
    print()
    if bocc != []:
        print ("    |Beta results\n____|________________")
        print("HOMO|  {}".format(bocc[len(bocc) - 1]))
        print("LUMO|  {}\n____|________________".format(bvirt[0]))
        print("μ= {0}  &  Ηη= {1}".format(((bvirt[0] + bocc[len(bocc) - 1]) * 0.5), ((bvirt[0] - bocc[len(bocc) - 1]))))


def to_element(atomic_number):
    for element in ELEMENTS:
        if int(element.number) == int(atomic_number):
            return element.symbol
        else: pass
    return atomic_number

def generate_multi_step(data, type='SP'):
    pass

def go():
    def process():
        global n
        global file
        global verbose
        global orbitals
        global itineration
        path=[{"Maximum points per path":None,"Step size":None, "Last path":None}]
        itineration = 0
        n = 0
        Line_number = 0
        Matrix_number = 0
        orbitals = []
        internal_angles = []
        zMatrix = []
        conflicted_lines=[]
        ubhflyp=[int(0)]
        if checkfile(file) == False:
            print("File Not Found\n=======================ABORTED======================")
            return
        Energy = (None, None, 0) #[E, type, found multiple HF values?] Types: 1=Hartree-Fock
        aocc = []
        bocc = []
        avirt = []
        bvirt = []
        condensed_matrix = []
        Standard_orientation = []
        steps = []
        points = []
        SP_points = []
        ub3lyp = [int(0)]
        rb3lyp = [int(0)]

        for line in f:
            Line_number += 1
            n += 1
            if int(checkBox2_v.get()) == 1:
                #
                if "Condensed to atoms (all electrons):" in line and int(checkBox12_v.get()) == 1:
                    import linecache
                    j = 2
                    while "Mulliken charges:" not in linecache.getline(file, Line_number + j):
                        matrix=[]
                        x, *y = linecache.getline(file, Line_number + j).split(" ")
                        tmp = y[:]
                        for item in tmp:
                            if item != '' and item != ' ':
                                y.append(item)
                        matrix.append(y) #Append the matrix element
                        j += 1
                    condensed_matrix.append(matrix) #Append all matrix elemnts to a list of matrixes
                #homo, lumo, orbitals, etc
                elif "Alpha  occ. eigenvalues" in line:
                    x, *y = line.split('--')
                    args = y[0].split(' ')
                    for item in args:
                        if item != '':
                            try: aocc.append(float(item))
                            except (ValueError) as err:
                                i = 0
                                tmp = ''
                                while i < len(item):
                                    if item[i] == '-':
                                        try:
                                            aocc.append(float(tmp))
                                        except(ValueError): pass
                                        tmp = ''
                                        tmp += item[i]
                                    else:
                                        try:
                                            if int(item[i]) >= 0 or item[i] == '.':
                                                tmp += item[i]
                                        except(ValueError):
                                            if item[i] == '.':
                                                tmp += item[i]
                                            else:
                                                print("Unknown error! {}".format(line))
                                    i += 1
                                aocc.append(float(tmp))
                elif "Alpha virt. eigenvalues" in line:
                    x, *y = line.split('--')
                    args = y[0].split(' ')
                    for item in args:
                        if item != '':
                            try: avirt.append(float(item))
                            except (ValueError) as err:
                                i = 0
                                tmp = ''
                                while i < len(item):
                                    if item[i] == '-':
                                        try:
                                            avirt.append(float(tmp))
                                        except(ValueError): pass
                                        tmp = ''
                                        tmp += item[i]
                                    else:
                                        try:
                                            if int(item[i]) >= 0 or item[i] == '.':
                                                tmp += item[i]
                                        except(ValueError):
                                            if item[i] == '.':
                                                tmp += item[i]
                                            else:
                                                print("Unknown error! {}".format(line))
                                    i += 1
                                avirt.append(float(tmp))
                elif "Beta  occ. eigenvalues" in line:
                    x, *y = line.split('--')
                    args = y[0].split(' ')
                    for item in args:
                        if item != '':
                            try: bocc.append(float(item))
                            except(ValueError) as err:
                                i = 0
                                tmp = ''
                                while i < len(item):
                                    if item[i] == '-':
                                        try:
                                            bocc.append(float(tmp))
                                        except(ValueError): pass
                                        tmp = ''
                                        tmp += item[i]
                                    else:
                                        try:
                                            if int(item[i]) >= 0 or item[i] == '.':
                                                tmp += item[i]
                                        except(ValueError):
                                            if item[i] == '.':
                                                tmp += item[i]
                                            else:
                                                print("Unknown error! {}".format(line))
                                    i += 1
                                bocc.append(float(tmp))
                elif "Beta virt. eigenvalues" in line:
                    x, *y = line.split('--')
                    args = y[0].split(' ')
                    for item in args:
                        if item != '':
                            try: bvirt.append(float(item))
                            except (ValueError) as err:
                                i = 0
                                tmp = ''
                                while i < len(item):
                                    if item[i] == '-':
                                        try:
                                            bvirt.append(float(tmp))
                                        except(ValueError): pass
                                        tmp = ''
                                        tmp += item[i]
                                    else:
                                        try:
                                            if int(item[i]) >= 0 or item[i] == '.':
                                                tmp += item[i]
                                        except(ValueError):
                                            if item[i] == '.':
                                                tmp += item[i]
                                            else:
                                                print("Unknown error! {}".format(line))
                                    i += 1
                                bvirt.append(float(tmp))
            if int(checkBox3_v.get()) == 1:
                if "Standard orientation" in line:
                    import linecache
                    j = 5
                    while "-----" not in linecache.getline(file, Line_number + j):
                        q, *r = linecache.getline(file, Line_number + j).split(' ')
                        q = []
                        for item in r:
                            if item != '': q.append(item)
                        coordinates = Coordinates(float(q[3]), float(q[4]), float(q[5]))
                        electron = Electron(int(q[0]), int(q[1]), int(q[2]), coordinates)
                        Standard_orientation.append(electron)
                        j = j + 1
            if int(checkBox4_v.get()) == 1:
                if "Z-Matrix orientation" in line:
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
            if int(checkBox5_v.get()) == 1:
                if "Interatomic angles" in line:
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
            if int(checkBox6_v.get()) == 1:
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
            if int(checkBox7_v.get()) == 1:
                if "Summary" in line and "potential" in line and "surface" in line and "scan"in line:
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
            if int(checkBox9_v.get()) == 1:
                if "INPUT DATA FOR" in line:
                    #Getting the input configuration
                    import linecache
                    j = 2
                    while path[0]["Maximum points per path"] == None and path[0]["Step size"] == None:
                        if "Maximum points per path" in linecache.getline(file, n + j):
                            x, *y = linecache.getline(file, n + j).split("=")
                            path[0]["Maximum points per path"] = int(y[0])
                        elif "Step size" in linecache.getline(file, n + j):
                            x, *y = linecache.getline(file, n + j).split("=")
                            path[0]["Step size"] = int(y[0])
                        j += 1
                elif "Cartesian Forces" in line:
                    tmp = None
                    import linecache
                    #Getting the point number and path number
                    j = 2
                    while True:
                        if "Point Number:" in linecache.getline(file, n + j) and "Path Number:" in linecache.getline(file, n + j):
                            tmp=''
                            for item in linecache.getline(file, n + j):
                                if item != '' and item != ' ' and item != ':':
                                    tmp += item
                            x, *y = tmp.split('PointNumber')
                            x, *y = y[0].split('PathNumber')
                            SP_point_number = int(x)
                            SP_path_number = int(y[0])
                            break
                        elif "Input orientation" in linecache.getline(file, n + j) or "IRC-IRC-IRC-IRC-IRC" in linecache.getline(file, n + j):
                            tmp = 'skip'
                            if verbose == True:
                                log.write('Warning: Incomplete point! Line: {}\n'.format(n))
                                print('Warning: Incomplete point! Line: {}'.format(n))
                                break
                            else:
                                print('Warning: Incomplete point! Line: {}'.format(n))
                                break
                        j += 1

                    if tmp != 'skip':
                        #Getting the input matrix
                        input_matrix = []
                        j = 0
                        while "Input orientation:" not in linecache.getline(file, n + j):
                            j -= 1
                        j += 5
                        while "--------" not in linecache.getline(file, n + j):
                            x, *y = linecache.getline(file, n + j).split(" ")
                            y2 = y
                            y = []
                            for item in y2:
                                if item != '' and item != ' ':
                                    y.append(item)
                            coords = Coordinates(y[3].rstrip('\n'), y[4].rstrip('\n'), y[5].rstrip('\n'))
                            electron = Electron(y[0], y[1], y[2], coords)
                            input_matrix.append(electron)
                            j += 1
                        #Getting the output matrix
                        output_matrix = []
                        j = -2
                        while "----" not in linecache.getline(file, n + j):
                            x, *y = linecache.getline(file, n + j).split(" ")
                            tmp = y
                            y = []
                            #Cleaning up the spaces
                            for item in tmp:
                                if item != '' and item != ' ': y.append(item)
                            coords = Coordinates(float(y[2]), float(y[3]), float(y[4]))
                            electron = Electron(y[0], y[1], None, coords)
                            output_matrix.append(electron)
                            j -= 1
                        if path[0]["Last path"] == SP_path_number:
                            SP_points.append({"Input":input_matrix, "Output":{"Matrix":output_matrix, "Point number":int(SP_point_number), "Path number":int(SP_path_number)}})
                        elif path[0]["Last path"] == None:
                            path[0]["Last path"] = int(SP_path_number)
                            SP_points.append({"Input":input_matrix, "Output":{"Matrix":output_matrix, "Point number":int(SP_point_number), "Path number":int(SP_path_number)}})
                        elif path[0]["Last path"] != SP_path_number:
                            path.append(SP_points)
                            SP_points = []
                            SP_points.append({"Input":input_matrix, "Output":{"Matrix":output_matrix, "Point number":int(SP_point_number), "Path number":int(SP_path_number)}})
                            path[0]["Last path"] = int(SP_path_number)
                    else:
                        tmp2 = None
            if "HF=" in line:
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

            elif ("GradGradGradGradGradGradGrad" in line or "Initial guess <" in line) and (int(checkBox2_v.get()) == 1 or int(checkBox3_v.get()) == 1):
                if verbose == True and (aocc != [] and avirt != []):
                    bocc.sort()
                    bvirt.sort()
                    aocc.sort()
                    avirt.sort()
                    bocc.reverse()
                    bvirt.reverse()
                    report_orbitals(aocc, bocc, avirt, bvirt)
                if Standard_orientation != [] and verbose == True:
                    report_Standard_orientation(Standard_orientation, Matrix_number)
                Standard_orientation = []
                aocc = []
                bocc = []
                avirt = []
                bvirt = []
            elif ("IRC-IRC-IRC-IRC-IRC" in line or "Initial guess <" in line) and ((int(checkBox4_v.get()) == 1) or (int(checkBox5_v.get()) == 1)):
                if internal_angles != [] or zMatrix != []: Matrix_number += 1
                if verbose == True:
                    if internal_angles != []: report_angles(internal_angles, Matrix_number)
                    if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)
                    print()
                zMatrix = []
                internal_angles = []
            elif "SCF Done:" in line and ("E(UB+HF-LYP)" in line or "E(UB3LYP)" in line or "E(RB3LYP)"):
                if "E(UB+HF-LYP)" in line:
                    x, *y = line.split("=")
                    x, *y = y[0].split(" ")
                    tmp = y
                    y = []
                    for item in tmp:
                        if item != '' and item != ' ': y.append(item)
                    e=UBHFLYP(ubhflyp[0] +1, float(y[0]))
                    ubhflyp.append(e)
                    ubhflyp[0] += 1
                elif "E(UB3LYP)" in line:
                    x, *y = line.split("=")
                    x, *y = y[0].split(" ")
                    for item in y:
                        if item == "": y.remove(item)
                    e=UB3LYP(ub3lyp[0] + 1, float(y[0]))
                    ub3lyp.append(e)
                    ub3lyp[0] += + 1
                elif "E(RB3LYP)" in line:
                    x, *y = line.split("=")
                    x, *y = y[0].split(" ")
                    for item in y:
                        if item == "": y.remove(item)
                    e=RB3LYP(rb3lyp[0] + 1, float(y[0]))
                    rb3lyp.append(e)
                    rb3lyp[0] += + 1
                else: pass

            if "Normal termination of Gaussian" in line and (int(checkBox2_v.get()) == 1 or
                                                                int(checkBox3_v.get()) == 1 or
                                                                int(checkBox4_v.get()) == 1 or
                                                               int(checkBox6_v.get()) == 1 or
                                                               int(checkBox7_v.get()) == 1 or
                                                               int(checkBox9_v.get()) == 1):

                print("matrix: {}".format(len(condensed_matrix)))
                print(condensed_matrix)
                if int(checkBox9_v.get()) == 1:
                    path.append(SP_points)
                    link_status = False
                    multi_step = open("results.com", "w")
                    multi_step.close()
                    multi_step = open("results.com", "a")
                    tmp = -1
                    for path_number in path:
                        if isinstance(path_number, list):
                            for point in path_number:
                                tmp += 1
                                if verbose == True:
                                    if link_status == True: multi_step.write("\n--Link1--\n")
                                    #multi_step.write("%chk=irc.chk\n")
                                    multi_step.write("%mem=6GB\n")
                                    multi_step.write("# b3lyp/6-31g(d) nosymm scf=qc test\n\n")
                                    multi_step.write("SP {}\n\n".format(tmp))
                                    multi_step.write("0 1\n")
                                    link_status = True
                                    #log.write("\n==============MATRIX FOR POINT {} AND PATH {}=============\n".format(point["Output"]["Point number"], point["Output"]["Path number"]))
                                for matrix_element in point["Input"]:
                                    if verbose == True:
                                        multi_step.write("   {}    {}\n".format(to_element(matrix_element.atomic_number),matrix_element.coordinates))
                    SP_points = []
                    if verbose == True:
                        tmp=''
                        log.write("=======================DATA=======================\n")
                        log.write("Type: ")
                        for item in path:
                            log.write("{} ".format(type(item)))
                            if len(tmp) > 0:
                                tmp += ("-{}".format(len(item)))
                            else:
                                tmp += ("{}".format(len(item)))
                        log.write("\nArgs: {}\n\n".format(tmp))
                        log.write("==================================================\n")
                        log.write("{}\n".format(path))
                if int(checkBox6_v.get()) == 1:
                    if steps != []:
                        if verbose == True:
                            log.write("=======================STEPS=======================\n")
                            log.write("Step     Var           Value\n")
                            print("=======================STEPS=======================")
                            print("Step     Var           Value")
                            for item in steps:
                                print("{}".format(item))
                                if verbose == True: log.write("{}\n".format(item))
                if int(checkBox7_v.get()) == 1:
                    if points != [] and ub3lyp != [0]:
                        for item in ub3lyp:
                            if isinstance(item, UB3LYP):
                                points[item.n-1].energy = item.E
                        if verbose == True: log.write("\n======================POINTS=======================\nPoint            Coordinate                Energy\n")
                        print("\n======================POINTS=======================")
                        print("Point            Coordinate                Energy\n")
                        for item in points:
                            print("{}".format(item))
                            if verbose == True: log.write("{}\n".format(item))
                if (int(checkBox6_v.get()) == 1 or int(checkBox7_v.get()) == 1) and verbose == True and ub3lyp != [0]:
                    log.write("\n=====================E(UB3LYP)=====================\n")
                    log.write("Point              Value\n")
                    print("\n=====================E(UB3LYP)=====================")
                    print("Point              Value")
                    for item in ub3lyp:
                        if isinstance(item, UB3LYP):
                            print("{}".format(item))
                            log.write("{}\n".format(item))
                        else: pass
                if aocc != [] and avirt != []:
                    bocc.sort()
                    bvirt.sort()
                    aocc.sort()
                    avirt.sort()
                    bocc.reverse()
                    bvirt.reverse()
                    report_orbitals(aocc, bocc, avirt, bvirt)
                if Standard_orientation != []:
                    report_Standard_orientation(Standard_orientation, Matrix_number)
                Standard_orientation = []
                aocc = []
                bocc = []
                avirt = []
                bvirt = []
                if int(checkBox4_v.get()) == 1 or int(checkBox5_v.get()) == 1:
                    if internal_angles != [] or zMatrix != []: Matrix_number += 1
                    if internal_angles != []: report_angles(internal_angles, Matrix_number)
                    if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)
                    if verbose == True and ubhflyp != [0]:
                        log.write("\n=====================E(UB+HF+LYP)=====================")
                        print("\n=====================E(UB+HF+LYP)=====================")
                        print("Step                 E")
                        log.write("Step                 E")
                        for item in ubhflyp:
                            if isinstance(item, UBHFLYP):
                                log.write("{}\n".format(item))
                                print("{}".format(item))
                            else: pass
                    elif verbose == False and ubhflyp != [0]:
                        print("\n=====================E(UB+HF+LYP)=====================")
                        print("Step                 E")
                        print("{}".format(ubhflyp[len(ubhflyp)-1]))
                    else:
                        pass

            #if Energy[1] == 1: print("Hartree-Fock= {}".format(Energy[0]))
            #if Energy[2] == 1: print("-Many HF found, see details in log file")

        if int(checkBox2_v.get()) == 1 and verbose == True:

            if rb3lyp != []:
                log.write("\n=====================E(RB3LYP)=====================\n")
                log.write("     Point                             Value\n")
                print("\n=====================E(RB3LYP)=====================")
                print("     Point                              Value")
                for item in rb3lyp:
                    if isinstance(item, RB3LYP):
                        print("     {}".format(item))
                        log.write("     {}\n".format(item))
                    else: pass

            log.write("=========================================================================\n=============================HOMO-LUMO SUMMARY===========================\n\n")
            if int(checkBox11_v.get()) == 1:
                log.write("No   Type     HOMO         LUMO         u                      Hn                    E(RB3LYP)\n")
            else:
                log.write("//No,HOMO,LUMO,u,Hn,E(RB3LYP);\n")
            i = 0
            for item in orbitals:
                log.write("{}".format(item["Itineration"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range(0, 6 - len(str(item["Itineration"]))-1):
                        log.write(" ")
                    log.write("{}".format(item["Type"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range (0, 8 - len(item["Type"]) - 1):
                        log.write(" ")
                else:
                    log.write(",")
                log.write("{}".format(item["HOMO"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range(0, 13 - int(len(str(item["HOMO"]))) - 1):
                        log.write(" ")
                else:
                    log.write(",")
                log.write("{}".format(item["LUMO"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range(0, 13 - len(str(item["LUMO"])) - 1):
                        log.write(" ")
                else:
                    log.write(",")
                log.write("{}".format(item["U"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range(0, 25 - len(str(item["U"])) - 1):
                        log.write(" ")
                else:
                    log.write(",")
                log.write("{}".format(item["HN"]))

                if int(checkBox11_v.get()) == 1:
                    for i in range(0, 22 - len(str(item["HN"]))):
                        log.write(" ")
                    log.write("{}\n".format(rb3lyp[item["Itineration"]].E))
                else:
                    log.write(",")
                    log.write("{};".format(rb3lyp[item["Itineration"]].E))
                i += 1

        if conflicted_lines != []:
            print("========================ERROR=======================")
            if verbose == True:
                print("Errors found in these lines:\n{}\n\nSee the log file to see more details.".format(conflicted_lines))
                log.write("Conflictive lines: {}\n".format(conflicted_lines))
            else:
                print("Errors found in these lines:\n{}\n\nTry turning on verbose mode to log the details.".format(conflicted_lines))

    try:
        with open(file, "r", encoding="latin-1") as f:
            process()
    except(FileNotFoundError):
        print("File Not Found\n=======================ABORTED======================")
    except(UnicodeDecodeError) as err:
        print("========================ERROR=======================")
        print(type(err))
        print(err)
        if verbose == True:
            log.write("{}\n\n".format(type(err)))
            log.write("{}\n\n".format(err))
            log.write("{}\n\n".format(err.args))
            print("Error details saved into the log file.")
        print("=======================ABORTED======================")

def dualm(index):
    #verbose_anything = False
    if verbose == True:
        log.write("\n\n{}\n\n".format(file))
    if int(checkBox8_v.get()) == 0:
        return
    global n
    global dual
    n = 0
    count = 0
    try:
        with open(file, "r") as f:
            for line in f:
                if " \n" == line: count = 1
                elif count == 1: count += 1
                elif count == 2: count += 1
                if "augmentation occupancies" not in line and count == 3:
                    x, *y = line.split(" ")
                    for item in y:
                        if item != ' ' and item != '':
                            dual[index] += 1
                    #if verbose_anything == True:
                    #    log.write("{}\n".format(y))
                elif "augmentation occupancies" in line and verbose == True:
                    #verbose_anything = True
                    log.write("{} - line {}\n".format(dual,n))
                    log.write("{}\n".format(line))
                else:
                    pass
                n+=1
        if verbose == True:
            log.write("\n\nTask completed with code {}\n\n".format(dual))
    except(FileNotFoundError):
        print("File Not Found\n=======================ABORTED======================")


def main():
    # ======= TOGGLE EVENTS ========
    def checkBox9_onClickEvent(event):
        checkBox2_v.set(0)
        checkBox3_v.set(0)
        checkBox4_v.set(0)
        checkBox5_v.set(0)
        checkBox6_v.set(0)
        checkBox7_v.set(0)
        checkBox8_v.set(0)
        checkBox2.set(state=DISABLED)

    def write(x): print (x)
    global file
    global textBox2
    def initialize():
        textBox2.delete(1.0, END)
        print("========================START=======================")
        global verbose
        global checkBox1_v
        if int(checkBox1_v.get()) == 1:
            verbose = True
            global log
            log = open("verbose.txt", "w")
            from datetime import datetime
            log.write("=============================Gauss09 sAWK=============================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
            log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
            log.close()
            log = open("verbose.txt", "a")
            print("NOTE: Logs are turned on")
            print("NOTE 2: Saving Gauss09 sAWK logs in 'verbose.txt'")
            if int(checkBox9_v.get()) == 1: print("NOTE 3: Saving multi-step file as 'results.com'")
        else: verbose = False
        if int(checkBox8_v.get()) == 1 and (
            int(checkBox2_v.get()) == 1 or
            int(checkBox3_v.get()) == 1 or
            int(checkBox4_v.get()) == 1 or
            int(checkBox5_v.get()) == 1 or
            int(checkBox6_v.get()) == 1 or
            int(checkBox7_v.get()) == 1):
            print("========================ERROR=======================")
            print("If you want to enable Dual mode in Settings tab, you have to unselect all other options.\nJust Dual and Verbose are allowed to work together\n=======================ABORTED======================")
            return False
        if checkfile(file) == False:
            print("File Not Found\n=======================ABORTED======================")
            return
        import time
        t0 = time.time()
        if int(checkBox8_v.get()) == 1:
            global dual
            dual = [0,0]
            dualm(0)
            dualStr = 'second'
            temp = file
            while temp == file:
                load_file()
            dualm(1)
            print("1º file has | {}\n2º file has | {}\n------------|\nDifference: | {}\n\n========================FINISHED======================".format(dual[0], dual[1], dual[0]-dual[1]))
            if verbose == True:
                log.write("\nFinal result: {}\nDifference: {}\n".format(dual, dual[0]-dual[1]))
            dual[0] == 0
        else:
            go()
            print("========================FINISHED======================")
        print ("Lapsed time: {0:.2f}s".format(time.time() - t0))
        if verbose == True: log.close()
        dualStr = 'first'
        return True
    # ======== MAIN WINDOW =========
    root = Tk()
    root.title("Gaussian 09 simple AWK (BETA 3)")
    root.resizable(0, 0)

    bar = TabBar(root, "Info")

    # ======== CONSOLE TAB =========
    tab1 = Tab(root, "Console")                # notice how this one's master is the root instead of the bar
    #FILE BROWSER
    global textBox1_v
    textBox1_v = StringVar()
    textBox1 = Entry(tab1, textvariable=textBox1_v, width=55, state="disabled").grid(padx=1, pady=5, sticky=NW, columnspan=40)
    textBox1_v.set(file)

    button1 = Button(tab1, text="Open...", command=lambda:load_file()).grid(pady=2, sticky=NW, row=0, column=41)
    button2 = Button(tab1, text="Execute", command=lambda: initialize()).grid(pady=0, sticky=NW, column=41, row=1)

    #button3.config(image=button3Image)
    #button3.image = button3Image
    #button3.configure(image=button3Image)

    # ========== CONSOLE ===========
    global textBox2
    textBox2_v = StringVar()
    textBox2 = Text(tab1, height=34, width=75, wrap='word', bd=1, relief=RIDGE, highlightthickness=0)
    textBox2.grid(padx=0, pady=0, row=2, column=0, sticky=NW, columnspan=50, rowspan=50)
    sys.stdout = StdoutRedirector(textBox2)

    #textBox2 = Text(tab1, width=63, height=1, bd=1, relief=RIDGE, highlightthickness=0)
    #textBox2.focus()
    #textBox2.grid(padx=0, pady=0, row=71, column=0, sticky=NW, columnspan=50, rowspan=34)
    #Button(tab1, text="Send", command=(lambda: write(textBox2.get('1.0', END).strip()))).grid(pady=0, sticky=NW, column=41, row=71)



    # ======== SETTINGS TAB =========
    tab2 = Tab(root, "Settings")

    global checkBox1_v
    checkBox1_v = IntVar()
    checkBox1_v.set(1)
    checkBox1 = Checkbutton(tab2, text="Verbose mode", variable=checkBox1_v, onvalue=1, offvalue=0, state = DISABLED).grid(padx=0, pady=0, sticky=NW, row=7, column=0, columnspan=1)

    global checkBox11_v
    checkBox11_v = IntVar()
    checkBox11_v.set(1)
    checkBox11 = Checkbutton(tab2, text="Report HOMO-LUMO summary as table", variable=checkBox11_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=8, column=0, columnspan=1)

    global checkBox12_v
    checkBox12_v = IntVar()
    checkBox12_v.set(1)
    checkBox12 = Checkbutton(tab2, text="Get 'Condensed to atoms' matrix", variable=checkBox12_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=9, column=0, columnspan=1)

    global checkBox10_v
    checkBox10_v = IntVar()
    checkBox10 = Checkbutton(tab2, text="Generate Multi-step", variable=checkBox10_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=6, column=1, columnspan=1)

    global checkBox2_v
    checkBox2_v = IntVar()
    checkBox2 = Checkbutton(tab2, text="HOMO and LUMO", variable=checkBox2_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=1, column=0, columnspan=1)

    global checkBox3_v
    checkBox3_v = IntVar()
    checkBox3 = Checkbutton(tab2, text="Standard orientation matrix", variable=checkBox3_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=1, column=1, columnspan=1)

    global checkBox4_v
    checkBox4_v = IntVar()
    checkBox4 = Checkbutton(tab2, text="zMatrix", variable=checkBox4_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=2, column=0, columnspan=1)

    global checkBox5_v
    checkBox5_v = IntVar()
    checkBox5 = Checkbutton(tab2, text="Interatomic Angles", variable=checkBox5_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=2, column=1, columnspan=1)

    global checkBox6_v
    checkBox6_v = IntVar()
    checkBox6 = Checkbutton(tab2, text="Steps", variable=checkBox6_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=3, column=1, columnspan=1)

    global checkBox7_v
    checkBox7_v = IntVar()
    checkBox7 = Checkbutton(tab2, text="Points", variable=checkBox7_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=3, column=0, columnspan=1)

    global checkBox8_v
    checkBox8_v = IntVar()
    checkBox8 = Checkbutton(tab2, text="Dual", state=DISABLED, variable=checkBox8_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=4, column=0, columnspan=1)

    global checkBox9_v
    checkBox9_v = IntVar()
    checkBox9 = Checkbutton(tab2, text="SP: IRC Input Orientation", variable=checkBox9_v, onvalue=1, offvalue=0).grid(padx=0, pady=0, sticky=NW, row=4, column=1, columnspan=1)


    tab3 = Tab(root, "Info")
    Label(tab3, bg='white', text="BETA version, report bugs in\nour GitHub page below.\n"+info).pack(side=LEFT, expand=YES, fill=BOTH)

    bar.add(tab1)                   # add the tabs to the tab bar
    bar.add(tab2)
    bar.add(tab3)


    bar.config(bd=1, relief=RIDGE)            # add some border

    bar.show()

    root.mainloop()

if __name__ == '__main__':
    main()


sys.stdout = old_stdout
