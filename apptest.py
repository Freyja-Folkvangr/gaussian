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
                        return "{0.center_number}       {0.atomic_number}     {0.atomic_type}     {0.coordinates}".format(self)
        def report_Standard_orientation(Standard_orientation, Matrix_number):
            if verbose == True: print("================STANDARD ORIENTATION {}================".format(Matrix_number))
            else: print("================STANDARD ORIENTATION================")
            print("NUMBER    Atom    Type        Coords")
            for item in Standard_orientation:
                print(item)
        def log_Standard_orientation(Standard_orientation, Matrix_number):
            if verbose == True:
                    log.write("================STANDARD ORIENTATION {}================\n".format(Matrix_number))
                    log.write("NUMBER    Atom    Type        Coords\n")
                    for item in Standard_orientation:
                            log.write("{}\n".format(item))
        def report_orbitals(aocc, bocc, avirt, bvirt):
                try:
                    global n
                    n = n + 1
                    print("=======================RESULTS======================")
                    if verbose == True: print("Itineration {}:".format(n))
                    if verbose == True:
                        log.write("\n==========================================================\n")
                        log.write("We organice the data, then we log it\n")
                        log.write("and finally we show final results from\n")
                        log.write("itineration {} as apptest output\n".format(n))
                        log.write("==========================================================\n")
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
                    print("mu= {0}   fh= {1}".format(((avirt[0] + aocc[len(aocc) - 1]) / 2), ((avirt[0] - aocc[len(aocc) - 1]) / 2)))
                    print()
                    if bocc != []:
                        print ("Beta results")
                        print("HOMO: {}".format(bocc[len(bocc) - 1]))
                        print("LUMO: {}".format(bvirt[0]))
                        print("mu= {0}   fh= {1}".format(((bvirt[0] + bocc[len(bocc) - 1]) / 2), ((bvirt[0] - bocc[len(bocc) - 1]) / 2)))
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
                        log = open("apptest.log", "w")
                        from datetime import datetime
                        log.write("=============================APPTEST=============================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
                        log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
                        log.close()
                        log = open("apptest.log", "a")
                        print("NOTE: Logs are turned on")
                        print("NOTE 2: Saving APPTEST logs in 'apptest.log'")    
                checkfile(file)
                with open(file, "r") as f:
                    if verbose == True:
                        log.write("==========================================================\n")
                        log.write("We read, take and log data we're interested in\n")
                        log.write("Empty items like '' will be ignored\n")
                        log.write("If there's another itineration, this step will start over\n")
                        log.write("==========================================================\n")
                    aocc = []
                    bocc = []
                    avirt = []
                    bvirt = []
                    Standard_orientation = []
                    Matrix_number = 0
                    Energy = (None, None) #[E, type] Types: 1=Hartree-Fock
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
                                #print("line: {}".format(linecache.getline(file, Line_number + j)))
                                q, *r = linecache.getline(file, Line_number + j).split(' ')
                                q = []
                                for item in r:
                                        if item != '':
                                                q.append(item)
                                coordinates = Coordinates(float(q[3]), float(q[4]), float(q[5]))
                                electron = Electron(float(q[0]), float(q[1]), float(q[2]), coordinates)
                                Standard_orientation.append(electron)
                                j = j + 1
                            
                        elif "HF=" in line:
                                if verbose == True:
                                        if n > 1: log.write("Last itineration should be the optimized one\n")
                                        log.write(line)
                                elif n > 1: print("Last itineration should be the optimized one")
                                else: pass
                                x, *y = line.split('HF=')
                                x, *y = y[0].split('\\')
                                try:
                                        Energy = (float(x), 1)
                                except (ValueError) as err:
                                    print("Error with HF (E)")
                                    print(type(err))
                                    print(err)
                                
                        elif "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad" in line:
                                if verbose == True and (aocc != [] and avirt != []):
                                        report_orbitals(aocc, bocc, avirt, bvirt)
                                if Standard_orientation != [] and verbose == True:
                                        report_Standard_orientation(Standard_orientation, Matrix_number)
                                        log_Standard_orientation(Standard_orientation, Matrix_number)
                                """
                                if aocc != []:print("\naocc\n{}\n".format(aocc))
                                if bocc != []:print("\nbocc\n{}\n".format(bocc))
                                if avirt != []:print("\navirt\n{}\n".format(avirt))
                                if bvirt != []:print("\nbvirt\n{}\n".format(bvirt))
                                """
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
                                """
                                if aocc != []:print("\naocc\n{}\n".format(aocc))
                                if bocc != []:print("\nbocc\n{}\n".format(bocc))
                                if avirt != []:print("\navirt\n{}\n".format(avirt))
                                if bvirt != []:print("\nbvirt\n{}\n".format(bvirt))
                                """
                                Standard_orientation = []
                                aocc = []
                                bocc = []
                                avirt = []
                                bvirt = []
                                
                                
                        else: pass
                if verbose == True:
                        log.write("\n#AO means alpha occ eigenvalues\n")
                        log.write("\n#AV means alpha virt eigenvalues\n")
                        log.write("\n#BO means beta occ eigenvalues\n")
                        log.write("\n#BV means beta virt eigenvalues\n")
                        log.close()
                if Energy[1] == 1:print("Hartree-Fock= {}".format(Energy[0]))
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
                        return "{0.center_number}       {0.atomic_number}     {0.atomic_type}     {0.coordinates}".format(self)
        def report_zMatrix(zMatrix, Matrix_number):
                if verbose == True: print("======================Z-MATRIX {}======================".format(Matrix_number))
                else: print("======================Z-MATRIX======================")
                print("NUMBER    Atom    Type        Coords")
                for item in zMatrix:
                        print(item)
        def log_zMatrix(zMatrix, Matrix_number):
            if verbose == True:
                    log.write("======================Z-MATRIX {}======================\n".format(Matrix_number))
                    log.write("NUMBER    Atom    Type        Coords\n")
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
                        log = open("apptest.log", "w")
                        from datetime import datetime
                        log.write("=============================APPTEST=============================\n#code's author: Giuliano Tognarelli Buono-core\n#{0}\n#Last run on ".format(file))
                        log.write(datetime.now().strftime("%A %d/%m/%Y at %H:%M (dd/mm/yyyy)\n"))
                        log.close()
                        log = open("apptest.log", "a")
                        print("NOTE: Logs are turned on")
                        print("NOTE 2: Saving APPTEST logs in 'apptest.log'")    
                checkfile(file)
                with open(file, "r") as f:
                    internal_angles = []
                    zMatrix = []
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
                                            electron = Electron(float(q[0]), float(q[1]), float(q[2]), coordinates)
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
                                    
                            elif "IRC-IRC-IRC-IRC-IRC" in line:
                                    if verbose == True:
                                        if zMatrix != []: log_zMatrix(zMatrix, Matrix_number)
                                        if internal_angles != []: log_angles(internal_angles, Matrix_number)
                                        log.write("\n\n")
                                        if internal_angles != []: report_angles(internal_angles, Matrix_number)
                                        if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)
                                        print()

                                    zMatrix = []
                                    internal_angles = []
                                    
                            elif "Normal termination of Gaussian" in line:
                                if verbose == True:
                                    if zMatrix != []: log_zMatrix(zMatrix, Matrix_number)
                                    if internal_angles != []: log_angles(internal_angles, Matrix_number)
                                    log.write("\n\n\n")
                                if internal_angles != []: report_angles(internal_angles, Matrix_number)
                                if zMatrix != []: report_zMatrix(zMatrix, Matrix_number)


                if verbose == True:
                        log.write("\n#AO means alpha occ eigenvalues\n")
                        log.write("\n#AV means alpha virt eigenvalues\n")
                        log.write("\n#BO means beta occ eigenvalues\n")
                        log.write("\n#BV means beta virt eigenvalues\n")
                        log.close()
                print("Finished")
        except (RuntimeError, TypeError, ValueError, IndexError) as inst:
                print ("There was an error on IRC")
                print (type(inst))
                print (args)

def main():
        global file
        
        def go():
                textBox1.delete(1.0, END) 
                print("========================START=======================")
                global verbose
                if int(checkBox1_v.get()) == 1: verbose = True
                else: verbose = False
                
                if int(radioButton_v.get()) == 0: op()
                else: irc()
                
        root = Tk()
        root.title("APPTEST")
        root.geometry("430x588")
        root.resizable(0, 0)

        global textBox1_v
        textBox1_v = StringVar()
        textBox1 = Entry(root, textvariable=textBox1_v, width=45, state="disabled").grid(padx=1, pady=5, sticky=NW, columnspan=40)
        textBox1_v.set(file)
        
        button1 = Button(root, text="...", command=lambda:load_file()).grid(padx=380, pady=2, sticky=NW, row=0, column=0)
        button2 = Button(root, text="Go", command=lambda: go()).grid(padx=380, pady=0, sticky=NW, column=0, row=1, columnspan=4)
        
        checkBox1_v = IntVar()
        checkBox1 = Checkbutton(root, text="Verbose", variable=checkBox1_v, onvalue=1, offvalue=0).grid(padx=300, pady=3, sticky=NW, row=1, column=0, columnspan=3)
        
        radioButton_v = IntVar()
        radioButton1 = Radiobutton(root, text="OPT", variable=radioButton_v, value=0).grid(padx=0, pady=3, sticky=NW, row=1, column=0, columnspan=2)
        radioButton2 = Radiobutton(root, text="IRC", variable=radioButton_v, value=1).grid(padx=52, pady=3, sticky=NW, row=1, column=0, columnspan=2)

        textBox1 = Text(root, height=34, width=60, wrap='word')
        textBox1.grid(padx=0, pady=0, row=2, column=0, sticky=NW, columnspan=20, rowspan=70)
        sys.stdout = StdoutRedirector(textBox1)
        
        root.mainloop()


if __name__ == '__main__':
    main()
sys.stdout = old_stdout
