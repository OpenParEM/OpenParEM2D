# -*- coding: utf-8 -*-
import FreeCAD,Draft
import PySide
from PySide import QtGui ,QtCore
from PySide.QtGui import *
from PySide.QtCore import *

def sequential_mode_check (modeList, checkVariable):

    # count the number of modes with type checkVariable
    count = 0
    i=0
    while (i < len(modeList)):
        if (modeList[i].type == checkVariable):
            count = count +1
        i = i+1

    # skip if there are none of this type
    if (count == 0):
        return False

    # check for too many
    i=0
    while (i < len(modeList)):
        if (modeList[i].type == checkVariable):
            if (modeList[i].name > count):
               App.Console.PrintError("ERROR: Invalid mode/line number of " + str(modeList[i].name) + " with type " + str(modeList[i].type) + ".\n")
               return True;
        i = i + 1

    # check for existence of each mode
    modeNumber=1
    while (modeNumber <= count):
       foundVariable = False
       i=0
       while (i < len(modeList)):
           if (modeList[i].type == checkVariable):
               if (modeList[i].name == modeNumber):
                   if (foundVariable):
                       App.Console.PrintError("ERROR: Duplicate " + str(checkVariable) + " definition found for mode/line number " + str(modeNumber) + ".\n")
                       return True
                   foundVariable = True
           i = i+1

       if (not foundVariable):
           App.Console.PrintError("ERROR: Missing mode/line number " + str(modeNumber) + " for type " + str(checkVariable) + ".\n")
           return True

       modeNumber = modeNumber + 1


    return False

def get_commandSymbol (label,name):
    if name[0:4] == "Text":
        return label[0:2]
    return "null"

def get_commandName (label,showMessage):
    splitText=label[2:].split('(')
    if (len(splitText) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"
    if (splitText[0] == ''):
        if (showMessage):
            App.Console.PrintError("ERROR: \"" + label + "\" is missing a name.\n")
        return "null"
    return splitText[0]

def get_optionCount (label):
    splitText1=label.split('(')
    if (len(splitText1) != 2):
        return 0

    splitText2=splitText1[1].split(')')
    if (len(splitText2) != 2):
        return 0

    splitText3=splitText2[0].split(',')
    return len(splitText3) 

def get_commandOption (label,index):
    splitText1=label.split('(')
    if (len(splitText1) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText2=splitText1[1].split(')')
    if (len(splitText2) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText3=splitText2[0].split(',')
    if (len(splitText3) < index):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    return splitText3[index-1]

def get_paths (label):
    splitText1=label[2:].split('{')
    if (len(splitText1) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText2=splitText1[1].split('}')
    if (len(splitText2) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    paths=splitText2[0].split(',')
    if (len(paths) == 0):
        App.Console.PrintError("ERROR: \"" + label + "\" is missing paths.\n")
        return "null"

    return paths


class Path:
    def __init__(self, name, closed):
        self.name = name 
        self.closed = closed
        self.x = []
        self.y = []
        self.z = []

    def add_point(self,x,y,z):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

    def print(self,file):
        file.write("Path\n")
        file.write("   name=" + self.name + "\n")

        i=0
        while (i < len(self.x)):
            file.write("   point=(" + str(self.x[i]) + "," + str(self.y[i]) + ")" + "\n")
            i=i+1

        file.write("   closed=" + self.closed + "\n")
        file.write("EndPath\n")
        file.write("\n")

class Boundary:
    def __init__(self, name, type, material):
        self.name = name
        self.type = type
        self.material = material
        self.paths = []
        self.directions = []

    def add_path(self,direction,path):
        self.directions.append(direction);
        self.paths.append(path)

    def print(self,file):
        file.write("Boundary\n")
        file.write("   name=" + self.name + "\n")
        file.write("   type=" + self.type + "\n")
        if self.type == "surface_impedance":
            file.write("   material=" + self.material + "\n")

        i=0
        while (i < len(self.paths)):
            if i == 0:
                file.write("   path=" + self.directions[0] + self.paths[0] + "\n")
            else:
                file.write("   path" + self.directions[i] + "=" + self.paths[i] + "\n")
            i=i+1

        file.write("EndBoundary\n")
        file.write("\n")

    def check(self,pathList):
        fail = False

        for path in self.paths:
            found = False
            for testpath in pathList:
                if path == testpath.name:
                    found = True
            if not found:
                App.Console.PrintError("ERROR: Boundary " + self.name + " calls for the undefined path \"" + path + "\".\n")
                fail = True

        if self.type != "surface_impedance" and self.type != "perfect_electric_conductor" and self.type != "perfect_magnetic_conductor" and self.type != "radiation":
            App.Console.PrintError("ERROR: Boundary " + self.name + " calls for an invalid type [must be \"SI\", \"PEC\", \"PMC\", or \"radiation\"].\n")
            fail = True

        return fail


class Mode:
    def __init__(self, name, type, scale):
        self.name = name
        self.type = type
        self.scale = scale
        self.paths = []
        self.directions = []
        self.category = "modal"

    def add_path(self,direction,path):
        self.directions.append(direction);
        self.paths.append(path)

    def set_category(self,category):
        self.category=category

    def print(self,file):
        if self.category == "modal":
            file.write("Mode\n")
            file.write("   mode=" + str(self.name) + "\n")
        if self.category == "line":
            file.write("Line\n")
            file.write("   line=" + str(self.name) + "\n")

        file.write("   type=" + self.type + "\n")
        if (self.scale != ""):
            file.write("   scale=" + self.scale + "\n")

        i=0
        while (i < len(self.paths)):
            if i == 0:
                file.write("   path=" + self.directions[0] + self.paths[0] + "\n")
            else:
                file.write("   path" + self.directions[i] + "=" + self.paths[i] + "\n")
            i=i+1

        if self.category == "modal":
            file.write("EndMode\n\n")
        if self.category == "line":
            file.write("EndLine\n\n")

    def check(self,pathList):
        fail = False

        if self.name < 1:
            App.Console.PrintError("ERROR: Mode \"" + str(self.name) + " uses an invalid mode number [must be > 0].\n")
            fail = True

        for path in self.paths:
            found = False
            for testpath in pathList:
                if path == testpath.name:
                    found = True
            if not found:
                App.Console.PrintError("ERROR: Mode \"" + str(self.name) + " calls for the undefined path \"" + path + "\".\n")
                fail = True

        if self.type != "voltage" and self.type != "current":
            App.Console.PrintError("ERROR: Mode \"" + str(self.name) + " calls for an invalid type [must be \"voltage\" or \"current\"].\n")
            fail = True

        return fail

pathList = []
boundaryList = []
modeList = []

doc = FreeCAD.ActiveDocument
objs = FreeCAD.ActiveDocument.Objects

file = doc.Name+"_modes.txt"
path = "./" + file

try:
    SaveName = QFileDialog.getSaveFileName(None,QString.fromLocal8Bit("Save the modes file"),path,"*.txt") # PyQt4
except Exception:
    SaveName, Filter = PySide.QtGui.QFileDialog.getSaveFileName(None, "Save the modes file", path,"*.txt") # PySide

if SaveName == "":
    App.Console.PrintMessage("Process aborted.\n")
else:
    try:
        fail = False

        file = open(SaveName, 'w')
        try:

            for obj in objs:

                recognizedCommand = False

                #if (not obj.ViewObject.Visibility):
                #    continue

                name = obj.Name
                label = obj.Label

                commandSymbol=get_commandSymbol(label,name)

                if (label[0:2] == "_P"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    nameText=label[2:]
                    if (nameText == ''):
                       App.Console.PrintError("ERROR: \"" + label +"\" does not include a name.\n")
                       fail = True
                       continue

                    if (name[0:4] == "Line"):
                       newPath=Path(nameText,'false')
                       newPath.add_point(obj.Start.x, obj.Start.y, obj.Start.z)
                       newPath.add_point(obj.End.x, obj.End.y, obj.End.z)
                       pathList.append(newPath)

                    if (name[0:9] == "Rectangle"):
                       newPath=Path(nameText,'true')
 
                       x1 = obj.Placement.Base.x
                       y1 = obj.Placement.Base.y
                       z1 = obj.Placement.Base.z
                       newPath.add_point(x1,y1,z1)

                       rotationMatrix=obj.Placement.Rotation.toMatrix();

                       x2 = x1 + rotationMatrix.A11*obj.Length.Value
                       y2 = y1 + rotationMatrix.A21*obj.Length.Value
                       z2 = z1 + rotationMatrix.A31*obj.Length.Value
                       newPath.add_point(x2,y2,z2)

                       x3 = x1 + rotationMatrix.A11*obj.Length.Value + rotationMatrix.A12*obj.Height.Value
                       y3 = y1 + rotationMatrix.A21*obj.Length.Value + rotationMatrix.A22*obj.Height.Value
                       z3 = z1 + rotationMatrix.A31*obj.Length.Value + rotationMatrix.A32*obj.Height.Value
                       newPath.add_point(x3,y3,z3)

                       x4 = x1 + rotationMatrix.A12*obj.Height.Value
                       y4 = y1 + rotationMatrix.A22*obj.Height.Value
                       z4 = z1 + rotationMatrix.A32*obj.Height.Value
                       newPath.add_point(x4,y4,z4)

                       pathList.append(newPath)

                    if (name[0:4] == "Wire"):
                       closed="false"
                       if (obj.Closed):
                           closed="true"

                       x1 = obj.Placement.Base.x
                       y1 = obj.Placement.Base.y
                       z1 = obj.Placement.Base.z
                       rotationMatrix=obj.Placement.Rotation.toMatrix();

                       newPath=Path(nameText,closed)                       
                       for point in obj.Points:
                           x2 = x1 + rotationMatrix.A11*point.x + rotationMatrix.A12*point.y
                           y2 = y1 + rotationMatrix.A21*point.x + rotationMatrix.A22*point.y
                           z2 = z1 + rotationMatrix.A31*point.x + rotationMatrix.A32*point.y
                           newPath.add_point(x2,y2,z2)
                       pathList.append(newPath)

                    if (name[0:7] == "Polygon"):
                       newPath=Path(nameText,"true")
                       for wire in obj.Shape.Wires:
                          for vertex in wire.Vertexes:
                             newPath.add_point(vertex.X,vertex.Y,vertex.Z)
                       pathList.append(newPath)

                if (commandSymbol == "_B"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    portName=get_commandName(label,True)
                    if (portName == "null"):
                        fail = True
                        continue

                    type = ""
                    fullTypeName = ""
                    material = ""
                    numOptions = get_optionCount(label)

                    if (numOptions == 0):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted with missing parameters.\n")
                        fail = True
                        continue

                    type=get_commandOption(label,1)

                    if (type == "null"):
                        fail = True
                        continue

                    if (type == "PEC"):
                        fullTypeName="perfect_electric_conductor"
                        if (numOptions > 1):
                           App.Console.PrintError("ERROR: \"" + label + "\" must not include a second parameter.\n")
                           fail = True
                           continue

                    if (type == "PMC"):
                        fullTypeName="perfect_magnetic_conductor"
                        if (numOptions > 1):
                           App.Console.PrintError("ERROR: \"" + label + "\" must not include a second parameter.\n")
                           fail = True
                           continue

                    if (type == "SI"):
                        fullTypeName="surface_impedance"
                        if (numOptions > 1):
                            material=get_commandOption(label,2)
                            if (material == "null"):
                                fail = True
                                continue

                            if (numOptions > 2):
                                App.Console.PrintError("ERROR: \"" + label + "\" includes too many parameters.\n")
                                fail = True
                                continue
                        else:
                            App.Console.PrintError("ERROR: \"" + label + "\" must include the material parameter.\n")
                            fail = True
                            continue

                    paths=get_paths(label)
                    if (paths == "null"):
                        fail = True
                        continue

                    newBoundary=Boundary(portName,fullTypeName,material)

                    i=0
                    while (i < len(paths)):

                        direction="+"
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            direction="-"
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        newBoundary.add_path(direction,pathName)
                        
                        i=i+1

                    boundaryList.append(newBoundary)

                if (commandSymbol == "_M" or commandSymbol == "_L"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    portName=get_commandName(label,True)
                    if (portName == "null"):
                        fail = True
                        continue

                    name = ""
                    type = ""
                    scale = ""
                    numOptions = get_optionCount(label)

                    if (numOptions == 0):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted with missing parameters.\n")
                        fail = True
                        continue

                    try:
                        name=int(portName)
                    except Exception:
                        App.Console.PrintError("ERROR: \"" + label + "\" does not use an integer mode/line number.\n")
                        fail = True
                        continue

                    if (numOptions > 0):
                        type=get_commandOption(label,1)
                        if (type == "null"):
                            fail = True
                            continue

                    if (numOptions > 1):
                        scale=get_commandOption(label,2)
                        if (scale == "null"):
                            fail = True
                            continue

                    if (numOptions > 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted with excess parameters.\n")
                        fail = True
                        continue

                    paths=get_paths(label)
                    if (paths == "null"):
                        fail = True
                        continue

                    newMode=Mode(name,type,scale)
                    if (commandSymbol == "_M"):
                        newMode.set_category("modal")
                    if (commandSymbol == "_L"):
                        newMode.set_category("line")

                    i=0
                    while (i < len(paths)):

                        direction="+"
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            direction="-"
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        newMode.add_path(direction,pathName)

                        i=i+1

                    modeList.append(newMode)

                if not recognizedCommand and label[0:1] == '_':
                    App.Console.PrintError("ERROR: \"" + label + "\" calls for the invalid command \"" + commandSymbol + "\".\n")
                    fail = True

            # error checks

            # modes
            for mode in modeList:
                if mode.check(pathList):
                    fail=True

            if sequential_mode_check(modeList,"voltage"):
                fail = True
            if sequential_mode_check(modeList,"current"):
                fail = True

            # boundaries
            for boundary in boundaryList:
                if boundary.check(pathList):
                    fail=True

            # write the file

            if not fail:

                # file type and version number
                file.write("#OpenParEMmodes 1.0\n\n")

                # File block
                file.write("File\n")
                file.write("   name=" + str(doc.FileName) + "\n")
                file.write("EndFile\n\n")

                # paths
                for path in pathList:
                    path.print(file)

                # boundaries
                for boundary in boundaryList:
                    boundary.print(file)

                # modes
                for mode in modeList:
                    mode.print(file)


        except Exception:
            App.Console.PrintError("Undetermined fatal error while saving.\n")
        finally:
            file.close()
            if not fail:
                App.Console.PrintMessage("Saved file " + SaveName + ".\n")
    except Exception:
        App.Console.PrintError("Error Open file "+SaveName+".\n")


