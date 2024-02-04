# -*- coding: utf-8 -*-
import FreeCAD,Draft
import PySide
from PySide import QtGui ,QtCore
from PySide.QtGui import *
from PySide.QtCore import *

doc = FreeCAD.ActiveDocument
objs = FreeCAD.ActiveDocument.Objects

file = doc.Name+"_modes.txt"
path = "./" + file

try:
    SaveName = QFileDialog.getSaveFileName(None,QString.fromLocal8Bit("Save the lines file"),path,             "*.txt") # PyQt4
#                                                                     "here the text displayed on windows" "here the filter (extension)"   
except Exception:
    SaveName, Filter = PySide.QtGui.QFileDialog.getSaveFileName(None, "Save the lines file", path,             "*.txt") # PySide
#                                                                     "here the text displayed on windows" "here the filter (extension)"   
if SaveName == "":                                                            # if the name file are not selected, then abort the process
    App.Console.PrintMessage("Process aborted\n")
else:                                                                         # if the name file are selected or created then 
    try:                                                                      # detect error ...
        file = open(SaveName, 'w')                                            # open the file selected to write (w)
        try:                                                                  # if error detected to write ...
            App.Console.PrintMessage("Process running ...\n")

            file.write("#OpenParEMmodes 1.0\n\n")

            file.write("File\n")
            file.write("   name=" + str(doc.FileName) + "\n")
            file.write("EndFile\n\n")

            for obj in objs:

                #if (not obj.ViewObject.Visibility):
                #    continue

                name = obj.Name                                             # list the Name  of the object  (not modifiable)
                label = obj.Label                                           # list the Label of the object  (modifiable)

                if (label[0:2] == "_P"):
                    print ("processing:" + label)

                    nameText=label[2:]
                    if (nameText == ''):
                       App.Console.PrintError("ERROR: \"" + label +"\" does not include a name.\n")
                       continue

                    if (name[0:4] == "Line"):

                       spos_x = obj.Start.x
                       spos_y = obj.Start.y
                       epos_x = obj.End.x
                       epos_y = obj.End.y

                       file.write("Path\n")
                       file.write("   name=" + str(nameText) + "\n")
                       file.write("   point=(" + str(spos_x) + "," + str(spos_y) + ")\n")
                       file.write("   point=(" + str(epos_x) + "," + str(epos_y) + ")\n")
                       file.write("   closed=false\n")
                       file.write("EndPath\n\n")

                    if (name[0:9] == "Rectangle"):
  
                       x1 = obj.Placement.Base.x
                       y1 = obj.Placement.Base.y
                       x2 = obj.Placement.Base.x + obj.Length.Value
                       y2 = obj.Placement.Base.y
                       x3 = obj.Placement.Base.x + obj.Length.Value
                       y3 = obj.Placement.Base.y + obj.Height.Value
                       x4 = obj.Placement.Base.x 
                       y4 = obj.Placement.Base.y + obj.Height.Value

                       file.write("Path\n")
                       file.write("   name=" + str(nameText) + "\n")
                       file.write("   point=(" + str(x1) + "," + str(y1) + ")\n")
                       file.write("   point=(" + str(x2) + "," + str(y2) + ")\n")
                       file.write("   point=(" + str(x3) + "," + str(y3) + ")\n")
                       file.write("   point=(" + str(x4) + "," + str(y4) + ")\n")
                       file.write("   closed=true\n")
                       file.write("EndPath\n\n")

                    if (name[0:4] == "Wire"):
                       print ("processing:" + label)
 
                       closed="false"
                       if (obj.Closed):
                           closed="true"

                       file.write("Path\n")
                       file.write("   name=" + str(nameText) + "\n")
                       for point in obj.Points:
                           x = obj.Placement.Base.x + point.x
                           y = obj.Placement.Base.y + point.y
                           file.write("   point=(" + str(x) + "," + str(y) + ")\n")
                       file.write("   closed=" + str(closed) + "\n")
                       file.write("EndPath\n\n")

                if (label[0:2] == "_B" and name[0:4] == "Text"):
                    print ("processing:" + label)

                    # get the name
                    splitText=label[2:].split('(')
                    if (len(splitText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    if (splitText[0] == ''):
                        App.Console.PrintError("ERROR: \"" + label + "\" is missing a name.\n")
                        continue
                    nameText=splitText[0]

                    # get the type
                    typeText=splitText[1].split(')')
                    if (len(typeText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    type=typeText[0]
                    if (type != 'SI' and type != 'PEC' and type != 'PMC'):
                        App.Console.PrintError("ERROR:\"" + label + "\" has an unsupported type of \"" + type + "\".\n")
                        continue

                    # get the path
                    pathText=typeText[1].split('{')
                    if (len(pathText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    pathText=pathText[1].split('}')
                    if (len(pathText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                     
                    paths=pathText[0].split(',')
                    if (len(paths) == 0):
                        App.Console.PrintError("ERROR: \"" + label + "\" is missing paths.\n")
                        continue

                    file.write("Boundary\n")
                    file.write("   name=" + str(nameText) + "\n")
                    if (type == "SI"):
                        file.write("   type=surface_impedance\n")
                    if (type == "PEC"):
                        file.write("   type=perfect_electric_conductor\n")
                    if (type == "PMC"):
                        file.write("   type=perfect_magnetic_conductor\n")

                    i=0
                    while (i < len(paths)):

                        is_positive=True
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            is_positive=False
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        if (i == 0):
                            if (is_positive):
                                file.write("   path=" + pathName + "\n")
                            else:
                                file.write("   path=-" + pathName + "\n")
                        else:
                            if (is_positive):
                                file.write("   path+=" + pathName + "\n")
                            else:
                                file.write("   path-=" + pathName + "\n")

                        i=i+1

                    file.write("EndBoundary\n\n")

                if ((label[0:2] == "_M" and name[0:4] == "Text") or (label[0:2] == "_L" and name[0:4])):
                    print ("processing:" + label)

                    # get the name
                    splitText=label[2:].split('(')
                    if (len(splitText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    if (splitText[0] == ''):
                        App.Console.PrintError("ERROR: \"" + label + "\" is missing a name.\n")
                        continue
                    nameText=splitText[0]
                    try:
                        mode=int(nameText)
                    except Exception:
                        App.Console.PrintError("ERROR: \"" + label + "\" does not use an integer mode indication.\n")
                        continue

                    # get the type
                    typeText=splitText[1].split(')')
                    if (len(typeText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    type=typeText[0]
                    if (type != 'voltage' and type != 'current'):
                        App.Console.PrintError("ERROR:\"" + label + "\" has an unsupported type of \"" + type + "\".\n")
                        continue

                    # get the path
                    pathText=typeText[1].split('{')
                    if (len(pathText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue
                    pathText=pathText[1].split('}')
                    if (len(pathText) != 2):
                        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
                        continue

                    paths=pathText[0].split(',')
                    if (len(paths) == 0):
                        App.Console.PrintError("ERROR: \"" + label + "\" is missing paths.\n")
                        continue

                    if (label[0:2] == "_M"):
                        file.write("Mode\n")
                        file.write("   mode=" + str(mode) + "\n")
                    if (label[0:2] == "_L"):
                        file.write("Line\n")
                        file.write("   line=" + str(mode) + "\n")
                    if (type == "voltage"):
                        file.write("   type=voltage\n")
                    if (type == "current"):
                        file.write("   type=current\n")

                    i=0
                    while (i < len(paths)):

                        is_positive=True
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            is_positive=False
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        if (i == 0):
                            if (is_positive):
                                file.write("   path=" + pathName + "\n")
                            else:
                                file.write("   path=-" + pathName + "\n")
                        else:
                            if (is_positive):
                                file.write("   path+=" + pathName + "\n")
                            else:
                                file.write("   path-=" + pathName + "\n")

                        i=i+1

                    if (label[0:2] == "_M"):
                        file.write("EndMode\n\n")
                    if (label[0:2] == "_L"):
                        file.write("EndLine\n\n")


        except Exception:                                                     # if error detected to write
            App.Console.PrintError("Error write file "+"\n")                  # detect error ... display the text in red (PrintError)
        finally:                                                              # if error detected to write ... or not the file is closed
            file.close()                                                      # if error detected to write ... or not the file is closed
            App.Console.PrintMessage("Saved file " + SaveName + ".\n")
    except Exception:
        App.Console.PrintError("Error Open file "+SaveName+"\n")      # detect error ... display the text in red (PrintError)


