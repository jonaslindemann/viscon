# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 11:02:04 2017

@author: Karin
"""

#-----------------------------IMPORTERING--------------------------------------

import consolidation_model as fm
import sys
import json
import numpy as np
import calfem.vis as cfv
import cfext as tg

from PyQt5.QtWidgets import QApplication, QTableWidgetItem, QMainWindow, QSlider, QMessageBox, QFileDialog
from PyQt5.uic import loadUi
from PyQt5.QtCore import pyqtSlot, QThread, pyqtSignal
from PyQt5 import QtGui, QtCore

    
    
        
# -----------------------MAINWINDOW--------------------------------------------
class MainWindow(QMainWindow):
    """MainWindow-klass som hanterar vårt huvudfönster"""

    def __init__(self):
        """Konstruktor"""
        
        super(MainWindow, self).__init__()

        # Load user interface from UI-file
        loadUi('viscon.ui',self)
         
        # Create figure class
        Figure = cfv.figureClass()  
        
        #Skapa figurer för förhandsvisning av geo, mesh, initialspänningar
        self.figGeometry = Figure(self)
        self.figMesh = Figure(self)
        self.figInitial = Figure(self)
        
        self.middleLayout.addWidget(self.figGeometry._widget, 20)
        self.middleLayout.addWidget(self.figMesh._widget, 20)
        self.middleLayout.addWidget(self.figInitial._widget,20)
        
        #Figurer för resultat
        self.figPress = Figure(self)
        self.figuvec = Figure(self)
        self.figsigvec = Figure(self) 

        
        self.resLayout1.addWidget(self.figPress._widget, 20)
        self.resLayout2.addWidget(self.figuvec._widget, 20)  
        self.resLayout2.addWidget(self.figsigvec._widget, 20)  
          
        #Koppla menyval till händelser
        self.actionNew.triggered.connect(self.onActionNew)
        self.actionOpen.triggered.connect(self.onActionOpen)
        self.actionSave.triggered.connect(self.onActionSave)
        self.actionSaveas.triggered.connect(self.onActionSaveas)
        self.actionExit.triggered.connect(self.onActionExit)

        
        #Koppla reglage, knappar etc till händelser
        self.comboBox.currentIndexChanged.connect(self.onDirectionChanged)
        self.timeSlider.valueChanged.connect(self.onTimeChanged)
        self.localSlider.valueChanged.connect(self.onLocalChanged)       
        self.meshSlider.valueChanged.connect(self.onMeshChanged)  
        
        self.meshSlider.setMinimum(2)
        self.meshSlider.setMaximum(10)
        self.meshSlider.setTickPosition(QSlider.TicksBelow)
        self.meshSlider.setTickInterval(1)
     
        # --- Visa fönster
        self.show()
        self.raise_()
        
        # --- Hämta värden från initModel och uppdatera kontrollerna med dessa        
        self.initModel()
        self.updateControls()
        

    def initModel(self):
        """ Skapa initiell modell"""
        
        self.inputData = fm.InputData()
        self.outputData = fm.OutputData()
        self.filename = None   
        
        
    # ----   UPPDATERA KONTROLL OCH MODELLER

    def updateControls(self):
        """Uppdatera kontrollerna med värden från modellen"""
        
        #Anpassa CheckBox efter grundens permeabilitet
        self.checkBox.setChecked(self.inputData.perm)
        
        #Uppdatera kontroller med givna värden på inputData
        self.bEdit.setText(str(self.inputData.b))
        self.qEdit.setText(str(self.inputData.q))
        self.BEdit.setText(str(self.inputData.B))
        self.gammaEdit.setText(str(self.inputData.rhow))

        self.TEdit.setText(str(self.inputData.T))        
        self.dtEdit.setText(str(self.inputData.dt))   
        self.TnrEdit.setText(str(self.inputData.Tnr))   
        
        self.XnrEdit.setText(str(self.inputData.Xnr))
        self.YnrEdit.setText(str(self.inputData.Ynr))
        
        self.meshSlider.setValue(self.inputData.elSizenr)
        self.meshLabel.setText(str(self.inputData.elSize))
        
        #Uppdatera tabellvärden
        self.tableLayer.setRowCount(self.inputData.hn)
        for i in range(self.inputData.hn):  
            self.tableLayer.setItem(i, 0, QTableWidgetItem(str(float(self.inputData.hvec[i]))))
            self.tableLayer.setItem(i, 1, QTableWidgetItem(str(float(self.inputData.kmat[i,0]))))
            self.tableLayer.setItem(i, 2, QTableWidgetItem(str(float(self.inputData.kmat[i,1]))))
            self.tableLayer.setItem(i, 3, QTableWidgetItem(str(float(self.inputData.Mvec[i]))))
            

    def updateModel(self):
       """Hämta värden från kontroller och uppdatera modellen"""
       
       #Justera grundens permeabilitet efter om grunden är permeabel
       if self.checkBox.isChecked():
            self.inputData.perm = True
       else:
           self.inputData.perm = False

       
       #Uppdatera modellen med givna värden
       self.inputData.b = float(self.bEdit.text())
       self.inputData.q = float(self.qEdit.text())
       self.inputData.B = float(self.BEdit.text())
       self.inputData.rhow = float(self.gammaEdit.text())
       
       
       self.inputData.elSizenr = float(self.meshSlider.value())
       self.inputData.elSize = self.inputData.elSizenr /10
       
       self.inputData.dt = float(self.dtEdit.text())
       self.inputData.T = float(self.TEdit.text())
       self.inputData.Tnr = float(self.TnrEdit.text())

           
       
       self.inputData.hvec = np.zeros([self.inputData.hn,1])
       self.inputData.kmat = np.zeros([self.inputData.hn,2])
       self.inputData.Mvec = np.zeros([self.inputData.hn,1])           
       for i in range (self.inputData.hn):
            self.inputData.hvec[i] = self.tableLayer.item(i, 0).text()
            self.inputData.kmat[i,0] = self.tableLayer.item(i, 1).text()
            self.inputData.kmat[i,1] = self.tableLayer.item(i, 2).text()
            self.inputData.Mvec[i] = self.tableLayer.item(i, 3).text()
        
        
       self.inputData.Xnr = float(self.XnrEdit.text())
       self.inputData.Ynr = float(self.YnrEdit.text())
    
        
       #Skapa varningar baserat på orimliga värden i inputdata
       warnings = {}
       if self.inputData.b <= 0:               
           warnings["b"]="b has to be > 0" 
                    
       if self.inputData.B <= 0:
           warnings["B"]="B has to be > 0"

       if self.inputData.B <= self.inputData.b:
           warnings["Bb"]="B has to be >b"

       if any(x<=0 for x in self.inputData.Mvec):
           warnings["Mvec"]="All values for M have to be > 0"

       if self.inputData.hn <= 0:               
           warnings["hn"]="hn has to be > 0" 
                    
       if any(x<=0 for x in self.inputData.hvec):
           warnings["hvec"]="All values in hvec have to be > 0"

       if self.inputData.T <= 0:
           warnings["T"]="T has to be > 0"

       if self.inputData.dt <= 0:
           warnings["dt"]="dt has to be > 0"
       elif self.inputData.dt >= self.inputData.T:
           warnings["dt"]="dt has to be < T"
       elif self.inputData.Tnr > self.inputData.T/self.inputData.dt:
           warnings["Tnr"]="Tnr has to be less than T/dt (number of time steps evaluated)" 
           
       if self.inputData.Tnr <= 0:
           warnings["Tnr"]="Tnr has to be > 0"

       if self.inputData.Xnr <= 0:
           warnings["Xnr"]="Xnr has to be > 0"
       
       if self.inputData.Ynr <= 0:
           warnings["Ynr"]="Ynr has to be > 0"

       if self.inputData.rhow <= 0:
           warnings["Rhow"]="The unit weight of pore water has to be > 0"

       if self.inputData.Xnr <= 0:
           warnings["Xnr"]="Xnr has to be > 0"

       if self.inputData.Ynr <= 0:
           warnings["Ynr"]="Ynr has to be > 0"

       if self.inputData.rhow <= 0:
           warnings["Ynr"]="Weigth of pore water has to be > 0"


       if any(x<=0 for x in self.inputData.kmat[:,0]):
           warnings["kmat"]="All kx-values in kmat has to be > 0"    
                   
       if any(x<=0 for x in self.inputData.kmat[:,1]):
           warnings["kmat"]="All ky-values in kmat has to be > 0"    
       self.inputData.warnings = warnings

     
          
       
     #-------------ACTIONS:  
       
    @pyqtSlot()              
    def on_addLayer_clicked(self):
        """Lägg till ett nytt layer"""
        #Lägg till rad i tabellen med givna värden
        rowPosition = self.tableLayer.rowCount()
        self.tableLayer.insertRow(rowPosition)
        self.tableLayer.setItem(rowPosition, 0, QTableWidgetItem(str(3)))
        self.tableLayer.setItem(rowPosition, 1, QTableWidgetItem(str(1e-9)))
        self.tableLayer.setItem(rowPosition, 2, QTableWidgetItem(str(1e-9)))
        self.tableLayer.setItem(rowPosition, 3, QTableWidgetItem(str(900e3)))
    
        #Uppdatera värde för antal lager samt modell.
        self.inputData.hn = self.inputData.hn + 1 
        self.updateModel   
        
    @pyqtSlot()   
    def on_removeLayer_clicked(self):
        """Ta bort ett layer""" 
        rownr = self.tableLayer.rowCount()
        if rownr<=1:
            return
        else:
            rowPosition = rownr - 1
            self.tableLayer.removeRow(rowPosition)

            #Uppdatera värde för lager samt modellen
            self.inputData.hn = self.inputData.hn - 1
            self.updateModel

    @pyqtSlot()
    def on_updateButton_clicked(self):
        """Update preview geo, mesh, stress"""
        
        self.updateModel()  
        warnings = self.inputData.warnings
        
        #Om det saknas varningar/felmeddelanden:
        if not warnings:
            self.clearfig()

            self.g = fm.InputData.geometry(self.inputData)
            fm.Solver.executeMesh(self)

            
            #Hämta värden som används vid plot
            q = self.inputData.q
            B = self.inputData.B
            b = self.inputData.b
            hvec = self.inputData.hvec
            coords = self.outputData.coords
            edof = self.outputData.edof
            dofsPerNode = self.outputData.dofsPerNode
            elType = self.outputData.elType
            ndof = edof.max()
            
            # Bestäm initialspänningarna
            a0 = tg.strdist(coords, ndof, q, B, b, hvec)
            
            # Figur geometri
            cfv.figure(self.figGeometry.nr) 
            cfv.drawGeometry(self.g, title="Geometry")
            
            # Figur mesh
            cfv.figure(self.figMesh.nr) 
            cfv.drawMesh(coords, edof, dofsPerNode, elType, filled=True, title="Mesh")
            
            # Figur initialspänningar 
            cfv.figure(self.figInitial.nr)
            cfv.drawNodalValues(a0, coords, edof, dofsPerNode, elType, clim=(0,self.inputData.q), doDrawMesh = False, title = "Initial stress distr.")
       
            self.tabWidget.setEnabled(True)
        
        #Om det finns felmeddelanden: visa ej preview utan aktivera felmeddelanden
        else:
            self.wmsg()

    @pyqtSlot()
    def on_updateTime_clicked(self):
        """Update time approximation"""
        
        self.updateModel()  
        warnings = self.inputData.warnings
        
        #Om det saknas varningar/felmeddelanden:
        if not warnings:
            #Töm figurer
            self.clearfig()
            
            self.tabWidget.setEnabled(False)
            #Genomför tidsapprox.
            self.solver = fm.Solver(self.inputData, self.outputData)
            self.progressBar = self.progressBar
            self.progresslabel = self.progresslabel
            self.solverThread = SolverThread(self.solver, self.progressBar, self.progresslabel, Mesh = False, Fem = False, Time = True, Extr = True)
            self.solverThread.finished.connect(self.onSolverFinished)   
            self.solverThread.start()
            

        # Om fel i inputdata - utfärda felmeddelande och genomför ej tidsapproximation
        else:
                self.wmsg()          
    
    
    
    @pyqtSlot()
    def on_oneDupdate_clicked(self):
        """Update 1D-plot"""
        
        self.updateModel()  
        warnings = self.inputData.warnings
        
        #Om det saknas varningar/felmeddelanden:
        if not warnings:
            #Töm figurer
            self.clearfig()
            
            self.tabWidget.setEnabled(False)
            #Genomför .
            self.solver = fm.Solver(self.inputData, self.outputData)
            self.progressBar = self.progressBar
            self.progresslabel = self.progresslabel
            self.solverThread = SolverThread(self.solver, self.progressBar, self.progresslabel, Mesh = False, Fem = False, Time = False, Extr = True)
            self.solverThread.finished.connect(self.onSolverFinished)   
            self.solverThread.start()
            

        # Om fel i inputdata - utfärda felmeddelande och genomför ej tidsapproximation
        else:
                self.wmsg()          
    
    
    
    @pyqtSlot(int)
    def on_tabWidget_currentChanged(self, tabIndex):
        """Handle tab change"""

        if tabIndex == 1:   

            #Töm figurer
            
            self.updateModel()
            self.updateControls()
            warnings = self.inputData.warnings
            
            
            #Om felmeddelanden saknas:
            if not warnings:
                
                # --- Starta en tråd för att köra beräkningen, så att gränssnittet inte fryser.
                self.clearfig()
                self.tabWidget.setEnabled(False)
                self.progressBar.setEnabled(True)

                #Genomför samtliga beräkningar
                self.progressBar = self.progressBar
                self.progresslabel = self.progresslabel
                self.solver = fm.Solver(self.inputData, self.outputData) 
                self.solverThread = SolverThread(self.solver, self.progressBar, self.progresslabel,Mesh = True, Fem = True, Time = True, Extr = True)
                self.solverThread.finished.connect(self.onSolverFinished)   
                self.solverThread.start()

    
            #Om felmeddelanden: kör ej och utfärda felmeddelanden   
            else:
                self.tabWidget.setCurrentIndex(0)
                self.wmsg()
    
        
        #Töm plottar när man återgår till tab för inputdata
        if tabIndex == 0:
            self.clearfig()
            self.inputData.Xv = 0
            self.inputData.Yv = 0
            self.inputData.Tv = 0
            
    def onDirectionChanged(self):
        self.inputData.comboindex = self.comboBox.currentIndex()

        if self.inputData.comboindex == 0:
            Xnr = self.inputData.Xnr
            self.localSlider.setMinimum(0)
            self.localSlider.setMaximum(Xnr-1)
            self.localSlider.setValue(0)
            self.localSlider.setTickPosition(QSlider.TicksBelow)
            self.localSlider.setTickInterval(1)
            self.combolabel.setText("Extract vertical 1D-plot at x =")

        
        else:
            Ynr = self.inputData.Ynr
            self.localSlider.setMinimum(0)
            self.localSlider.setMaximum(Ynr-1)
            self.localSlider.setValue(0)
            self.localSlider.setTickPosition(QSlider.TicksBelow)
            self.localSlider.setTickInterval(1)
            self.combolabel.setText("Extract horisontal 1D-plot at z =")
        fm.Visualisation.show(self)   
        
       
    def onMeshChanged(self):

        self.meshLabel.setText(str((self.meshSlider.value())/10))
        
    def onTimeChanged(self, value):
        """Uppdatera värde för slidern """ 
        
        #Plocka upp vilket värde som valts för slidern
        Tvec =np.append(np.array([0]), self.outputData.Tvec)*self.inputData.dt
        Tv = self.timeSlider.value()
        q = self.inputData.q
        Tvalue = Tvec[Tv]
        self.timelabel.setText(str(np.round(Tvalue,2)))
        self.maxu.setText(str(np.round(np.max(self.outputData.avec[:,Tv]),1)))
        self.maxeff.setText(str(np.round(q-np.max(self.outputData.avec[:,Tv]),1)))
        #self.maxtot.setText(str(np.round(np.max(self.outputData.avec[:,0]),1)))
        self.inputData.Tv = Tv
        
        #Se till att visa figurer för resultat för vald slidervärde
        fm.Visualisation.show(self)       
        
    def onLocalChanged(self, value):
        """Uppdatera värde för slidern """ 
        
        #Plocka upp vilket värde som valts för slidern
        if self.inputData.comboindex == 0:
            xlvec = self.outputData.xlvec
        
            Xv = self.localSlider.value()
            Xvalue = xlvec[Xv]
            self.locallabel.setText(str(np.round(Xvalue,2)))
            self.inputData.Xv = Xv

        else:
            ylvec = self.outputData.ylvec
        
            Yv = self.localSlider.value()
            Yvalue = ylvec[Yv]
            self.locallabel.setText(str(np.round(Yvalue,2)))
            self.inputData.Yv = Yv

            
            
        #Se till att visa figurer för resultat för vald slidervärde
        fm.Visualisation.show(self)               
            
        
        
            
    def waitingmsg(self):
       waitmsg = QMessageBox()
       waitmsg.setIcon(QMessageBox.Information)       
       waitmsg.exec_()
    
    def wmsg(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        warnings = self.inputData.warnings
        newline = ""
        for key in warnings:
            newline = newline + "\n" + warnings[key]    
        msg.setText(newline)
        msg.exec_()
        
    def clearfig(self):
        
        cfv.figure(self.figPress.nr)
        cfv.clf() 
        cfv.figure(self.figuvec.nr)
        cfv.clf()
        cfv.figure(self.figsigvec.nr)
        cfv.clf()  

        cfv.figure(self.figGeometry.nr)
        cfv.clf()
        cfv.figure(self.figMesh.nr)
        cfv.clf()
        cfv.figure(self.figInitial.nr)
        cfv.clf()
        
   
# -------------------SOLVERTHREAD----------------------------------------------------------
        
        
        
    def onSolverFinished(self):
         """Anropas när beräkningstråden avslutas"""

         #Uppdatera slider med tidsteg

         Tvec = self.outputData.Tvec
         Tsnap = np.shape(Tvec)[0]
         self.timeSlider.setMinimum(0)
         self.timeSlider.setMaximum(Tsnap)
         self.timeSlider.setValue(0)
         self.timeSlider.setTickPosition(QSlider.TicksBelow)
         self.timeSlider.setTickInterval(1)
         self.maxu.setText("-")
         self.maxeff.setText("-")
         #self.maxtot.setText("-")

        #Uppdatera slider med lägen
         if self.inputData.comboindex == 0:
             Xnr = self.inputData.Xnr
             self.localSlider.setMinimum(0)
             self.localSlider.setMaximum(Xnr-1)
             self.localSlider.setValue(0)
             self.localSlider.setTickPosition(QSlider.TicksBelow)
             self.localSlider.setTickInterval(1)
             self.combolabel.setText("Extract vertical 1D-plot at x =")
 
         
         else:
             Ynr = self.inputData.Ynr
             self.localSlider.setMinimum(0)
             self.localSlider.setMaximum(Ynr-1)
             self.localSlider.setValue(0)
             self.localSlider.setTickPosition(QSlider.TicksBelow)
             self.localSlider.setTickInterval(1)
             self.combolabel.setText("Extract horisontal 1D-plot at z =")
        
        
         self.progressBar.setValue(0)
         self.progresslabel.setText("Ready")
         self.progressBar.setEnabled(False)
         self.tabWidget.setEnabled(True)      
        
        
        
    #--------- KOPPLAT TILL FILHANTERAREN
            
        
    def save(self, filename):
        """Spara indata till fil"""                   
        inputData = {}
        inputData["hn"] = self.inputData.hn
        inputData["hvec"] = self.inputData.hvec.tolist()
        inputData["b"] = self.inputData.b
        inputData["B"] = self.inputData.B
        inputData["q"] = self.inputData.q    
        inputData["kmat"] = self.inputData.kmat.tolist()
        inputData["Mvec"] = self.inputData.Mvec.tolist()
        inputData["elSize"] = self.inputData.elSize
        inputData["elSizenr"] = self.inputData.elSizenr
        inputData["rhow"] = self.inputData.rhow
        inputData["dt"] = self.inputData.dt
        inputData["T"] = self.inputData.T
        inputData["Tnr"] = self.inputData.Tnr
        inputData["Xnr"] = self.inputData.Xnr
        inputData["Ynr"] = self.inputData.Ynr

        
        if self.inputData.perm == True:
            inputData["perm"] = 1
        else:
            inputData["perm"] = 0
      
        ofile = open(filename, "w")
        json.dump(inputData, ofile, sort_keys = True, indent = 4)
        ofile.close()   
        
    def load(self, filename):
        """Hämtning av värden från fil"""
        ifile = open(filename, "r")
        inputData = json.load(ifile)
        ifile.close()
        
        self.inputData.hn = inputData["hn"]
        self.inputData.hvec = np.asarray(inputData["hvec"])
        self.inputData.b = inputData["b"]
        self.inputData.B = inputData["B"]
        self.inputData.q = inputData["q"]
        self.inputData.kmat = np.asarray(inputData["kmat"])
        self.inputData.Mvec = np.asarray(inputData["Mvec"])
        self.inputData.elSize = inputData["elSize"]
        self.inputData.elSizenr = inputData["elSizenr"]
        self.inputData.rhow = inputData["rhow"]
        self.inputData.dt = inputData["dt"]
        self.inputData.T = inputData["T"]
        self.inputData.Tnr = inputData["Tnr"] 
        self.inputData.Xnr = inputData["Xnr"]
        self.inputData.Ynr = inputData["Ynr"]

        

        test = inputData["perm"]
        if test == 1:
            self.inputData.perm = True
        else:
            self.inputData.perm = False
    
    def onActionNew(self):
        """Skapa en ny modell"""
        
        # --- Hämta initmodel och uppdatera kontrollerna med dessa värden
        self.initModel()
        self.updateControls()  
        
        #Töm figurer
        self.clearfig()
        

         
    def onActionOpen(self):
        """Öppna indatafil"""
        
        # --- Fil som ska öppnas väljs i filhanterare 
        self.filename, _ = QFileDialog.getOpenFileName(self,
        "Öppna modell", "", "Modell filer (*.json *.jpg *.bmp)")
        
        # --- Om en fil har valts i filhanteraren laddas värden och kontroller uppdateras
        if self.filename != "":
            self.load(self.filename)
           
            self.tabWidget.setCurrentIndex(0)
            self.updateControls()
            self.updateModel()
            
            #Töm figurer
            self.clearfig()
            
    def onActionSave(self):
        """Spara fil"""
    
        self.updateModel()

        
       # --- Om filen saknar namn öppnas ett fönster för att ange namn och vart filen ska sparas 
        if self.filename == None:
            self.filename, _  = QFileDialog.getSaveFileName(self,
            "Spara modell", "", "Modell filer (*.json)")           
           
        # --- Om filen tilldelats ett namn sparas den i detta
        if self.filename != "":
            self.save(self.filename)
        else:
            self.filename = None
    
    def onActionSaveas(self):
        """Spara fil som...""" 
        
        # --- Skapa ett temporärt filnamn
        self.filenameny = None
        

        self.updateModel()
        
        # --- Öppna filhanterare för att tilldela namn och plats att spara filen på
        self.filenameny, _  = QFileDialog.getSaveFileName(self, "Spara modell", "", "Modell filer (*.json)")
        
        # --- Har filen tilldelats ett nytt namn sparas filen i detta
        if self.filenameny != "":
           
            
            # --- Det riktiga filnamnet ersätts med det tillfälliga
            self.filename = self.filenameny
            self.save(self.filename)

    def onActionExit(self):
        """Stäng programmet"""        
        self.close()
        
        
        
        
        
class SolverThread(QtCore.QThread):
    """Klass för att hantera beräkning i bakgrunden"""

    def __init__(self, solver, progressBar, progresslabel, Mesh = True, Fem = True, Time = True, Extr = True ):
        """Klasskonstruktor"""
        
        QtCore.QThread.__init__(self)
 
        #self.progressBarLabel = progressBarLabel
        self.progressBar = progressBar
        self.progresslabel = progresslabel
        self.solver = solver
        self.Mesh = Mesh
        self.Fem = Fem
        self.Time = Time
        self.Extr = Extr


    def __del__(self):
        """Väntar"""
        self.wait()

    def run(self):
        """Genomför beräkning"""
 
        if self.Mesh == True:
            self.progresslabel.setText("Generating geometry and mesh")
            self.solver.executeMesh()
            self.progressBar.setValue(30)

         


        if self.Fem == True:
            self.progresslabel.setText("Generating FE-model")
            self.solver.executeFem()
            self.progressBar.setValue(60)


        if self.Time == True:
            self.progresslabel.setText("Time approximation")
            self.solver.executeTime()
            self.progressBar.setValue(80)
     
        if self.Extr == True:
            self.progresslabel.setText("Preparing 1D-plot")
            self.solver.executeExtr()
            self.progressBar.setValue(100)    
                

if __name__ == '__main__':    
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec_())  
    

 