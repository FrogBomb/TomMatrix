# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:07:26 2013

@author: Tom
"""

from Tkinter import *
from TomMatrix import *

class TomMatrixEditor():
    def __init__(self, inMatrix):
        self.matrix = inMatrix
        
    def __call__(self, *args, **kwargs):
        
        self.root = Tk(*args, **kwargs)
        self.createEntryMatrix()
        self.buttonBox()
        
        self.root.focus_set()
        
        self.root.mainloop()
    
    def buttonBox(self):
        
        box = Frame(self.root)

        w = Button(box, text="Save", width=10, command=self.save, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Save & Quit", width=10, command=self.quit)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        self.root.bind("<Return>", self.save)
        self.root.bind("<Escape>", self.cancel)
        
        box.pack()
    
    def createEntryMatrix(self):
        box = Frame(self.root)
        numRows, numCols = self.matrix.dims()
        self.EntryMatrix = [Entry(box) for i in range(numRows*numCols)]
        
        for row in range(numRows):
            for col in range(numCols):
                self.EntryMatrix[col + row*numCols].grid(row=row, column=col)
                self.EntryMatrix[col + row*numCols].insert(0, str(float(self.matrix[row, col])))
                
        box.pack()
                
    
    def save(self):
        numRows, numCols = self.matrix.dims()
        for row in range(numRows):
            for col in range(numCols):
                try:
                    self.matrix[row, col] = float(self.EntryMatrix[col + row*numCols].get())
                except ValueError:
                    self.matrix[row, col] = 0
    
    def quit(self):
        self.save()
        self.cancel()
    
    def cancel(self):
        self.root.destroy()