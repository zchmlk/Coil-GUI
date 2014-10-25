# -*- coding: utf-8 -*-
"""
Created on Sat Oct 04 14:05:40 2014

@author: zchmlk
"""
from coil_object import Coil
from coil_set_object import CoilSet

from Tkinter import *
import pylab as P
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure



import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

###############################################################################

fdtype = P.float64
idtype = P.int16

widthFactor = 275.0
heightFactor = 400.0
rectScaleFactor = 10.0 # 1 / rsf
idCoords = {}
canvasMade = False
contourf = False
    
def makeCanvas():
    
    print "Creating mesh..."
    global canvasMade
    canvasMade = True    
    createCanvas.config(state='disabled')
    global root
    root = Tk.Tk()
    root.wm_title("Coil Code GUI")
    root.resizable(0,0)
    
    global width
    width = (r_range.get() * widthFactor)
    global height
    height = (z_range.get() * heightFactor)

    global cs   
    cs = CoilSet(coils=[])    
    
    global window
    window = Figure(figsize=(width / 100.0, height / 100.0), dpi=100)
    
    global canvas
    canvas = FigureCanvasTkAgg(window, master=root)
    window.canvas.get_tk_widget().bind("<Button-3>", addCoilEvent)
    window.canvas.get_tk_widget().bind("<Button-1>", removeCoilEvent)   
    
    makeSubPlot()
    
    canvas.show()
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    root.protocol("WM_DELETE_WINDOW", on_quit)
    
    print "Right click to add coil with selected parameters"
    print "Left click to remove coil(s)"
    print "Flux contours will be plotted dynamically"
    print "Exit mesh window in order to create new mesh"
    
def makeSubPlot(): 
    global subPlot
    subPlot = window.add_subplot(111)
    
    subPlot.set_title('Coil Code GUI')
    subPlot.set_xlabel('R Coil')
    subPlot.set_xlim([0, r_range.get()])
    subPlot.set_ylabel('Z Coil')
    subPlot.set_ylim([-z_range.get(), z_range.get()])

def addCoilEvent(event):
    inv = subPlot.transData.inverted()
    coordinates = inv.transform((event.x,event.y))
    rPlot = coordinates[0]
    zPlot = -coordinates[1]
    if rPlot < r_range.get() and rPlot > 0.0 and zPlot > -z_range.get() and zPlot < z_range.get():
        addCoil(event.x, event.y)
    else: 
        pass

def addCoil(r = 0.0, z = 0.0):
    inv = subPlot.transData.inverted()
    coordinates = inv.transform((r, z))
    rPlot = coordinates[0]
    zPlot = -coordinates[1]
    coil = Coil(r_coil = rPlot, z_coil = zPlot, i_coil = 
    i_coil.get(), r_width = r_width.get(), z_width = z_width.get(), 
    n_r_filaments = n_r_filaments.get(), n_z_filaments = n_z_filaments.get())
    
    cs.add_Coil(coil)     
    
    plotFieldContours()
    
    coordinates = subPlot.transData.transform((rPlot, zPlot))
    
    coilNum = canvas.get_tk_widget().create_rectangle(coordinates[0] - (width / rectScaleFactor) * r_width.get() / 2, z - 
                                (height / rectScaleFactor) * z_width.get() / 2 , coordinates[0] + 
                                (width / rectScaleFactor) * r_width.get() / 2, z + (height / rectScaleFactor)
                                * z_width.get() / 2, fill="red")
                            
    global coilCoords
    coilCoords = str(rPlot) + "," + str(zPlot) + "," + str(r_width.get()) + "," + str(z_width.get())
    global idCoords
    idCoords[coilNum] = coilCoords 
    print "Coil added"

def removeCoilEvent(event):
    removeCoil(event.x, event.y)

def removeCoil(r = 0.0, z = 0.0):
    coilsToRemove = canvas.get_tk_widget().find_overlapping(r - r_width.get() 
                                           / 200.0, z -
                                           z_width.get() / 200.0 , r +
                                           r_width.get() / 200.0, z +
                                           z_width.get() / 200.0)                           
    for coil in coilsToRemove:
        if idCoords.has_key(coil):
            canvas.get_tk_widget().delete(coil)
            coords = idCoords[coil].split(",")
            idCoords[coil] = None
            r = coords[0]
            z = coords[1]
            for setcoil in cs.coils:
                coilR = str(setcoil.r_coil)
                coilZ = str(setcoil.z_coil)                
                if r == coilR and z == coilZ:
                    print "Coil removed"
                    cs.remove_Coil(setcoil)
                    setcoil = None;
                else: 
                    pass
    
    plotFieldContours()
                    
def plotFieldContours():
    window.clear()
    makeSubPlot()
    cs.calc_fields(r_range=P.array([1e-6, r_range.get()],
                                               dtype=fdtype),
                                      z_range=P.array([-z_range.get(), z_range.get()],
                                               dtype=fdtype), r_dim=r_dims.get(), z_dim=z_dims.get()) 
    if contourf_var.get() and cs.n_coils != 0:
        cp = subPlot.contourf(cs.r_arr, cs.z_arr, cs.psi_cumulative) 
        makeColorBar(cp)        

        
    elif cs.n_coils != 0:
        cp = subPlot.contour(cs.r_arr, cs.z_arr, cs.psi_cumulative) 
        makeColorBar(cp)
       
    canvas.show()
    
def makeColorBar(cp = None):
    cb = window.colorbar(cp)
    cb.set_label("Psi Levels")
    

def on_quit():
    createCanvas.config(state='normal')
    global canvasMade    
    canvasMade = False
    root.destroy()

def on_quit_parameters():
    if "cs" in globals():
        print "Access coil set attributes using cs. notation"
        print "Number of coils: cs.n_coils = " + str(cs.n_coils)
    if canvasMade:    
        root.destroy()    
    master.destroy()
    
###############################################################################
'''Code for mesh/coil parameters window'''

print "Input mesh and coil parameters"

master = Tk.Tk()
master.wm_title("Parameters")

hit_pop_var = IntVar()
contourf_var = IntVar()
scrollBarSize = 250.0

def hit_pop():
    if hit_pop_var.get():
        r_range.set(4.0)
        z_range.set(1.5)
        r_dims.set(100)
        z_dims.set(100)
        i_coil.set(1.0)
        r_width.set(0.1)
        z_width.set(0.1)
        n_r_filaments.set(10)
        n_z_filaments.set(10)  

hit_pop_button = Checkbutton(master, text="Hit-Pop Settings?", command=hit_pop, variable=hit_pop_var)

hit_pop_button.pack()

contourf_button = Checkbutton(master, text="Filled Contours?", variable=contourf_var)

contourf_button.pack()

r_range = Scale(master, length= scrollBarSize, label="r_range (0 to r_range)", from_=0, to=10, 
               resolution=0.05, orient=HORIZONTAL, fg="black", relief=RAISED, 
               highlightbackground="black", 
               font="helvetica", activebackground="yellow")
r_range.pack()

z_range = Scale(master, length= scrollBarSize, label="z_range (-z_range to z_range)", from_=0, to=5, 
                resolution=0.05, orient=HORIZONTAL, fg="black", relief=RAISED,
                highlightbackground="black", activebackground="yellow",  
                font="helvetica")
z_range.pack()

r_dims = Scale(master, length= scrollBarSize, label="r_dim", from_=0, to=300, 
               resolution=1, orient=HORIZONTAL, fg="black", relief=RAISED, 
               highlightbackground="black", 
               font="helvetica", activebackground="yellow")
r_dims.pack()

z_dims = Scale(master, length= scrollBarSize, label="z_dim", from_=0, to=300, 
                resolution=1, orient=HORIZONTAL, fg="black", relief=RAISED,
                highlightbackground="black", activebackground="yellow",  
                font="helvetica")
z_dims.pack()

createCanvas = Button(master, text="Make Mesh", command=makeCanvas, fg="black",
                      highlightbackground="black", relief=RAISED,
                      font="helvetica", activebackground="yellow")
                      
createCanvas.pack()

i_coil = Scale(master, length= scrollBarSize, label="current", from_=0, to=10, 
               resolution=0.05, orient=HORIZONTAL, fg="black", relief=RAISED, 
               highlightbackground="black", 
               font="helvetica", activebackground="yellow")
i_coil.pack()

r_width = Scale(master, length= scrollBarSize, label="r_width", from_=0, to=1, 
                resolution=0.05, orient=HORIZONTAL, fg="black", relief=RAISED,
                highlightbackground="black", activebackground="yellow",  
                font="helvetica")
r_width.pack()

z_width = Scale(master, length= scrollBarSize, label="z_width", from_=0, to=1, 
                resolution=0.05, orient=HORIZONTAL, fg="black", relief=RAISED, 
                highlightbackground="black",
                font="helvetica", activebackground="yellow")
z_width.pack()

n_r_filaments = Scale(master, length= scrollBarSize, label="n_r_filaments",
                      from_=0, to=100, orient=HORIZONTAL, fg="black", 
                      relief=RAISED, highlightbackground="black", 
                      font="helvetica", activebackground="yellow")
n_r_filaments.pack()

n_z_filaments = Scale(master,  length= scrollBarSize, label="n_z_filaments",
                      from_=0, to=100, orient=HORIZONTAL, fg="black",
                      highlightbackground="black", relief=RAISED,
                      font="helvetica", activebackground="yellow")
n_z_filaments.pack()

master.protocol("WM_DELETE_WINDOW", on_quit_parameters)

Tk.mainloop()


