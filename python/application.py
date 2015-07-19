#!/usr/bin/env python

# Copyright (C) 2015  Stefano Martina

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


import exercises as exe

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

import Tkinter as tk
import sys
from functools import partial

class E(tk.Tk):
    def __init__(self,parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent


        self.protocol("WM_DELETE_WINDOW", self.dest)
        self.attributes('-zoomed', True)
        self.main()

    def main(self):
        self.menu = {
            'Bezier' : {
                'exercise 1.2' : exe.exer1_2,
                'exercise 1.3' : exe.exer1_3,
                'exercise 1.4' : exe.exer1_4,
                'exercise 1.5' : exe.exer1_5,
                'exercise 1.6' : exe.exer1_6,
                'exercise 1.7' : exe.exer1_7,
                'exercise 2.1' : exe.exer2_1,
                'exercise 2.2' : exe.exer2_2,
                'exercise 2.3' : exe.exer2_3,
                'exercise 2.4' : exe.exer2_4,
            },
            'B-Splines' : {
                #'t1' : exe.test1,
                'example 1' : exe.example1,
                'example 2' : exe.example2,
                'example 3' : exe.example3,
                'example 4' : exe.example4,
                'example 5' : exe.example5,
                'exercise 1' : exe.exerB_1,
                'exercise 2' : exe.exerB_2,
                'exercise 3' : exe.exerB_3,
                'exercise 4' : exe.exerB_4,
                'exercise 5' : exe.exerB_5,
                'exercise 6a' : exe.exerB_6a,
                'exercise 6b' : exe.exerB_6b,
            }
        }

        self.fig = plt.figure(figsize=(14,7))

        self.upFrame = tk.Frame(self)
        self.upFrame.pack(fill=tk.BOTH)

        self.mainFrame = tk.Frame(self.upFrame)
        self.mainFrame.pack(side=tk.RIGHT, fill=tk.BOTH)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.mainFrame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.plot = self.fig.add_subplot(111)

        self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.mainFrame )
        self.toolbar.update()
        self.toolbar.pack()

        self.controlFrame = tk.Frame(self.upFrame)
        self.controlFrame.pack(side=tk.LEFT)

        tk.Label(self.controlFrame, text='Choose the curve type:', padx = 20).pack()
        self.exerType = tk.StringVar()
        for currType in sorted(self.menu.keys()):
            tk.Radiobutton(self.controlFrame, text=currType, variable=self.exerType, value=currType, command=self.radioClick).pack(anchor=tk.W)

        tk.Label(self.controlFrame, text='Choose the exercise:', padx = 20).pack()

        self.exerButtons = []

        self.textFrame = tk.Frame(self)
        self.textFrame.pack(fill=tk.X)
        self.scrollbar = tk.Scrollbar(self.textFrame)
        self.scrollbar.pack(side=tk.LEFT, fill=tk.Y)
        
        self.textArea = tk.Text(self.textFrame, height=4)

        self.textArea.pack(side=tk.RIGHT, fill=tk.X, expand=True)
        self.textArea.config(yscrollcommand=self.scrollbar.set, state=tk.DISABLED)
        self.scrollbar.config(command=self.textArea.yview)

    def radioClick(self):
        for butt in self.exerButtons:
            butt.destroy()
        for exer in sorted(self.menu[self.exerType.get()].keys()):
            self.exerButtons.append(tk.Button(self.controlFrame,text=exer,command=partial(self.execExer, exer)))
            self.exerButtons[-1].pack(fill=tk.BOTH, expand=True)


    def execExer(self, exer):
        self.plot.clear()
        plt.clf()
        self.fig.suptitle(self.exerType.get()+' - '+exer)

        self.menu[self.exerType.get()][exer]()

        self.canvas.draw()

        helpText = self.menu[self.exerType.get()][exer].__doc__
        self.textArea.config(state=tk.NORMAL)
        self.textArea.delete(1.0, tk.END)
        if helpText != None:
            self.textArea.insert(tk.END, helpText)
        self.textArea.config(state=tk.DISABLED)
        
    def dest(self):
        self.destroy()
        sys.exit()



if __name__ == "__main__":
    app = E(None)
    app.title('Exercises')
    app.mainloop()
    
