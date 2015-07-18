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
        self.main()

    def main(self):
        self.menu = {
            'Bezier' : {
                '1.2' : exe.exer1_2,
                '1.3' : exe.exer1_3,
                '1.4' : exe.exer1_4,
                '1.5' : exe.exer1_5,
                '1.6' : exe.exer1_6,
                '1.7' : exe.exer1_7,
                '2.1' : exe.exer2_1,
                '2.2' : exe.exer2_2,
                '2.3' : exe.exer2_3,
                '2.4' : exe.exer2_4,
            },
            'B-Splines' : {
                't1' : exe.test1,
                'e1' : exe.example1,
                'e2' : exe.example2,
                'e3' : exe.example3,
                'e4' : exe.example4,
                'b1' : exe.exerB_1,
                'b2' : exe.exerB_2,
                'b3' : exe.exerB_3,
                'b4' : exe.exerB_4,
                'b5' : exe.exerB_5,
                'b6a' : exe.exerB_6a,
                'b6b' : exe.exerB_6b,
            }
        }

        self.fig = plt.figure(figsize=(13,7))

        self.mainFrame = tk.Frame(self)
        self.mainFrame.grid(row=0, column=1)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.mainFrame)
        self.canvas.get_tk_widget().pack(side='top', fill='both')
        self.canvas._tkcanvas.pack(side='top', fill='both', expand=1)

        self.plot = self.fig.add_subplot(111)

        self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.mainFrame )
        self.toolbar.update()
        self.toolbar.pack()

        self.controlFrame = tk.Frame(self)
        self.controlFrame.grid(row=0, column=0)

        tk.Label(self.controlFrame, text='Choose the curve type:', padx = 20).pack()
        self.exerType = tk.StringVar()
        for currType in self.menu.keys():
            tk.Radiobutton(self.controlFrame, text=currType, variable=self.exerType, value=currType, command=self.radioClick).pack(anchor=tk.W)

        tk.Label(self.controlFrame, text='Choose the exercise:', padx = 20).pack()

        self.exerButtons = []


    def radioClick(self):
        for butt in self.exerButtons:
            butt.destroy()
        for exer in self.menu[self.exerType.get()].keys():
            self.exerButtons.append(tk.Button(self.controlFrame,text=exer,command=partial(self.execExer, self.menu[self.exerType.get()][exer])))
            self.exerButtons[-1].pack()


    def execExer(self, exer):
        self.plot.clear()
        plt.clf()

        exer()

        self.canvas.draw()
        
    def dest(self):
        self.destroy()
        sys.exit()



if __name__ == "__main__":
    app = E(None)
    app.title('Exercises')
    app.mainloop()
    
