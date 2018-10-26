#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Tkinter as tk
import subprocess
import os




#variable declarations

nwx=10#default x size for a new window =500
nwy=10#default y size for a new window =500



def getColor(value,rangeMax):
	intervalSize=(rangeMax)/3
	Rint=255
	Gint=255
	Bint=255
	RGB='#%02x%02x%02x' % (Rint, Gint, Bint)


	if(value<=intervalSize):
		Rint=(200/intervalSize)*value
		Gint=Rint
		Bint=255
		RGB='#%02x%02x%02x' % (Rint, Gint, Bint)
	elif(value>intervalSize and value<2*intervalSize):
		Rint=200+((50/intervalSize)*(value-intervalSize))
		Gint=Rint
		Bint=250-(200/intervalSize)*(value-intervalSize)
		RGB='#%02x%02x%02x' % (Rint, Gint, Bint)
	elif (value>=2*intervalSize):
		Rint=250
		Gint=250-(200/intervalSize)*(value-2*intervalSize)
		Bint=50
		RGB='#%02x%02x%02x' % (Rint, Gint, Bint)
	return RGB

def getGreater(values):
	great=0
	for vec in values:
		N=len(vec)
		for i in range(0,N-1):
			val=vec[i]
			if(val!="\n" and val!="\0"):
				tval=float(val)
				if (tval>great):
					great=tval
	return great









#main interface exe
class mainGUI(tk.Tk):
	def __init__(self):
		tk.Tk.__init__(self)
		self.geometry("400x250")
		self.configure(bg="#fff")
		self.instruction=tk.Entry(self,font=("Helvetica",12),width=40)
		self.labelInst=tk.Label(self,text="Digite la instrucción aquí",
			bg="#fff",font=("Helvetica",24),fg="#32467D")
		self.labelInst.grid(row=1,padx=20)
		self.labelInst.configure()
		self.instruction.grid(row=2,pady=15,padx=15)
		self.startButton=tk.Button(self,bg="#32467D",fg="#fff",
			font=("Helvetica",30),text="EJECUTAR",
			activebackground="#0096FF",activeforeground="#fff",
			command=lambda: self.startButtonFunction(20))
		self.startButton.grid(pady=20)
		

	def startButtonFunction(self,sF):#sF:size factor
		strLine=self.instruction.get()
		strList=strLine.split(" ")
		subprocess.call(["chmod", "+x", "scriptPython"])	
		subprocess.call("./scriptPython")
		wd=os.getcwd()
		os.chdir(wd+"/build")
		strList[0]="./src/"+strList[0]
		print(strList[0])
		print(strList[1])
		print(strList[2])
		subprocess.call(strList)
		os.chdir(wd)
		

		values=[]
		f=open("Matrix.txt","r")

		f2=open("flags.txt","r")
		flags=f2.read().split("\t")

		if (flags[0]=="True"):
			lines=f.readlines()
			for line in lines:
				instructLst=line.split("\t")
				values.append(instructLst)


			ancho=len(values[0])
			alto=len(values)
			otherWindow=tk.Tk()
			otherWindow.geometry("{0}x{1}".format(ancho*sF,alto*sF))
			

			canv=tk.Canvas(otherWindow,bg="#fff",height=ancho*sF,width=alto*sF)
			for x1 in range(0,len(values)):
				for x2 in range(0,len(values[0])-1):
					coord=x2*sF,x1*sF,x2*sF+sF,x1*sF+sF
					color=getColor(float(values[x1][x2]),getGreater(values))
					canv.create_rectangle(coord,fill=color,outline=color)
					
			canv.pack()
			otherWindow.mainloop()

app=mainGUI()
app.mainloop()

	


# import Tkinter
# import tkMessageBox

# top = Tkinter.Tk()

# C = Tkinter.Canvas(top, bg="blue", height=250, width=300)
# for i in range(0,50,10):
# 	coord = i, i, i+10, i+10
# 	C.create_rectangle(coord,fill="black",outline="black")

# C.pack()
# top.mainloop()
