# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 19:26:46 2021

@author: elekp
"""
from PIL import Image

from datetime import datetime
from math import log

import struct


now = datetime.now() # current date and time

im=[]

num_of_frames=100

for num in range(num_of_frames):
    filename="./gif_material/kep_double"+str(num)+".dat"#"data"+str(num//100)+str((num%100)//10)+str((num%100)%10)+".dat"#_double
    print(filename)
    file = open(filename, 'rb')
    a=int.from_bytes(file.read(4),byteorder='little', signed=True)
    b=int.from_bytes(file.read(4),byteorder='little', signed=True)

    print(a)
    print(b)



    xmax=0.0
    img = Image.new( 'RGB', (b,a), "black") # Create a new black image
    pixels = img.load() # Create the pixel map
    for i in range(img.size[0]):    # For every pixel:
        for j in range(img.size[1]):
            x=struct.unpack('f', file.read(4))[0]
            if(xmax<x):
                xmax=x
            if (x==0):
                pixels[i,j] = (0, 0, 0)
            elif(x==-1):
                pixels[i,j] = (0, 0, 0)#(255, 0, 0)#(255, 255, 255)#
            else:
                
                T=x/0.5*256*4
                r=1.1
                if((T//256)%3==0):
                    T=T%256
                    pixels[i,j] = (int(255-T/r), 0, 0)
                elif((T//256)%3==1):
                    T=T%256
                    pixels[i,j] = (0, int(255-T/r), 0)
                elif((T//256)%3==2):
                    T=T%256
                    pixels[i,j] = (0, 0, int(255-T/r))
                #pixels[i,j] = (0, 255, 0) ha 1 szinure akarjuk
                    
                
    im.append(img)
    
    #nevex="blackholle"+now.strftime("%m_%d_%Y_%H_%M_%S")+"__"+str(num)
    #windows#"C:\\Users\\elekp\\source\\repos\\kerr_class_szines_kep\\kerr_class_szines_kep\\kepekx\\blackholle"+now.strftime("%m_%d_%Y_%H_%M_%S")+"__"+str(num)
    #img.save(nevex+".png")
        
print(xmax)
#img.show()
neve="blackhole"
#"C:\\Users\\elekp\\source\\repos\\kerr_class_szines_kep\\kerr_class_szines_kep\\kepekx\\blackholle"+now.strftime("%m_%d_%Y_%H_%M_%S")

agif=Image.new( 'RGB', (b,a), "black")
agif=im[0]

for num in range(num_of_frames):
    im.append(im[num_of_frames-num-1])



agif.save(neve+".gif",save_all=True, append_images=im,loop=0)

print("kesz")
