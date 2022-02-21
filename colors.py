# get colors for the plots

import matplotlib.pyplot as plt
import numpy as np

# make viridis available
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='magma', cmap=cmaps.magma)


def get_plasma(num):
    plasma_custom = {}
    for i in range(num):
        plasma_custom[i] = cmaps.plasma(0+i*(int(256./float(num-1))))[:3]
    return plasma_custom

def get_inferno(num):
    inferno_custom = {}
    for i in range(num):
        inferno_custom[i] = cmaps.inferno(0+i*(int(256./float(num-1))))[:3]
    return inferno_custom

def get_viridis(num):
    viridis_custom = {}
    for i in range(num):
        viridis_custom[i] = cmaps.viridis(0+i*(int(256./float(num-1))))[:3]
    return viridis_custom

def get_magma(num):
    magma_custom = {}
    for i in range(num):
        magma_custom[i] = cmaps.magma(0+i*(int(256./float(num-1))))[:3]
    return magma_custom

num_colours = 5 # how many colours do you want 

viridis_green = {}
for i in range(num_colours):
    viridis_green[i] = cmaps.viridis(70+i*37)[:3]
viridis_purple = {}
for i in range(num_colours):
    viridis_purple[i] = cmaps.viridis(0+i*40)[:3]
viridis_all = {}
for i in range(num_colours):
    viridis_all[i] = cmaps.viridis(0+i*64)[:3]
viridis_more = {}
for i in range(10):
    viridis_more[i] = cmaps.viridis(0+i*(int(256./10.)))[:3]

plasma_pink = {}
for i in range(num_colours):
    plasma_pink[i] = cmaps.plasma(50+i*30)[:3]
plasma_orange = {}
for i in range(num_colours):
    plasma_orange[i] = cmaps.plasma(140+i*20)[:3]
plasma_all = {}
for i in range(num_colours):
    plasma_all[i] = cmaps.plasma(0+i*64)[:3]
plasma_more = {}
for i in range(10):
    plasma_more[i] = cmaps.plasma(0+i*(256/10))[:3]

inferno_all = {}
for i in range(num_colours):
    inferno_all[i] = cmaps.inferno(20+i*55)[:3]

magma_pink = {}
for i in range(num_colours):
    magma_pink[i] = cmaps.magma(50+i*30)[:3]
magma_orange = {}
for i in range(num_colours):
    magma_orange[i] = cmaps.magma(140+i*20)[:3]
magma_all = {}
for i in range(num_colours):
    magma_all[i] = cmaps.magma(20+i*55)[:3]
magma_more = {}
for i in range(10):
    magma_more[i] = cmaps.magma(0+i*(256/10))[:3]


standard_green = viridis_green[2]
standard_red   = magma_orange[2]


def rgb2hex(rgb):
    r,g,b = rgb
    r=int(r*256)
    g=int(g*256)
    b=int(b*256)
    hex = "#{:02x}{:02x}{:02x}".format(r,g,b)
    return hex


def change_color_names(color):
    # change format of colors from RGB to Hex
    try:
        col = '#%02x%02x%02x' % tuple(int(255*x) for x in color)
    except:
        col = '#%02x%02x%02x%02x' % tuple(int(255*x) for x in color)
    return col

