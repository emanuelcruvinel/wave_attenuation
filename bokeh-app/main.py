from os.path import join, dirname
import numpy as np
 
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure
import os
 
#os.chdir(r'C:\Users\user\Desktop\POSMOT\DISSERTACAO\ARTIGO_PERIODIC_BEAM')
 
 
def beam(lC1, lC2, wi, hC1, hC2):
 
    Ia = wi*hC1**3/12
    Ib = wi*hC2**3/12
    Aa = wi*hC1
    Ab = wi*hC2
    E = 69e9
    La = lC1/5
    Lb = lC2/5
    rho = 2700
    
    N = 1000
    omega = np.linspace(0,50000,N)
    #freq = omega/(2*np.pi)
    
    def periodic_beam(w, E, I, L, A, rho):
    
      K =  (E*I/L**3)* np.array([[12,6*L,-12,6*L],[6*L,4*L**2,-6*L,2*L**2],[-12,-6*L,12,-6*L],[6*L,2*L**2,-6*L,4*L**2]])
    
      M = (rho*A*L/420)*np.array([[156,22*L,54,-13*L],[22*L,4*L**2,13*L,-3*L**2],[54,13*L,156,-22*L],[-13*L,-3*L**2,-22*L,4*L**2]])
    
      D = K - w**2*M
      Dlr = np.array([[ D[0][2], D[0][3]], [D[1][2], D[1][3]]])
      Drl = np.array([[ D[2][0], D[2][1]], [D[3][0], D[3][1]]])
      Dll = np.array([[ D[0][0], D[0][1]], [D[1][0], D[1][1]]])
      Drr = np.array([[ D[2][2], D[2][3]], [D[3][2], D[3][3]]])
      T1 = np.matmul(-np.linalg.inv(Dlr) ,Dll)
      T2 = np.linalg.inv(Dlr)
      T3 = np.matmul(np.matmul(Drr ,np.linalg.inv(Dlr)) , Dll)  -  Drl
      T4 = np.matmul(-Drr , np.linalg.inv(Dlr))
      T = np.zeros(shape=(4,4), dtype=complex)
      T[0:2, 0:2] = T1
      T[0:2, 2:4] = T2
      T[2:4, 0:2] = T3
      T[2:4, 2:4] = T4
      
      return T
    
    final = np.zeros(shape=(N,4), dtype=complex)
    i = 0
    for w in omega:
      eq1= periodic_beam(w, E, Ia, La, Aa, rho) 
      eq2 = periodic_beam(w, E ,Ib, Lb, Ab,rho)
      Tfinal = eq1@eq1@ eq1 @ eq1 @ eq1@ eq2 @ eq2@eq2@eq2@eq2
      eig_final = np.linalg.eigvals(Tfinal)
      eig_final = np.sort(eig_final)
      final[i] = eig_final
      i += 1
    
    mi = np.arccosh((final[:,1] + 1/final[:,1])/2)
    delta = np.real(mi)
    epsilon = np.abs(np.imag(mi))
    
    return delta, epsilon
 
 
# Set up data
lengthC1 = 0.05
lengthC2 = 0.015
width = 0.05
heightC1 = 0.001
heightC2 = 0.003
 
 
N = 1000
omega = np.linspace(0,50000,N)
x = omega/(2*np.pi)
# x = np.linspace(0, 4*np.pi, N)
# y = np.sin(x)
y = beam(lengthC1, lengthC2, width, heightC1, heightC2)[0]
source = ColumnDataSource(data=dict(x=x, y=y))
 
 
# Set up plot
plot = figure(plot_height=400, plot_width=400, title="Attenuation factor",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 8000], y_range=[0, 0.8])
 
plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
 
 
# Set up widgets
text = TextInput(title="title", value='Attenuation factor')
lc1 = Slider(title="lc1", value=5, start=4.75, end=5.25, step=0.0005)
lc2 = Slider(title="lc2", value=1.5, start=1, end=5.25, step=0.0005)
wi = Slider(title="wi", value=5, start=3, end=6, step=0.0005)
hc1 = Slider(title="hc1", value=1, start=0.8, end=5, step=0.0005)
hc2 = Slider(title="hc2", value=3, start=0.8, end=5, step=0.0005)
 
 
 
# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value
 
text.on_change('value', update_title)
 
def update_data(attrname, old, new):
 
    # Get the current slider values
    l1 = lc1.value
    l2 = lc2.value
    w = wi.value
    h1 = hc1.value
    h2= hc2.value
 
 
    # Generate the new curve
    x = omega/(2*np.pi)
    y = beam(l1*0.01, l2*0.01, w*0.01, h1*0.001, h2*0.001)[0]
 
    source.data = dict(x=x, y=y)
 
for w in [lc1, lc2, wi, hc1, hc2]:
    w.on_change('value', update_data)
 
 
# Set up layouts and add to document
inputs = column(text,lc1, lc2, wi, hc1, hc2)
 
curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Wave Attenuation Periodic Beam"
