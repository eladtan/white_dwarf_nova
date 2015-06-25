import plotly.plotly as py
from plotly.graph_objs import *
import pickle
import os
import numpy

if not os.path.isfile('init_cond_density.pkl'):
    print 'retrieving density data'
    pickle.dump(py.get_figure('https://plot.ly/~almog.yalin/190').get_data(),
                open('init_cond_density.pkl','wb'))
if not os.path.isfile('radius_list.txt'):
    density_data = pickle.load(open('init_cond_density.pkl','rb'))
    print 'writing radius_list.txt'
    numpy.savetxt('radius_list.txt',density_data[0]['x'])
if not os.path.isfile('density_list.txt'):
    density_data = pickle.load(open('init_cond_density.pkl','rb'))
    print 'writing density list'
    density_data = pickle.load(open('init_cond_density.pkl','rb'))
    numpy.savetxt('density_list.txt',density_data[0]['y'])

if not os.path.isfile('init_cond_temperature.pkl'):
    print 'retrieving temperature data'
    pickle.dump(py.get_figure('https://plot.ly/~almog.yalin/194').get_data(),
                open('init_cond_temperature.pkl','wb'))
if not os.path.isfile('temperature_list.txt'):
    print 'writing temperature list'
    temperature_data = pickle.load(open('init_cond_temperature.pkl','rb'))
    numpy.savetxt('temperature_list.txt',temperature_data[0]['y'])

if not os.path.isfile('init_cond_pressure.pkl'):
    print 'retrieviing pressure data'
    pickle.dump(py.get_figure('https://plot.ly/~almog.yalin/192').get_data(),
                open('init_cond_pressure.pkl','wb'))
if not os.path.isfile('pressure_list.txt'):
    print 'writing pressure list'
    pressure_data = pickle.load(open('init_cond_pressure.pkl','rb'))
    numpy.savetxt('pressure_list.txt',pressure_data[0]['y'])        

if not os.path.isfile('init_cond_velocity.pkl'):
    print 'retrieving velocity data'
    pickle.dump(py.get_figure('https://plot.ly/~almog.yalin/196').get_data(),
                open('init_cond_velocity.pkl','wb'))
if not os.path.isfile('velocity_list.txt'):
    print 'writing velocity list'
    velocity_data = pickle.load(open('init_cond_velocity.pkl','rb'))
    numpy.savetxt('velocity_list.txt',velocity_data[0]['y'])

if not os.path.isfile('init_cond_composition.pkl'):
    print 'retrieving composition data'
    pickle.dump(py.get_figure('https://plot.ly/~almog.yalin/188').get_data(),
                open('init_cond_composition.pkl','wb'))
composition_data = pickle.load(open('init_cond_composition.pkl','rb'))
for datum in composition_data:
    text_file_name = 'tracer_'+datum['name']+'.txt'
    if not os.path.isfile(text_file_name):
        print 'writing '+text_file_name
        numpy.savetxt(text_file_name,datum['y'])
