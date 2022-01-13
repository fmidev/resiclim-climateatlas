#!/usr/bin/env python





import matplotlib.pyplot as plt
import allas_utils as allas


files = allas.ls_s3('resiclim/2m_temperature/')



for file in files:
    
    ds = allas.read_nc_s3(file)
    ds['t2m'].mean('time').plot(); plt.show()
