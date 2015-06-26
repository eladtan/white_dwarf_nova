import os

env = Environment(ENV = os.environ,
                  CXX='clang++',
                  #CXX='g++',
                  CPPPATH=[os.environ['RICH_ROOT']+'/source',
                           os.environ['RICH_ROOT']],
                  LIBPATH=[os.environ['RICH_ROOT'],'.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp','gfortran'],
                  #LINKFLAGS=' -g -Og ',
                  LINKFLAGS='',
                  #F90FLAGS=' -g -Og -ffpe-trap=invalid,overflow,underflow,zero,denormal -ffpe-summary=all ',
                  F90FLAGS=' -O2 ',
                  #CXXFLAGS='-Weverything -Werror -ferror-limit=1 -Wno-error=padded -g -O0 ')
                  #CXXFLAGS=' -g -Og -Wfatal-errors ')
                  CXXFLAGS=' -Weverything -Werror -ferror-limit=1 -Wno-error=padded -O2 ')
                  
env.Program(['nuclear_burn.cpp','rich.cpp','fermi_table.cpp','eos_fermi.f90','tables_module.f90',
             'rho_enr.f90','rho_tmp.f90','rho_prs.f90','tmp_prs.f90','enr_prs.f90',
             'burn_step.f90','nse.f90','nse_equilib.f90','net_step.f90',
             'sigmav_rates.f90','alpha_rates.f90','linsys_burn.f90','sigmav.f90',
             'screening.f90'])
