import os

debug = ARGUMENTS.get('debug',0)
compiler = ARGUMENTS.get('compiler','clang++')

linkflags = ''
if compiler=='g++':
    cflags = ' -Wfatal-errors '
    if int(debug):
        cflags += ' -Og -g -pg '
        linkflags = ' -g -pg '
    else:
        cflags += ' -O3 '
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg'
        linkflags = ' -g -pg'
    else:
        cflags += ' -O3 -march=native'
else:
    raise NameError('unsupported compiler')

if int(debug):
    f90flags = ' -g -Og '
else:
    f90flags = ' -O3 '

env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=[os.environ['RICH_ROOT']+'/source',
                           os.environ['RICH_ROOT']],
                  LIBPATH=[os.environ['RICH_ROOT'],'.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp','gfortran'],
                  LINKFLAGS=linkflags,
                  F90FLAGS=f90flags,
                  CXXFLAGS=cflags)
                  
env.Program('rich',
            ['nuclear_burn.cpp','rich.cpp','fermi_table.cpp','eos_fermi.f90','tables_module.f90',
             'rho_enr.f90','rho_tmp.f90','rho_prs.f90','tmp_prs.f90','enr_prs.f90',
             'burn_step.f90','nse.f90','nse_equilib.f90','net_step.f90',
             'sigmav_rates.f90','alpha_rates.f90','linsys_burn.f90','sigmav.f90',
             'screening.f90','initnet.f90'])
