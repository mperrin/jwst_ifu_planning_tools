# For debugging: set to debug to 'true' and optimization to 'g'
project('jwst_planning_tools',
    'c',
    # version must be updated here as well
    version : '0.0.2',
    license: 'MIT',
    meson_version: '>=1.1.0',
    default_options : [
        'debug=false',
        'optimization=3',
        'warning_level=1',
        'c_std=c99'
    ],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Install python sources
subdir('jwst_planning_tools')
