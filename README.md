# jwst_planning_tools
Convenience functions to plan JWST observations

## Install 
Clone the repository
```
git clone https://github.com/patapisp/jwst_planning_tools.git
```
Dependencies
- webbpsf development branch: https://github.com/spacetelescope/webbpsf.git
- miricoord package: https://github.com/STScI-MIRI/miricoord

For miricoord you need to follow the instructions on the gitgub repo. In your python environment, clone the repository 
and run
```commandline
python setup.py install
```
Then continue the jwst planning tools installation and run
```
cd jwst_planning_tools/
pip install -e .
```
Be sure to also point to the webbpsf data by downloading them from https://webbpsf.readthedocs.io/en/stable/installation.html#installing-the-required-data-files
and exporting the "WEBBPSF_PATH" environment variable.