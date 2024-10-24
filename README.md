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
You may need to remove the `__init__.py` file from the top level miricoord directory.
Be sure to also point to the webbpsf data by downloading them from https://webbpsf.readthedocs.io/en/stable/installation.html#installing-the-required-data-files and exporting the "WEBBPSF_PATH" environment variable.

Then continue the jwst planning tools installation and run
```
cd jwst_planning_tools/
pip install -e .
```
You may need to add `jwst_planning_tools` to your python path, and add this to your bashrc file:
```
export PYTHONPATH=/path/to/jwst_planning_tools/:$PYTHONPATH
```
## Example
![HR8799](https://github.com/patapisp/jwst_planning_tools/blob/main/_HR%208799_1A_V3_75.0_offsetXYidl_%5B%200.07%20-0.49%5D.png)
