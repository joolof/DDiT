# DDiT
## Debris DIsks Tool

![screenshot](screenshots/DDiT.png)

A Python code to (quickly) produce synthetic images of debris disks, in total and polarized intensity using the Henyey-Greenstein approximation. The disk model is not infinitely flat and in principle should work for any inclinations. At the moment, there are some weird issues for large inclinations, and this is on the todo list.

To use the module, download the python files above and install it with (using the `develop` option so that any further changes will be included).
```python
python3 setup.py develop
```

Given that the code only uses matplotlib and numpy, it works with python2 and python3.




```python
from DDiT import Disk
```

