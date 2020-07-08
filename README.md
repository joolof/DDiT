# DDiT
## Debris DIsks Tool

![screenshot](screenshots/DDiT.png)

A Python code to (quickly) produce synthetic images of debris disks, in total and polarized intensity using the Henyey-Greenstein approximation. In the model, the disk is not infinitely flat and while in principle it should work for any inclinations, at the moment, there are some *very* weird issues for large inclinations, and needs to be fixed as soon as possible.

To use the module, download the python files above and install it with the following (using the `develop` option so that any further changes are automatically included):
```python
python3 setup.py develop
```

Given that the code only uses matplotlib and numpy, it should work with python2 and python3. To compute a model, you can have a look at the end of the "DDiT.py" file, but in a nutshell, it is as simple as

```python
from DDiT import Disk
disk = Disk()
disk.compute_model(a = 0.6, i = 60., e = 0.3)
```
Once the model is computed, you will find the images in the variables `disk.intensity` and `disk.polarized` for the total and polarized intensity.

## References

If you use this code for your research, feel free to contact me if you have any doubt or question, and please cite the following paper:
```
Olofsson et al. (2020)
The challenge of measuring the phase function of debris disks. Application to HR\,4796
```

