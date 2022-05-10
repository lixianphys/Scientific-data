# Scientific-data (SciData)

This package serves as a tool box to prepare raw data for later plotting
in matplotlib, largely smoothing the process from reading raw data to generating
processed data then to plotting them. 

Beyond general data handling, we also offer a toy-band model to calculate Landau
levels in energy with a changing magnetic field in multiple-band material system
(user-defined). Check it out [Here](docs/readme.md).

### Introduction to ```SciData.py```
In this ```.py``` file, we define all the classes for structing raw data and pre-processing
them. In many experiments, we take data for a changing parameter A (sweep A) at a fixed
parameter B in one ```.dat``` file, and then change the parameter B and sweep A. All the data
will be stored in a folder with files named by the specific parameter B. Then in this case,
our defined class can turn all the data in the folder into an instance of this class.
By adding methods to this class, we can pre-process the data (```getdata(self)```), do Hall
fit if one dimension is the magnetic field (```hallfit(self,fitrange)```) and quick plot data 
(```plotdata(self,label_value)```). The functionality of parent class ```Datajungle``` can be 
extended by creating more child class.

| Child class | parameter A | parameter B|
|-------------|------------|-----------|
| Databs      |magnetic field|gate-voltage|
| Datags      | gate-voltage|magnetic field|
| Datafc      | gate-voltage|gate-voltage|
| DataX| unknown|unknown|

The first three are most useful when handling (quantum) Hall measurement. For general flexibility,
we also provide a less-defined class ```DataX``` to handle many other types of data.

### How to use this package to speed up your work flow
```startnb.py``` is a script to import default setting of your ```jupyter notebook```, especially
importing all the packages, scripts and functions to be ready to use and also
plotting setting under ```plt.rc```. 
```python
%run [dir to SciData]/startnb.py
```
Add this line in your first cell of opened notebook to start your data processing.