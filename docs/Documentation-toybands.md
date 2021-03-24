 # Table of Contents
### Project Overview
### Prerequisites
- Install Python package (v3.8.7 recommended)
- Pip install dependencies list in requirements.txt `pip install -r requirements.txt`
### Main tasks
- build a band of arbitrary parameters in E-B relationship
- build a system consisting of bands specified above
- calculate the density of state (DOS) at arbitrary B-field and energy E given the position of band bottom/a specific set of carrier densities. (They are equivalent in our model)
- Based on DOS results, draw the chemical potentail line at each B-field
- plot the E-B relation of Landau levels
- plot the n(density)-B relation of Landau levels. Here we need to use the DOS information to map energy (E) to density (n)
- plot the n-B relation of Landau levels at the chemical potential to show the simulation of experimentally obtained Landau fan chart
### Command-line-interface
Use a teminal (powershell, gitbash for Win user) and switch current directory to `cd ..\SciData>` and start the command with `python` to invoke your python program and the `toybands` module follows. 
An example is `..\SciData> python toybands/xxx.py`.  To access a brief help doc, just add argument `-h` by the end like `python toybands/addband.py -h`. As of now, I will introduce the CLI and its usage, covering common scenarios. Note that you can always interrupt the script by `Ctrl+C` but nothing will be saved after this interruption.

##### module `addband` 
`python toybands/addband.py [-option][parameters]`
- `density`: positional, density for this band in unit 1/m2
- `-is_cond`: optional, conduction (present)/valence (absent) band
- `-gfactor`: optional, gfactor for this band
- `-dp DP1 DP2`: optional, DP1: M (eV) DP2: vf(m/s) for a Dirac-like band. M is the mass term for massive Dirac dispersion. For M=0, it represents for the massless Dirac dispersion
- `-cp CP1 CP2`: optional, CP1: meff (me) CP2: spin (+1/-1) for a conventional band. meff is the effective mass in unit of the rest mass of electron, spin stands for spin-up (+1) and spin-down (-1) cases.
- Initialization of a system:
A system will be initiated at the same time its first band is created.

	
#####  module `peeksys`
Peek into the system you created
 Once you'd like to see what is in your system, how many bands and their parameters, simply type: `python toybands/peaksys.py`.  It will print out a summary of bands in table-like format.
 
 ##### module `delband`
 Remove bands from a system
 This can be done by `python toybands/delband.py`. It will prompt a dialogue `which band to delete? Input the index number: Or press 'e' to exit`. Type the number in front of all the parameters within a row for that band you want to delete. Or quit by typing `e`. 
 
 ##### module `run`
 `python toybands/addband.py [-option][parameters]`
 - `-enplot`: optional, plot the energy versus bfield (yes/no)
 - `-denplot`: optional, plot the density versus bfield (yes/no)
 - `-simu`: optional, dynamically generate relationship between the density and the bfield at steps of input density (yes/no)
 - `-dos`: optional, plot dos versus bfield (yes/no)
 - `--allden`: optional, densities for each band: start1 end1 start2 end2 ...
 - `-nos`: optional, number of steps in the simulation
 - `-dir`: optional, relative output directory
 - `-fnm`: optional, filename
 - `--enrange`: optional, energy range: start end numofpoints
 - `--bfrange`: optional, magnetic field range: start end numofpoints
 - `-nmax`: optional, number of Landau levels involved (default=20)
 - `-angle`: optional, angle in degree made with the sample plane norm by the external field (default=0)
 
 - Plot the E-B relationship for the configured system
 `python toybands/run.py [-enplot] [-dir DIR] [-fnm FNM] [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] [-nmax NMAX] [-angle ANGLE]` If no `-dir` and `-fnm` is specified, the output figure will be stored in `./output/[auto]default.pdf`  
 - Plot the n-B relationship for the configured system
 `python toybands/run.py [-denplot] [-dir DIR] [-fnm FNM] [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] [-nmax NMAX] [-angle ANGLE]`
 - plot the n-B relationship for a set of given densities as a simulation of Landau fan chart.
 `python toybands/run.py [-simu] [--allden "NS1 NE1 NS2 NE2 ..."] [-nos NOS] [-dir DIR] [-fnm FNM] [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] [-nmax NMAX] [-angle ANGLE]`
 - plot the DOS-B relationship for the configured system
 `python toybands/run.py [-dos] [-fnm FNM] [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum]`
 ### Features in the future
 - plot the DOS mapping onto (n,B) for a set of given densities [***]
 - Input customized E-B relationship [*]
 - label the output figure and allow separate figures in a single pdf file for a batch input [***]
 - Is it possible to abort the calculation until I continue it? I mean, keep the memory usage but free its CPU usage. It might be helpful when you just want your PC to handle other CPU-heavy tasks during a long calculation. [*]
 - plot the data from csv file directly from the CLI [**]
 - allow multiple stage input of densities in simu[**]
 - connect the same LL in simu [**]
