 # Table of Contents
[Project Overview](#project-overview)

[Prerequisites](#prerequisites)

[Main tasks](#main-tasks)

[Command-line-interface CLI](#command-line-interface-(CLI))

[Usage Examples below](#usage-examples-below)

[Module introduction](#module-introduction)

[Features in the future](#features-in-the-future)

### Project Overview
This toybands model allows one to get some fundamentals of quantum Hall effect and Landau levels in high magnetic field if you happens to have some Dirac or non-Dirac materals, which trends in condensed matter physics in the past ten years or so. 
You can use some parameters as inputs to describe the properties of electron system of your materials, usually in a simplified form (Fermi velocity, spin projection in the magnetic field direction, effective mass, etc).  In the end, you are enabled to carry out multiple tasks with this model, like plotting the Landau levels in energy or monitoring the effect of density change in Landau levels as the magnetic field ramps. The output will be stored in a separate folder for later check-out. 
All the physics behind is explained in a separate documentation  [[Physics here]]. 

### Prerequisites
- Install Python package (v3.8.7 recommended)
- Pip install dependencies list in requirements.txt `pip install -r requirements.txt`
### Main tasks
- build a band with preset parameters and add it to an existing system (`addband.py`)
- build a system consisting of bands specified above , the system can contain arbitrary number of bands with arbitrary parameters. (`addband.py`)
- delete band(s) from a system. (`delband.py`)
- snapshot the system (`peeksys.py`)
- do some calculations on the system (`run.py`)
	- calculate the density of state (DOS) at arbitrary B-field and energy E given the position of band bottom/a specific set of carrier densities. (They are equivalent in our model)
	- Based on DOS results, draw the chemical potentail line at each B-field
	- plot the E(energy)-B(field) relation of Landau levels
	- plot the n(density)-B relation of Landau levels. Here we need to use the DOS information to map energy (E) to density (n)
	- plot the n-B relation of Landau levels at the chemical potential to show the simulation of experimentally obtained Landau fan chart
	- map the density of state at each (B,n)
- replot the results stored in csv file(`plotfile.py`)
### Command-line-interface (CLI)
- If you are the first time to use a CLI:
	Use a teminal (powershell, gitbash for Win user) and switch current directory to `cd ..\SciData>` and start the command with `python` to invoke your python program and the `toybands` module follows. 
	An example is `..\SciData> python toybands/xxx.py`.  To access a brief help doc, just add argument `-h` by the end like `python toybands/addband.py -h`. As of now, I will introduce the CLI and its usage, covering common scenarios. Note that you can always interrupt the script by `Ctrl+C` but nothing will be saved after this interruption.
- else: Go to `SciData\` directory.
 ### Usage (Examples below)
 Note that `[]` means a "must", `{}` means optional, and `[{case1},{case2}]` means one of them must exist. All the arguments can be in any order within the command line. 
 - add a band to the system
 `python toybands/addband.py [-density] {-is_dirac} {-is_cond} [-gfactor GFACTOR] [{-dp MASS VF} {-cp MEFF SPIN}]` 
 - delete band(s) from the system without UI (useful when you run a script)
 `python toybands/delband.py [-i INDEX or all]`
 - have a peek into the system
 `python tobybands/peeksys.py`
 - plot the E-B relationship of Landau levels
 `python toybands/run.py [-enplot] {-dir DIR} {-fnm FNM} [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] {-nmax NMAX} {-angle ANGLE}` If no `-dir` and `-fnm` is specified, the output figure will be stored in `./output/[auto]default.pdf`  
 - plot the n-B relationship of Landau levels
 `python toybands/run.py [-denplot] {-dir DIR} {-fnm FNM} [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] {-nmax NMAX} {-angle ANGLE}`
 - plot the n-B relationship for a set of given densities (specified by  `--allden` and `-nos`) as a simulation of Landau fan chart.
 `python toybands/run.py [-simu] [--allden "NS1 NE1 NS2 NE2 ..."] [-nos NOS] {-dir DIR} {-fnm FNM} [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] {-nmax NMAX} {-angle ANGLE}`
 	Alternatively, one can load densitites stored in columns of a csv file (no header) elsewhere like:
 `python toybands/run.py [-simu] [--loadden path-to-csvfile] {-dir DIR} {-fnm FNM} [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] {-nmax NMAX} {-angle ANGLE}`
 - plot the DOS-B relationship at a fixed chemical potential
 `python toybands/run.py [-dos] {-dir DIR} {-fnm FNM} [--enrange Estart Eend Enum] [--bfrange Bstart Bend Bnum] {-nmax NMAX} {-angle ANGLE}`
 - plot the DOS mapping onto (n,B) for a set of given densities
 substitute `[-simu]` with `[-dosm]` in simulation.
 - plot the data from csv file directly from the CLI
 `python toybands/plotfile.py [-f path-to-csvfile]`

### Module  introduction
##### module `addband` 
`python toybands/addband.py [-option][parameters]`
- `density`: positional, density for this band in unit 1/m2
- `-is_cond`: optional, conduction (present)/valence (absent) band
- `-gfactor`: optional, gfactor for this band
- `-dp DP1 DP2`: optional, DP1: M (eV) DP2: vf(m/s) for a Dirac-like band. M is the mass term for massive Dirac dispersion. For M=0, it represents for the massless Dirac dispersion
- `-cp CP1 CP2`: optional, CP1: meff (me) CP2: spin (+1/-1/0) for a conventional band. meff is the effective mass in unit of the rest mass of electron, spin stands for spin-up (+1), spin-down (-1) and spinless (0) cases.
- Initialization of a system:
A system will be initiated at the same time its first band is created.

	
#####  module `peeksys`
Peek into the system you created
 Once you'd like to see what is in your system, how many bands and their parameters, simply type: `python toybands/peaksys.py`.  It will print out a summary of bands in table-like format.
 
 ##### module `delband`
 Remove bands from a system
 This can be done by `python toybands/delband.py`. It will prompt a dialogue `which band to delete? Input the index number: Or press 'e' to exit`. Type the number in front of all the parameters within a row for that band you want to delete. Or quit by typing `e`. 
 
 ##### module `run`
 `python toybands/run.py [-option][parameters]`
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

 
 ### Features in the future
 - Input customized E-B relationship [*]
 - allow separate figures in a single pdf file for a batch input [*]
 - Is it possible to abort the calculation until I continue it? I mean, keep the memory usage but free its CPU usage. It might be helpful when you just want your PC to handle other CPU-heavy tasks during a long calculation. [*]
 - connect the same LL in simu [**] realize via the output .csv to plot lines with the same band, N at once.
