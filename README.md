# MastersThesis

This is the repository of the specific code used for my Master's thesis.\
The bulk was done in python3, while the actual report was written in $\LaTeX2e$.\

A few enviroment variables are defined in the virtual enviroment's activate script, although they could also be defined in the computer itself.\
Under the directory 'Code' is where the different scripts live, sorted by programing language. 
Python has two directories for itself, one storing the executable scripts, while the other stores the packages I made for usability.
Julia was hoped to play a bigger role early on, but I've been remaking the simulation packages in the form of Julia modules.
Nothing can beat the versatility of matplotlib however, so all results are stored in `.npz` files to be opened and plotted by 
`Code/python/plotter.py` using my plotting package for consistency.

If there is one thing I have learned, if nothing else from all this, it is that I hate object oriented programing.\
And that I wish I commited to Julia earlier on, rather than prototyping so heavily in python that it just did not make sense to 
try and recreate everything in Julia, within the allocated timeline that is.
