Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

Go back

> cd -

Tkinter must be installed beforehand

> sudo apt install python-tk

Open Python and run the project in the file location

> python GUI.py

Then type in the textfield all the conditions 

Example:

> placa -ib -t30 -p perfil1.txt

The standard for the edges conditions is

> [Val[0], Val[n]]

Example:

> [25, 100]

Where:

> Val[0] = 25 and Val[n] = 100

Yo can try the test scrips:

> scriptPrueba1

Or

> scriptPrueba2

Or
> scriptPrueba3