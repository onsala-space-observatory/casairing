# CASAIRING

## Installation

Steps to install the `casairing` task into `casa`

 1. Clone the git repository into a directory of your choice
 (e.g. `$HOME/.casa/NordicTools`):

``` shell
cd $HOME/.casa/NordicTools
git clone <repository url>
cd casairing
buildmytasks --module casairing casairing.xml
```
 2. Inside `casa` add the folder to your `PYTHONPATH`:

``` python
CASA <1>: sys.path.insert(0, os.getenv("HOME") + "/.casa/NordicTools")
CASA <2>: from casairing.gotasks.casairing import casairing
CASA <3>: inp(casairing)

```

That's it! Enjoy!
