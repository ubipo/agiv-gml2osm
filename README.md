# agiv-gml2osm
Python script to convert AGIV's gml files to the osm format.

## Installation
```
$ pip install -r requirements.txt
```

## Usage
```
$ ./agiv-gml2osm.py -h
usage: agiv-gml2osm.py [-h] [-v VERBOSE] src [des]

Convert AGIV GRB gml XML files with lambert72 coordinates to the osm XML
format.

positional arguments:
  src                   the directory to scan or .gml file to convert
  des                   the .osm file to write to

optional arguments:
  -h, --help            show this help message and exit
  -v VERBOSE, --verbose VERBOSE
                        run in verbose (debug) mode
```
