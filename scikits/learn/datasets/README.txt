Last Change: Tue Jul 17 04:00 PM 2007 J

This packages datasets defines a set of packages which contain datasets useful
for demo, examples, etc... This can be seen as an equivalent of the R dataset
package, but for python.

Each subdir is a python package, and should define the function load, returning
the corresponding data. For example, to access datasets data1, you should be able to do:

>> from datasets.data1 import load
>> d = load() # -> d contains the data of the datasets data1

load can do whatever it wants: fetching data from a file (python script, csv
file, etc...), from the internet, etc... Some special variables must be defined
for each package, containing a python string:
    - COPYRIGHT: copyright informations
    - SOURCE: where the data are coming from
    - DESCHOSRT: short description
    - DESCLONG: long description
    - NOTE: some notes on the datasets.

For the datasets to be useful in the learn scikits, which is the project which initiated this datasets package, the data returned by load has to be a dict with the following conventions:
    - 'data': this value should be a record array containing the actual data.
    - 'label': this value should be a rank 1 array of integers, contains the
      label index for each sample, that is label[i] should be the label index
      of data[i].
    - 'class': a record array such as class[i] is the class name. In other
      words, this makes the correspondance label index <> label name.
