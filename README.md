# bustools

__bustools__ is a program for manipulating [__BUS__](https://github.com/BUStools/BUS) files for single cell 
RNA-Seq datasets. 


The design and motivation for the BUS format and BUStools are described in detail in 

P Melsted, V Ntranos, L Pachter, [The Barcode, UMI, Set format and BUStools](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz279/5487510), Bioinformatics, btz279, 2019.


## BUS format

__bustools__ works with __BUS__ files which can be generated efficiently from raw sequencing data, e.g. using [__kallisto__](http://pachterlab.github.io/kallisto).

## Installation

Download bustools with

`git clone https://github.com/BUStools/bustools.git`

Navigate to the bustools directory

`cd bustools`

Make a build directory and move there:

`mkdir build`

`cd build`

Run cmake:

`cmake ..`

Build the code:

`make`

The bustools executable is now located in build/src. To install bustools into the cmake install prefix path type:

`make install`

## Usage

To see a list of available commands type `bustools` in the terminal

~~~
> bustools 
Usage: bustools <CMD> [arguments] ..

Where <CMD> can be one of:

sort            Sort bus file by barcodes and UMI
text            Output as tab separated text file

Running bustools <CMD> without arguments prints usage information for <CMD>
~~~

### Sorting

Raw BUS output from pseudoalignment programs may be unsorted. To simply and accelerate downstream processing BUS files can be sorted using `bustools sort`

~~~
> bustools sort 
Usage: bustools sort [options] bus-files

Options:
-t, --threads         Number of threads to use
-o, --output          File for sorted output
~~~

This will create a new BUS file where the BUS records are sorted by barcode first, UMI second, and equivalence class third.

### Text

BUS files can be converted to a tab-separated format for easy inspection and processing using shell scripts or high level languages. `bustools text` 

~~~
> bustools text
Usage: bustools text [options] bus-files

Options: 
-o, --output          File for text output
~~~
