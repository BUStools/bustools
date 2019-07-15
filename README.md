# bustools

__bustools__ is a program for manipulating [__BUS__](https://github.com/BUStools/BUS) files for single cell 
RNA-Seq datasets. It can be used to error correct barcodes, collapse UMIs, produce gene count or transcript compatbility count matrices, and is useful for many other tasks. See the [__kallisto &#124; bustools website__](https://www.kallistobus.tools/) for examples and instructions on how to use __bustools__ as part of a single-cell RNA-seq workflow.

The design and motivation for the __BUS__ format and __bustools__ are described in detail in 

P Melsted, V Ntranos, L Pachter, [The Barcode, UMI, Set format and BUStools](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz279/5487510), Bioinformatics, btz279, 2019.


## BUS format

__bustools__ works with __BUS__ files which can be generated efficiently from raw sequencing data, e.g. using [__kallisto__](http://pachterlab.github.io/kallisto).

## Installation

Binaries for Mac, Linux, Windows, and Rock64 can be downloaded from the [__bustools__ website](https://bustools.github.io/download). 

To compile bustools download the source code with

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

The bustools executable will be located in build/src. To install bustools into the cmake install prefix path type:

`make install`

## Usage

To see a list of available commands, type `bustools` in the terminal

~~~
> bustools 
Usage: bustools <CMD> [arguments] ..

Where <CMD> can be one of:

sort            Sort bus file by barcodes and UMI
text            Output as tab separated text file
correct         Error correct bus files
count           Generate count matrices from bus file
capture         Capture reads mapping to a transcript capture list

Running bustools <CMD> without arguments prints usage information for <CMD>
~~~

### Sorting

Raw BUS output from pseudoalignment programs may be unsorted. To simply and accelerate downstream processing BUS files can be sorted using `bustools sort`

~~~
> bustools sort 
Usage: bustools sort [options] bus-files

Options: 
-t, --threads         Number of threads to use
-m, --memory          Maximum memory used
-T, --temp            Location and prefix for temporary files 
                      required if using -p, otherwise defaults to output
-o, --output          File for sorted output
-p, --pipe            Write to standard output
~~~

This will create a new BUS file where the BUS records are sorted by barcode first, UMI second, and equivalence class third.

### Text

BUS files can be converted to a tab-separated format for easy inspection and processing using shell scripts or high level languages with `bustools text`.

~~~
> bustools text
Usage: bustools text [options] bus-files

Options: 
-o, --output          File for text output
~~~

### Correct
BUS files can be barcode error corrected with respect to a technology-specific whitelist of barcodes using `bustools correct`.

~~~
> bustools correct
Usage: bustools correct [options] bus-files

Options: 
-o, --output          File for corrected bus output
-w, --whitelist       File of whitelisted barcodes to correct to
-p, --pipe            Write to standard output
~~~

### Count
BUS files can be converted into a barcode-feature matrix, where the feature can be TCCs (Transcript Compatibility Counts) or genes using `bustools count`.

~~~
> bustools count
Usage: bustools count [options] bus-files

Options: 
-o, --output          File for corrected bus output
-g, --genemap         File for mapping transcripts to genes
-e, --ecmap           File for mapping equivalence classes to transcripts
-t, --txnames         File with names of transcripts
--genecounts          Aggregate counts to genes only
~~~

### Capture
`bustools capture` can separate BUS files into multiple files according to the capture criteria.

~~~
Usage: bustools capture [options] bus-files

Options: 
-o, --output          Directory for output 
-c, --capture         List of transcripts to capture
-e, --ecmap           File for mapping equivalence classes to transcripts
-t, --txnames         File with names of transcripts
~~~


