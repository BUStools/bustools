# bustools

__bustools__ is a program for manipulating [__BUS__](https://github.com/BUStools/BUS) files for single cell 
RNA-Seq datasets. It can be used to error correct barcodes, collapse UMIs, produce gene count or transcript compatibility count matrices, and is useful for many other tasks. See the [__kallisto &#124; bustools website__](https://www.kallistobus.tools/) for examples and instructions on how to use __bustools__ as part of a single-cell RNA-seq workflow.

If you use __bustools__ please cite

Melsted, Páll, Booeshaghi, A. Sina et al. [Modular and efficient pre-processing of single-cell RNA-seq.](https://www.biorxiv.org/content/10.1101/673285v2) BioRxiv (2019): 673285, doi.org/10.1101/673285.

For some background on the design and motivation for the __BUS__ format and __bustools__ see 

Melsted, Páll, Ntranos, Vasilis and Pachter, Lior [The Barcode, UMI, Set format and BUStools](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz279/5487510), Bioinformatics, btz279, 2019.


## BUS format

__bustools__ works with __BUS__ files which can be generated efficiently from raw sequencing data, e.g. using [__kallisto__](http://pachterlab.github.io/kallisto).

## Installation

Binaries for Mac, Linux, Windows, and Rock64 can be downloaded from the [__bustools__ website](https://bustools.github.io/download). Binary installation time is less than two minutes.

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

capture         Capture records from a BUS file
correct         Error correct a BUS file
count           Generate count matrices from a BUS file
inspect         Produce a report summarizing a BUS file
linker          Remove section of barcodes in BUS files
project         Project a BUS file to gene sets
sort            Sort a BUS file by barcodes and UMIs
text            Convert a binary BUS file to a tab-delimited text file
whitelist       Generate a whitelist from a BUS file

Running bustools <CMD> without arguments prints usage information for <CMD>
~~~

### capture
`bustools capture` can separate BUS files into multiple files according to the capture criteria.

~~~
Usage: bustools capture [options] bus-files

Options: 
-o, --output          Directory for output 
-c, --capture         List of transcripts to capture
-e, --ecmap           File for mapping equivalence classes to transcripts
-t, --txnames         File with names of transcripts
~~~


### correct
BUS files can be barcode error corrected with respect to a technology-specific whitelist of barcodes using `bustools correct`.

~~~
> bustools correct
Usage: bustools correct [options] bus-files

Options: 
-o, --output          File for corrected bus output
-w, --whitelist       File of whitelisted barcodes to correct to
-p, --pipe            Write to standard output
~~~

### count
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

### inspect
A report summarizing the contents of a sorted BUS file can be output either to standard out or to a JSON file for further analysis using `bustools inspect`.

~~~
> bustools inspect
Usage: bustools inspect [options] sorted-bus-file

Options: 
-o, --output          File for JSON output (optional)
-e, --ecmap           File for mapping equivalence classes to transcripts
-w, --whitelist       File of whitelisted barcodes to correct to
-p, --pipe            Write to standard output
~~~

`--ecmap` and `--whitelist` are optional parameters; `bustools inspect` is much faster without them, especially without the former.

Sample output (to stdout):
~~~
Read in 3148815 BUS records
Total number of reads: 3431849

Number of distinct barcodes: 162360
Median number of reads per barcode: 1.000000
Mean number of reads per barcode: 21.137281

Number of distinct UMIs: 966593
Number of distinct barcode-UMI pairs: 3062719
Median number of UMIs per barcode: 1.000000
Mean number of UMIs per barcode: 18.863753

Estimated number of new records at 2x sequencing depth: 2719327

Number of distinct targets detected: 70492
Median number of targets per set: 2.000000
Mean number of targets per set: 3.091267

Number of reads with singleton target: 1233940

Estimated number of new targets at 2x seuqencing depth: 6168

Number of barcodes in agreement with whitelist: 92889 (57.211752%)
Number of reads with barcode in agreement with whitelist: 3281671 (95.623992%)
~~~

### linker
`bustools linker` removes specified section of barcode in BUS files.

~~~
Usage: bustools linker [options] bus-files

Options: 
-s, --start           Start coordinate for section of barcode to remove (0-indexed, inclusive)
-e, --end             End coordinate for section of barcode to remove (0-indexed, exclusive)
-p, --pipe            Write to standard output
~~~

If `--start` is -1, the removed section begins at beginning of barcode. Likewise, if `--end` is -1, the removed section ends at the end of the barcode. BUS files should contain barcodes of the same length.

### project
The `kallisto bus` command maps reads to a set of transcripts. `bustools project` takes as input kallisto's (sorted) output and a transcript to gene map (tr2g file), and outputs a BUS file, a matrix.ec file, and a list of genes, which collectively map each read to a set of genes.

~~~
Usage: bustools project [options] sorted-bus-file

Options: 
-o, --output          File for project bug output and list of genes (no extension)
-g, --genemap         File for mapping transcripts to genes
-e, --ecmap           File for mapping equivalence classes to transcripts
-t, --txnames         File with names of transcripts
-p, --pipe            Write to standard output
~~~

### sort

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

### text

BUS files can be converted to a tab-separated format for easy inspection and processing using shell scripts or high level languages with `bustools text`.

~~~
> bustools text
Usage: bustools text [options] bus-files

Options: 
-o, --output          File for text output
~~~

### whitelist
`bustools whitelist` generates a whitelist based on the barcodes in a sorted BUS file.

~~~
Usage: bustools whitelist [options] sorted-bus-file

Options: 
-o, --output        File for the whitelist
-f, --threshold     Minimum number of times a barcode must appear to be included in whitelist
~~~

`--threshold` is a (highly) optional parameter. If not provided, `bustools whitelist` will determine a threshold based on the first 200 to 100,200 records.
