# bustools

__bustools__ is a program for manipulating [__BUS__](https://github.com/BUStools/BUS) files for single cell 
RNA-Seq datasets. It can be used to error correct barcodes, collapse UMIs, produce gene count or transcript compatibility count matrices, and is useful for many other tasks. See the [__kallisto &#124; bustools website__](https://www.kallistobus.tools/) for examples and instructions on how to use __bustools__ as part of a single-cell RNA-seq workflow.

If you use __bustools__ please cite

Melsted, Páll, Booeshaghi, A. Sina et al. [Modular, efficient and constant-memory single-cell RNA-seq preprocessing.](https://doi.org/10.1038/s41587-021-00870-2) Nature Biotechnology, 2021.

For some background on the design and motivation for the __BUS__ format and __bustools__ see 

Melsted, Páll, Ntranos, Vasilis and Pachter, Lior [The Barcode, UMI, Set format and BUStools](https://doi.org/10.1093/bioinformatics/btz279), Bioinformatics, 2019.


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

sort            Sort a BUS file by barcodes and UMIs
correct         Error correct a BUS file
count           Generate count matrices from a sorted BUS file
inspect         Produce a report summarizing a sorted BUS file
allowlist       Generate an on-list from a sorted BUS file
capture         Capture records from a BUS file
text            Convert a binary BUS file to a tab-delimited text file
fromtext        Convert tab-delimited text file to a binary BUS file
extract         Extracts reads from input FASTQ files based on BUS file
umicorrect      Use a UMI correction algorithm
compress        Compress a sorted BUS file
decompress      Decompress a compressed BUS file

Running bustools <CMD> without arguments prints usage information for <CMD>
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


### correct
BUS files can be barcode error corrected with respect to a technology-specific "on list" of barcodes using `bustools correct`.

~~~
> bustools correct
Usage: bustools correct [options] bus-files

Options: 
-o, --output          File for corrected bus output
-w, --onlist          File of on-listed barcodes to correct to
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
--cm                  Count multiplicites instead of UMIs
~~~


### inspect
A report summarizing the contents of a sorted BUS file can be output either to standard out or to a JSON file for further analysis using `bustools inspect`.

~~~
> bustools inspect
Usage: bustools inspect [options] sorted-bus-file

Options: 
-o, --output          File for JSON output (optional)
-e, --ecmap           File for mapping equivalence classes to transcripts
-w, --onlist          File of on-listed barcodes to correct to
-p, --pipe            Write to standard output
~~~

`--ecmap` and `--onlist` are optional parameters; `bustools inspect` is much faster without them, especially without the former.

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

Number of barcodes in agreement with on-list: 92889 (57.211752%)
Number of reads with barcode in agreement with on-list: 3281671 (95.623992%)
~~~

### allowlist
`bustools allowlist` generates an on-list based on the barcodes in a sorted BUS file.

~~~
Usage: bustools allowlist [options] sorted-bus-file

Options: 
-o, --output        File for the on-list
-f, --threshold     Minimum number of times a barcode must appear to be included in on-list
~~~

`--threshold` is a (highly) optional parameter. If not provided, `bustools allowlist` will determine a threshold based on the first 200 to 100,200 records.


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



### text

BUS files can be converted to a tab-separated format for easy inspection and processing using shell scripts or high level languages with `bustools text`.

~~~
> bustools text
Usage: bustools text [options] bus-files

Options: 
-o, --output          File for text output
~~~


### fromtext

Converts a plaintext tab-separated representation of a BUS file to a binary BUS file.

~~~
> bustools text
Usage: bustools text [options] bus-files

Options: 
-o, --output          File for text output
~~~


### compress

Compresses a sorted BUS file.

~~~
> bustools compress
Usage: bustools compress [options] sorted-bus-file
Note: BUS file should be sorted by barcode-umi-ec

Options: 
-N, --chunk-size    Number of rows to compress as a single block
-o, --output        File to write compressed output
-p, --pipe          Write to standard output.
~~~

Reference for bustools compression:

Einarsson, P and Melsted, Páll [BUSZ: compressed BUS files](https://doi.org/10.1093/bioinformatics/btad295), Bioinformatics, 2023.


### decompress

Decompresses (inflates) a compressed BUS file.

~~~
> bustools compress
Usage: bustools compress [options] sorted-bus-file
Note: BUS file should be sorted by barcode-umi-ec

Options: 
-o, --output        File to write decompressed output
-p, --pipe          Write to standard output.
~~~


