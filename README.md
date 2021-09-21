# pacbio-preprocess

Preprocessing of Pacbio subread data for delivery to customers

## Run options

### `--bam` 

Path to a subread BAM file, or a list of subread BAM files given as `--bam /path/to/*.subreads.bam`

### `--qc` [ false | true (default) ]

Whether to generate QC plots using Nanoplot. 

### `--hifi` [ true | false (default) ]

whether to convert subreads into HiFi reads.

### `--demux` [ true | false (default) ]

Whether to perform demuxing on the data. If "--hifi" is specified, demuxing is performed on the HiFi reads only. Otherwise, the raw subread file is used. 

## Expert options

### `--min_passes` [ default = 3 ]

Number of passes required for a read to be considered for polishing. Requires "--hifi"

### `--min_length` [ default = 1500 ]

Minimum length a read must have after polishing. Requires "--hifi"

### `--max_length`[ default = false ]

Maximum length a read can have after polishing. Requires "--hifi". 



