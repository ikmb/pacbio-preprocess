#!/usr/bin/env nextflow

/*

||  ||
Pacbio Preprocessing
||  ||

*/

def helpMessage() {
  log.info"""
======================================
>> Processing of PacBio reads <<
======================================

Usage:

A typical command to run this pipeline would be:

nextflow run ikmb/pacbio-preprocessing --bam movie.bam --qc --demux --ccs

Mandatory arguments:

--bam			A PacBio movie file in BAM format

""".stripIndent()
}

// Show help message
if (params.help){
        helpMessage()
        exit 0
}

// Barcode file
barcodes = file("$baseDir/assets/Sequel_16_Barcodes_v3.fasta")

// Enables splitting of CCS read generation into 10 parallel processes
def chunks = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]

if (!params.bam) {
	exit 1, "Must provide a BAM movie file (--bam)"
} else {
	Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "Could not find an input bam file" }
		.map { b -> [ file(b).getBaseName(), file(b), file("${b}.pbi") ] }
		.into { bamFile; bamQC }
}

// requested generation of HiFi reads
if (params.hifi) {

	process BamToCCS {
		
		scratch true

		publishDir "${params.outdir}/${sample}/CCS", mode: 'copy',
			saveAs: { filename ->
				if( filename.indexOf("report.txt") > 0 ) filename
				else null
			}

		input:
		set val(sample),file(bam),file(bam_index) from bamFile
		each chunk from chunks


		output:
		set val(sample),file(reads) into ReadChunks
		file(report)

		script:
		reads = bam.getBaseName() + "." + chunk + ".ccs.bam"
		report = bam.getBaseName() + "." + chunk + ".ccs.report.txt"
		options =""
		if (params.max_length) {
			options = "--max-length=${params.max_length}"
		}

		"""
			ccs $bam $reads --min-passes ${params.min_passes} --min-length=${params.min_length} --chunk $chunk/10 -j ${task.cpus}
			mv *ccs_report.txt $report
		"""

	}

	ReadChunksGrouped = ReadChunks.groupTuple()

	process CcsMerge {

		scratch true 

	        publishDir "${params.outdir}/${sample}/CCS", mode: 'copy'

		input:
		set val(sample),file(read_chunks) from ReadChunksGrouped

		output:
		set val(sample),file(bam),file(pbi) into mergedReads

		script:
		bam = sample + "." + params.min_passes + "_min_passes" + "." + params.min_length + "_min_length.reads.ccs.bam"
		pbi = bam + ".pbi"

		

		"""
			pbmerge -o $bam $read_chunks
			pbindex $bam
		"""
	}

	process GetHiFiReads {

                publishDir "${params.outdir}/${sample}/HiFi", mode: 'copy'

		input:
		set val(sample),file(bam),file(pbi) from mergedReads

		output:
		set val(sample),file(hifi),file(hifi_pbi) into HiFiBam

		script:
		hifi = bam.getBaseName() + ".hifi.bam"
		hifi_pbi = hifi + ".pbi"

		"""
			extracthifi $bam $hifi
			pbindex $hifi
		"""

	}

} else {
	HiFiBam = bamFile
}

if (params.demux) {

	process PbDemux {

		publishDir "${params.outdir}/${sample}/Demux", mode: 'copy'

		input:
		set val(sample),file(bam),file(pbi) from HiFiBam

		output:
		set val(sample),file("*-*.bam") into final_bams
		set val(sample),file("*.lima.*") into lima_reports

		script:
		def options = ""
		if (params.hifi) {
			options = "--ccs --min-score 80"
		}
		demux = bam.getBaseName() + ".demux.bam"

		"""
			lima $bam $barcodes $demux --same $options --split
		"""

	}

} else {

	final_bams = HiFiBam
}

if (params.qc) {

	process bam2fasta {

		input:
		set val(sample),file(bam),file(pbi) from bamQC

		output:
		set val(sample),file(fasta) into fasta_reads

		script:
		fasta = bam.getBaseName() + ".fasta.gz"

		"""
			bam2fasta -o ${bam.getBaseName()} $bam
		"""
		
	}

	process nanoplot {

        	publishDir "${params.outdir}/${sample}/QC", mode: 'copy'

	        input:
        	set val(sample),file(fasta) from fasta_reads

	        output:
        	file("nanoplot")

	        script:

        	"""
                	NanoPlot -t ${task.cpus} -o nanoplot --fasta $fasta
	        """
	}
}
