nextflow.enable.dsl=2

/* parameters
*/
params.fasta = "$projectDir/dm6.fa"
params.gtf = "$projectDir/Drosophila_melanogaster.BDGP6.28.102.chr.gtf"
params.sample = "$projectDir/reads/*.fastq.gz"
params.readLength = 74
params.rmsk = "$projectDir/dm6_rmsk_short.bed"

//Fasta_ch = channel.fromPath(params.Fasta, checkIfExists:true)
    // this does not generate a tuple containing the .fa prefix, which might be useful later... 
fasta_ch = channel.fromFilePairs(params.fasta, size:1) // returns [dm6, [home/jk/nf-dm/dm6.fa]]
gtf_ch = channel.fromPath(params.gtf) // I don't need the prefix from this file
//sample_ch = channel.fromPath(params.sample)
//sample_ch = channel.fromPath(params.sample,checkIfExists : true ).map{
    //tuple(it.name.split('\\.')[0],it)
    //}// this returns a tuple with sample name and paths --- looks complicateds
sample_ch = channel.fromFilePairs(params.sample, size: 1) // probably the easiest way to make a tuple even though there are no paired reads? 
rmsk_ch = channel.fromFilePairs(params.rmsk, size:1)

//starFasta_ch = channel.fromPath("$baseDir/index/genomes/mm10.fa")
//gtf_ch = channel.fromPath("$baseDir/index/genes/gencode.vM22.annotation.gtf")
///reads_ch = channel.of(params.reads)

process STAR_INDEX {
    publishDir "$projectDir/index/STAR"
    cpus 8

    input:
        tuple val(build), path(fasta)
            // instead of path fasta from channel.path()
        path gtf
        val read_length
    output: // The script will generate a folder called STARindex and files within. Nextflow will check if such directory was indeed created
        path "$build", emit: index
        path "version.yml"
    script:
    """
    
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir $build --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $read_length --limitGenomeGenerateRAM 12000000000 --genomeSAindexNbases 12
    cat << EOF >> version.yml
    $task.process:
        STAR: \$(STAR --version 2>&1)
    EOF
    """
    // STAR expects the output sub-directory to be there
}

process STAR_ALIGN {
    publishDir "$projectDir/aligned", mode: "copy"
    cpus 1


    input:
        path index //  Exit code 137 caused by lack of memory if there are too many samples...
        tuple val(sample_id), path(sample_path)
        
    
    output: // All of the outputs from STAR must be listed here or else they won't be kept
            // Nextflow checks to see if these output files were indeed generated 
            // A minor problem is that STAR takes in output prefix, and automatically appends the rest 
            // This means that the code must be run separately to know what outputs (and their names) will be generated. 
            
        tuple val(sample_id), path("*Aligned.out.sam"), emit : bam
        //tuple val(sample_id), path("*Log.final.out")
        //tuple val(sample_id), path("*Log.out")
        //tuple val(sample_id), path("*Log.progress.out")
        //tuple val(sample_id), path("*SJ.out.tab")

        // not really important to check if other files were generated. Makes dag look simpler.





    script:
  //  def prefix = task.ext.prefix ?: "${sample_id}[0]"
    """
    STAR --runMode alignReads --runThreadN $task.cpus --genomeDir $index --readFilesCommand zcat --readFilesIn ${sample_path[0]} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix $sample_id
    """

}


process SAMTOOLS_VIEW_SORT {
    publishDir "$projectDir/aligned"
    cpus 2

    input:
        tuple val(sample_id), path(sample_path)
    output: 
        tuple val(sample_id), path("*.sorted.bam"), emit : sortedBam

    script:
    """
    samtools view -b -@ $task.cpus $sample_path | samtools sort -@ $task.cpus - -o ${sample_id}.sorted.bam 
    """
    // for this command, -o must indicate the exact output file name, matching output name in output:
}

process HTSEQ_COUNT {
    publishDir "$projectDir/countTable", mode: "copy"
    cpus 2

    input:
        tuple val(sample_id), path(sample_path)
        path gtf
    output:
        tuple val(sample_id), path("*.csv")
    script:
    """
    htseq-count -f bam -s no -q $sample_path $gtf > ${sample_id}.csv
    """

}

process BT2_INDEX {
    publishDir "$projectDir/index/bt2"
    memory 7.GB


    input:
        tuple val(build), path(fasta)
    output:
    /*
        tuple val(build), path("*.1.bt2") //output is a prefix
        tuple val(build), path("*.2.bt2")
        tuple val(build), path("*.3.bt2")
        tuple val(build), path("*.4.bt2")
        tuple val(build), path("*.rev.1.bt2")
        tuple val(build), path("*.rev.2.bt2")
    */
        tuple val("$build"), path("$build*")
        // BT2ALIGN expects index files' prefix, but bowtie2-build outputs various .bt2 files; $build* is required to check for .bt2 files

    script:
    """
    bowtie2-build  $fasta $build
    """
}


process BT2_ALIGN {
    publishDir "$projectDir/aligned"
    cpus 4

    input:
        tuple val(sample_id), path(fasta)
        tuple val(prefix), path(bt2) //path(bt2) is never actually used
    output:
        tuple val(sample_id), path("*.sam")
    script:
    """
    bowtie2 -q --very-fast-local -p $task.cpus -x $prefix -U $fasta -S ${sample_id}.sam
    """

}

process SAMTOOLS_VIEW_BT2 {
    publishDir "$projectDir/aligned"
    cpus 4

    input:
        tuple val(sample_id), path(sam)
    output:
        tuple val(sample_id), path("*.bam")
    script:
    """
    samtools view -b -@ $task.cpus $sam -o ${sample_id}.bam
    """
}

process REPENRICH_SETUP {
    publishDir "$projectDir/RepEnrich2_setup", mode: "copy"
    conda  "$projectDir/py2.yml"
    memory 3.GB
    cpus 2

    input:
        tuple val(build_ignore), path(rmsk)
        tuple val(build), path(fasta) // remember the channel is a tuple
    output:
        path "RepEnrich2Index"
    script:
    """
    RepEnrich2_setup.py  $rmsk $fasta RepEnrich2Index --is_bed TRUE 
    """
}

process REPENRICH_SUBSET {
    publishDir "$projectDir/aligned", mode: "copy"
    conda  "$projectDir/py2.yml"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*_unique.bam"), path("*_multimap.fastq")
    script:
    """
    RepEnrich2_subset.py $bam 30 $sample_id
    """
}

process REPENRICH_SORT_INDEX {
    publishDir "$projectDir/aligned", mode: "copy"
    conda  "$projectDir/py2.yml"

    input:
        tuple val(sample_id), path(unique), path(multimap)
    output:
        tuple val(sample_id), path("*.sorted_unique.bam"), path("*.sorted_unique.bam.bai")
    script:
    """
    samtools sort -@ $task.cpus $unique -o ${sample_id}.sorted_unique.bam
    samtools index -b ${sample_id}.sorted_unique.bam
    """
}

process REPENRICH {
    publishDir "$projectDir/RepEnrich2", mode:"copy"
    conda  "$projectDir/py2.yml"
    
    input:
        tuple val(build_ignore), path(rmsk)
        path setup  
        tuple val(sample_id), path(sorted_unique), path(samtools_index)
        tuple val(sample_id), path(unique), path(multimap)
    output:
        tuple val(sample_id),  path("${sample_id}_fraction_counts.txt") // why is csv the output?
    script:
    """
    RepEnrich2.py $rmsk . $sample_id $setup $multimap $sorted_unique --is_bed TRUE
    """
    // simply output to . (current directory i.e. work) which will be copied to publishDir anyway. This process can't find the output files if it's another directory
}

workflow {
    //sample_ch.view()
    STAR_INDEX(fasta_ch, gtf_ch, params.readLength)
    STAR_ALIGN(STAR_INDEX.out.index.collect(),sample_ch)
        // .collect() is the KEY. Or else it returns a queue channel which gets used up after one task 
    SAMTOOLS_VIEW_SORT(STAR_ALIGN.out.bam)
        // assigning each workflow item to another channel doesn't work i.e.
        // ch_STARgenomeGenerate = STARgenomeGenerate(...)
        // STARalign(ch_STARgenomeGenerate.out.collect(), sample_ch) fails
        // I guess it's too redundant anyway?
    HTSEQ_COUNT(SAMTOOLS_VIEW_SORT.out.sortedBam, gtf_ch.collect()) //.collect seems to be the magic operator to iterate over multiple files when there are many vs. one inputs
    BT2_INDEX(fasta_ch)
    BT2_ALIGN(sample_ch, BT2_INDEX.out.collect())
    SAMTOOLS_VIEW_BT2(BT2_ALIGN.out)
    REPENRICH_SETUP(rmsk_ch, fasta_ch)
    REPENRICH_SUBSET(SAMTOOLS_VIEW_BT2.out)
    REPENRICH_SORT_INDEX(REPENRICH_SUBSET.out)
    //REPENRICHSETUP(rmsk_ch,fasta_ch)
    REPENRICH(rmsk_ch.collect(), REPENRICH_SETUP.out.collect(), REPENRICH_SORT_INDEX.out, REPENRICH_SUBSET.out) // anythingl
    


                                                               

}
