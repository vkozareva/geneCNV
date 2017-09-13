task CreateMatrix {
  String sample_name
  File input_bam
  File input_bam_index

  File baseline_intervals

  Int preemptible_tries
  Int disk_size

  command {
    cnv \
      create-matrix \
        --outputFile ${sample_name}.csv \
        ${baseline_intervals} \
        ${input_bam}
  }

  runtime {
    docker: "us.gcr.io/genepeeks-bioinformatics/genecnv:head"
    memory: "500 MB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }

  output {
    File output_matrix = "${sample_name}.csv"
  }
}

workflow CreateMatrixWf {
  File samples_file
  File baseline_intervals

  Int preemptible_tries
  Int disk_size

  Array[Object] samples = read_objects(samples_file)

  scatter (sample in samples) {
    call CreateMatrix {
      input:
        sample_name = sample.sample,
        input_bam = sample.bam,
        input_bam_index = sample.bam_index,
        baseline_intervals = baseline_intervals,
        preemptible_tries = preemptible_tries,
        disk_size = disk_size
    }
  }

  output {
    Array[File] output_matrices = CreateMatrix.output_matrix
  }
}

