from mando import command, main


@command('get-samples-files')
def get_samples_files(alignment_index_file, output_prefix, lines=250, path_prefix=''):
    """Get input samples files for the WDL script.

    Splits and reformats alignment index file from 1000 Genomes FTP:
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.exome.alignment.index

    :param alignment_index_file: Alignment index file from 1000 Genomes.
    :param output_prefix: Output prefix.
    :param lines: How many lines per file.
    :param path_prefix: Path prefix to add if any."""

    with open(alignment_index_file, 'r') as aln_fh:
        line_counter = 0
        file_index = 0
        output_fh = open('{prefix}.{index:02d}.csv'.format(prefix=output_prefix, index=file_index), 'w')
        output_fh.write('\t'.join(['sample', 'bam', 'bam_index']))
        output_fh.write('\n')
        for line in aln_fh:
            if '.mapped.' not in line:
                continue
            bam, _, bam_index, _, _, _ = line.split()
            sample = bam.split('/')[1]
            output_fh.write('\t'.join([sample, path_prefix + bam, path_prefix + bam_index]))
            output_fh.write('\n')
            line_counter += 1
            if line_counter >= lines:
                output_fh.close()
                file_index += 1
                line_counter = 0
                output_fh = open('{prefix}.{index:02d}.csv'.format(prefix=output_prefix, index=file_index), 'w')
                output_fh.write('\t'.join(['sample', 'bam', 'bam_index']))
                output_fh.write('\n')

if __name__ == '__main__':
    main()
