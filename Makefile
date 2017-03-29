# Generate emacs TAGS files that indexes relevant source files.

SOURCE_FILE_MATCH = \
	"(" -type d -a "(" -name "bams" -o -name ".ipynb_checkpoints" -o -name "build" -o -name "dist" \
			-o -name "obj" -o -name "bin" -o -name ".git" \
			")" -a \! -prune ")" \
	-o \! "(" -type d -o -type l \
		-o -name "*.odt" \
		-o -name "*.pyc" -o -name "*.a" -o -name "*.o" \
		-o -name "*.jar" \
		-o -name "*.rrd" -o -name "*.out" -o -name "*.log" -o -name "*.txt" -o -name "*.csv" \
		-o -name "*.png" -o -name "*.gif" -o -name "*.jpg" -o -name "*\#" -o -name "*~" \
		-o -name "*.gz" -o -name "*.tgz" -o -name "*.tsv" \
		-o -name "TAGS" -o -name "TESTTAGS" \
		-o -name "optimize_vgd_new.py" -o -name "optimize_gene_liability.py" \
	")"

.PHONY: tags
tags:
	rm -f TAGS
	find . $(SOURCE_FILE_MATCH) -a -exec etags -a {} \;

	# rm -f TESTTAGS
	# find . $(SOURCE_FILE_MATCH) -a -path "*/tests/*" -a -exec etags --output=TESTTAGS -a {} \;
