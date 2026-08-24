/*
 * GET_OUTGROUP{} will download 1 genome for specified tax id
 */

process GET_OUTGROUP {
	conda "${projectDir}/envs/get_data.yaml"

	input:
	val outgroup_ID
	
	output:
	path '*fna', emit: outgroupFasta

	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${outgroup_ID}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter] OR scaffold level[filter] OR contig level[filter]))" \\
		| esummary \\
		| xtract -pattern DocumentSummary -element RefSeq_category AssemblyStatus FtpPath_RefSeq \\
		| awk -F'\\t' '{
			category=\$1
			status=\$2
			cat_priority=3
			if (category == "reference genome") cat_priority=1
			else if (category == "representative genome") cat_priority=2

			stat_priority=5
			if (status == "Complete Genome") stat_priority=1
			else if (status == "Chromosome") stat_priority=2
			else if (status == "Scaffold") stat_priority=3
			else if (status == "Contig") stat_priority=4

			print cat_priority "\\t" stat_priority "\\t" \$0
		}' \\
		| sort -t \$'\\t' -k1,1n -k2,2n \\
		| cut -f5 \\
		| while read url; do
        
		if [ "\$counter" -ge 1 ]; then
			break
		fi
	
		if [ -z "\$url" ]; then
			continue
		fi
        
		fname="\$(basename "\$url")"
		downloadMeAndCheck() {
			wget "\$url/\${fname}_genomic.gbff.gz"
	        	wget "\$url/\${fname}_genomic.fna.gz"

			gunzip -f "\${fname}_genomic.gbff.gz"
			gunzip -f "\${fname}_genomic.fna.gz"


			if [ -f "\${fname}_genomic.gbff.gz" ] || [ -f "\${fname}_genomic.fna.gz" ]; then	
				rm -f "\${fname}_genomic.gbff.gz" "\${fname}_genomic.fna.gz"
				return 1
			fi

			return 0
		}

		while ! downloadMeAndCheck; do
			echo "Files were corrupted. Retrying"
			sleep 3
		done
		# Adding to the counter
	        counter="\$((counter + 1))"     
                fna_count=\$(ls -1 *.fna 2>/dev/null | wc -l)
                gbff_count=\$(ls -1 *.gbff 2>/dev/null | wc -l)

                # check that there is already one file
                if [ "\$fna_count" -lt 1 ] || [ "\$gbff_count" -lt 1 ]; then
                        echo "It seems a file was already downloaded. Moving on"
                        break
                else
                        continue
                fi

	done
	"""
}
