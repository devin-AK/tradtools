# First, build hg19 genome database following https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/build_genome_database.md
# e.g. bash scripts/download_genome_data.sh hg19 [DESTINATION_DIR]
# Then run the pipeline
cd /work/devin/trad/lamin
caper run /work/devin/trad/chip-seq-pipeline2/chip.wdl -i /work/devin/trad/lamin/LMNA_GSE54334.json --singularity --max-concurrent-tasks 8
caper run /work/devin/trad/chip-seq-pipeline2/chip.wdl -i /work/devin/trad/lamin/LMNA_HSF_CTL-1_GSE81671.json --singularity --max-concurrent-tasks 8
caper run /work/devin/trad/chip-seq-pipeline2/chip.wdl -i /work/devin/trad/lamin/LMNA_HSF_CTL-2_GSE81671.json --singularity --max-concurrent-tasks 8
caper run /work/devin/trad/chip-seq-pipeline2/chip.wdl -i /work/devin/trad/lamin/LMNA_HSF_CTL-3_GSE81671.json --singularity --max-concurrent-tasks 8
caper run /work/devin/trad/chip-seq-pipeline2/chip.wdl -i /work/devin/trad/lamin/LMNB1_GSE49341.json --singularity --max-concurrent-tasks 8
# Next, organize outputs using Croo
# Lastly, use signal tracks (.bigwig format) as indicated in chromatin.R script