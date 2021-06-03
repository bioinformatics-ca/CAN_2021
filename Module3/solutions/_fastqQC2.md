Here is the command to generate the new quality graphs:

```{.bash}
# Generate post-trimmed QC
mkdir -p postTrimQC/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 33 \
  --read1 reads/normal/runBD06UFACXX_4/normal.t30l50.pair1.fastq.gz \
  --read2 reads/normal/runBD06UFACXX_4/normal.t30l50.pair2.fastq.gz \
  --threads 2 --regionName normalrunBD06UFACXX_4 --output postTrimQC/
```
