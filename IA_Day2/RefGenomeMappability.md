To generate a reference genome mappability file yourself, the following two commands are needed:

```
cd ref;
$HMMCOPY_DIR/util/mappability/internal/fastaToRead -w 35 Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa | bowtie2 -x $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa -f /dev/stdin -p 16 -N 0 -k 5 --quiet | grep "20:" | cut -f 1 | uniq -c | awk '{print $2"\tign\tign\tign\tign\tign\t"$1}' | $HMMCOPY_DIR/util/mappability/internal/readToMap.pl -m 4 | $HMMCOPY_DIR/util/bigwig/wigToBigWig stdin Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.sizes Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.bw

$HMMCOPY_DIR/bin/mapCounter \
    -w 1000 \
    $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.bw \
    > $IA_HOME/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr20_adjusted.fa.map.ws_1000.wig
```