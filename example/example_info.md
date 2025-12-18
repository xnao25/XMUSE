## Input files
- intron_vs_exon_ASC.fa: A fasta file containing 1000 intron sequences that cover the region from upstream 100bp of annotated 3'SS to the 3'SS, and 1000 exon sequences that cover the region from the base after the 3'SS (the beginning of the exons) to 100bp downstream from the 3'SS. This file is made for the ASC analysis.
- intron_vs_exon_AWSC.fa: A fasta file containing 1000 intron sequences that cover the region from 50bp to 5bp upstream of annotated 3'SS and 1000 exon sequences that cover the region from 5bp to 50bp downstream of annoated 3'SS. This file is made for the model fitting step in AWSC analysis.
- intron_vs_exon_AWSC_test.fa: A fasta file containing 300 intron sequences and 300 exon sequences from the same region in `intron_vs_exon_AWSC.fa` file. This file is made for the model testing and new data transformation in AWSC analysis.
- intron_vs_exon_UWSC.fa: A fasta file containing 1000 intron sequences that cover the region from 80bp to 10bp upstream of annotated 3'SS and 1000 exon sequences that cover the region from 10bp to 80bp downstream of annotated 3'SS. This file is made for the model fitting step in UWSC analysis.
- intron_vs_exon_UWSC_test.fa: A fasta file containing 300 intron sequences and 300 exon sequences from the same region in `intron_vs_exon_UWSC.fa` file. This file is made for the model testing and new data transformation in UWSC analysis.

- pickle files from below are removed from this repo due to size limit
## Commands:
### Aligned Sequence Characterization:
#### model fitting:
- `python3 ../xmuse.py -if intron_vs_exon_ASC.fa -it intron_vs_exon_label.txt -m ASC -k 4 -t 4 -sw 4 -o ASC_result/`

### Aligned Whole Sequence Characterization:
#### model fitting:
- `python3 ../xmuse.py -if intron_vs_exon_AWSC.fa -it intron_vs_exon_label.txt -m AWSC -k 4 -t 2 -sg 50 -o AWSC_result/`
#### model testing:
- `python3 ../xmuse.py -if intron_vs_exon_AWSC_test.fa -it intron_vs_exon_label_test.txt -fm AWSC_result/AWSC_model.pickle -o AWSC_result/`

### Unaligned Whole Sequence Characterization:
#### model fitting:
- `python3 ../xmuse.py -if intron_vs_exon_UWSC.fa -it intron_vs_exon_label.txt -m UWSC -k 4 -t 4 -sg 50 -o UWSC_result/`
#### model testing:
- `python3 ../xmuse.py -if intron_vs_exon_UWSC_test.fa -it intron_vs_exon_label_test.txt -fm UWSC_result/UWSC_model.pickle -o UWSC_result/`