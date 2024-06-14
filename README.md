## <b>*Alu*-mediated DNA Transduction Tracker and Verifier</b>

---

### <b>Overview</b>
Transduction Tracker (TransductionTracker.V2.py) and Verifier (TransductionVerifier.V2.py) are tools designed to track and verify DNA transductions caused by non-LTR retroelements, optimized for "*Alu*s" and especially "*Alu*Ys."  Therefore, the current version is not guaranteed for other non-LTR retroelements such as L1s and SVAs. The tools are implemented in Python 3.10.12 and require specific Python packages and external tools to function correctly.

---
### <b>Prerequisites</b>
Before running the tools, ensure you have the following installed:

#### <span style="border-bottom: 1.5px solid; text-decoration: none;">Python packages</span>
- biopython version: 1.81 (or later)
- colorama version: 0.4.6 (or later)
- joblib version: 1.3.2 (or later)
- regex
- pandas 

You can install these packages using pip3 as follows:

```bash
pip3 install biopython==1.81 colorama==0.4.6 joblib==1.3.2 regex pandas
```

#### <u>External tools</u>
- bedtools version: 2.31.0
- pblat version: 2.5.1
- [pslScore.pl](https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl) script from UCSC
  
---

### <b>Installation</b>
1. Download the code: go to the GitHub repository and click the "Code" button, then select "Download ZIP".  Extract the downloaded ZIP file to your desired location.

2. Navigate to the project directory:
```bash
cd path/to/transduction-tracker-verifier
```
3.  Install the required Python packages:
```bash
pip3 install -r requirements.txt
```

> <b>Tip:</b> Create a requirements.txt file with the following content:

```
biopython==1.81
colorama==0.4.6
joblib==1.3.2
regex
pandas
```

4. Ensure `bedtools (version 2.31.0 or later)` and `pblat (version 2.5.1 or later)` are installed and available in your PATH. Additionally, [pslScore](https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl) script needs to be downloaded from UCSC.

---

### <b>Getting started</b>

The software package comprises two main programs that need to be executed in the following order:
1. TransductionTracker.V2.py
2. TransductionVerifier.V2.py

#### <span style="border-bottom: 1.5px solid; text-decoration: none;">Running transduction tracker </span>
```bash
python TransductionTracker.V2.py
usage: python TransductionTracker.V2.py -r <repeat_annotation_path> -g <reference_genome> -s <segmental_duplications (optional)> -o <output_directory> -u <upstream_flank> -d <downstream_flank> -a <Alu_type> -c <cpu_number>
```
Parameters:

```bash
  -h, --help            show this help message and exit

Inputs Information:
  -r , --repeat_annotation 
                        Path to the repeat annotation file.
  -g , --reference_genome 
                        path to the reference genome in fasta format.
  -s , --segmental_duplication 
                        Path to the Segmental Duplication (bed format) file (optional).

Output Directory:
  -o , --output_directory 
                        Path to the directory where output files will be saved.

Analysis Parameters:
  -u , --upstream_base_count 
                        Number of upstream bases to include in the search for 5' TSDs. The default is 100.
  -d , --downstream_base_count 
                        Number of downstream bases to include in the search for 3' TSDs. The default is 4500.
  -a , --alu_type       Type of Alu element to search for. The default is 'AluY'.

Performance:
  -c , --cpus           Number of CPUs to use for processing. The default is 4.
```
Example: 

```bash
python TransductionTracker.V2.py -r repeats.bed -g genome.fasta -s seg_dup.bed -o output_directory -u 100 -d 4500 -a AluY -c 4
```

#### <span style="border-bottom: 1.5px solid; text-decoration: none;">Running transduction verifier</span>
```bash 
python TransductionVerifier.V2.py
usage: python TransductionVerifier.V2.py -i <potential_transd.tsv> -r <RepeatMasker.out> -g <reference_genome.fa> -t <Alu family> -s <path/to/pslScore.pl> -c <cpu_number> 
```
Parameters:

```bash
  -h, --help            show this help message and exit
  -i , --transduction_file 
                        Path to the potential transduction file (TSV format) detected in by TransductionTracker.py.
  -r , --repeat_annotation 
                        Path to the segmental duplication-free RepeatMasker output file (.out format). Use "*_NoSegDup_RM.out" if SegDup annotation was provided earlier with "TransductionTracker.V2.py".
  -g , --reference_genome 
                        Path to the reference genome file.
  -t , --alu_type       Type of Alu element to process. The default is 'AluY'.
  -s , --pslScore       Path to the pslScore.pl script downloaded from UCSC.
  -c , --cpus           Number of CPUs for the alignment. The default is 4.
```
Example:

```bash
python TransductionVerifier.V2.py -i alu_transd_with_pA.tsv -r RM.out -g ref_genome.fa -t AluY -s pslScore.pl -c 4
```

---
### <b>Outputs</b>

1. `TransductionTracker.V2.py` generates several files, including:
    * *_heteroTSD.tsv: contains detected TSDs and related information, excluding those consisting solely of T, A, C, or G homopolymers. This file includes TSDs for *Alu* elements that do not overlap with segmental duplications, if such annotation has been provided to "TransductionTracker.V3.py".
    * *_CHM13v2_potential_alu_transd_with_pA.tsv: contains unfiltered identified transductions. This file serves as input for "TransductionVerifier.V3.py".
      
The script also produces other files with self-explanatory names that can be used for further investigation.

2. `TransductionVerifier.V2.py` produces multiple files, including:
    * Results of alignment and post-processing, which are placed in the alignment folder.
    * *_verified_transductions.tsv: comprises transductions for which the potential source was identified and met our other verification criteria. A manual inspection of the results is recommended.
  







