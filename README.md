RNA-seq Isoform Analysis Pipeline
=================================

This pipeline identifies and quantifies transcript isoforms for user-specified genes (default: MS4A4A, MS4A6A) from RNA-seq data. It is designed for single-condition studies (e.g., microglia) and is easily configurable for new genes or datasets.


Quick Start
-----------
1. **Set up the environment:**
    - Install [conda](https://docs.conda.io/en/latest/miniconda.html) if not already installed.
    - Create and activate the environment:
      ```sh
      conda env create -f environment.yml
      conda activate splice-isoform-pipeline
      ```

2. **Prepare your data:**
    - Place raw FASTQ files in `data/raw/`.
    - Download reference files (genome.fa, transcripts.fa, reference.gtf) to `data/reference/`.
    - Edit `config/samples.csv` to list your samples and FASTQ paths.
    - Edit `config/config.yaml` to set the genes you want to analyze (see below).

2. **Prepare your data:**
    - Place raw FASTQ files in `data/raw/`.
    - Download reference files (genome.fa, transcripts.fa, reference.gtf) to `data/reference/`.
    - Edit `config/samples.csv` to list your samples and FASTQ paths.
    - Edit `config/config.yaml` to set the genes you want to analyze (see below).

3. **Set genes of interest:**
    - Open `config/config.yaml` and edit the `genes_of_interest` list under `r_analysis:`. Example:
      ```yaml
      r_analysis:
         genes_of_interest:
            - MS4A4A
            - MS4A6A
      ```

4. **Build indices (one-time):**
    - Salmon: `salmon index -t data/reference/transcripts.fa -i data/reference/salmon_index -k 31`
    - HISAT2 (optional): `hisat2-build -p 8 data/reference/genome.fa data/reference/hisat2_index/GRCh38`

5. **Run the pipeline:**
    - From the project root, run: `snakemake --cores 8`
    - This will perform QC, trimming, Salmon quantification, and (optionally) alignment/assembly.

6. **Analyze isoforms:**
    - Use the Python module in `analysis/isoform_analysis.py` to extract and plot isoform proportions for your selected genes.
    - Example usage is provided in the module docstring.

7. **Add new samples or genes:**
    - Add new FASTQs to `data/raw/` and update `config/samples.csv`.
    - Edit `genes_of_interest` in `config/config.yaml` as needed.
    - Re-run the pipeline.

File Structure
--------------
See comments in each folder and config file for details. Main files to edit:
- `config/samples.csv`: List of samples and FASTQ files
- `config/config.yaml`: Reference paths and genes of interest
- `analysis/isoform_analysis.py`: Python code for isoform quantification/plots

Outputs
-------
- `salmon_quant/`: Salmon quantification results (quant.sf)
- `alignments/`: BAM files (if running alignment)
- `analysis/`: Downstream results, merged GTFs, and logs

For more details on each step, see the comments in the Snakefile and scripts/.
HISAT2 alignment:
hisat2 -p 8 --dta -x hisat2_index/GRCh38 \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -S alignments/sampleX.sam
samtools view -bS sampleX.sam | samtools sort -@4 -o alignments/sampleX.sorted.bam
samtools index alignments/sampleX.sorted.bam
StringTie assembly:
stringtie alignments/sampleX.sorted.bam -p 8 \
    -G reference.gtf -o stringtie/sampleX.stringtie.gtf -A stringtie/sampleX.gene_abund.tab -e -B
Merge all sample assemblies:
ls stringtie/*.gtf > stringtie/mergelist.txt
stringtie --merge -p 8 -G reference.gtf -o stringtie/merged.stringtie.gtf stringtie/mergelist.txt
gffcompare -r reference.gtf -o stringtie/gffcompare merged.stringtie.gtf
Step 4: Analysis
Python (primary quantification & dominance)
Read Salmon quant.sf files.
Compute transcript proportions per gene.
Identify dominant isoforms per sample and per condition.
Generate plots (stacked bar charts or proportion plots for MS4A4A and MS4A6A).
R (optional / DTU analysis)
Use tximport to import Salmon counts.
Run DRIMSeq or IsoformSwitchAnalyzeR to detect differential transcript usage.
Step 5: Validation & Visualization
Inspect BAMs in IGV for junction support.
Optional sashimi plots using IGV or ggsashimi.
Experimental validation of novel isoforms if detected.
Step 6: Adding New Samples
Place new FASTQs in data/raw/.
Update samples.csv.
Run Salmon quantification and append to existing analysis.
If novel isoform discovery is required, re-run HISAT2 + StringTie on the new samples or the full cohort.
Notes
Default genome: GRCh38 (human).
Default transcriptome: GENCODE or RefSeq.
For low-RAM machines (16 GB), use Salmon for quantification; HISAT2 + StringTie is feasible for novel isoforms. STAR requires ~32 GB RAM.
Maintain reproducibility using scripts in scripts/ and keep config/samples.csv updated.
