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

