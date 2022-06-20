### Download pool 1 barcode linkage sequencing data

Download raw sequencing reads from SRA using one of the links below

| Type | Size         | Location | Name                                                                      | Free Egress | Access Type  |
| :--- | :----------- | :------- | :------------------------------------------------------------------------ | :---------- | :----------- |
| run  | 1,491,610 Kb | NCBI     | https://sra-download.ncbi.nlm.nih.gov/traces/sra75/SRR/016724/SRR17125831 | worldwide   | anonymous    |
|      |              | AWS      | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17125831/SRR17125831      | worldwide   | anonymous    |
|      |              | GCP      | gs://sra-pub-run-7/SRR17125831/SRR17125831.1                              | gs.US       | gcp identity |

___

Alternatively, download using the [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) as follows:

`prefetch -v SRR17125831`

Convert this file to fastq using 

`fastq-dump --outdir /opt/fastq/ --split-files /home/[USER]/ncbi/public/sra/SRR17125831.sra` 

___

Rename file `100p10.R1.fastq.gz`

Move this file to `~/git/burke_2022/data/fastq`

