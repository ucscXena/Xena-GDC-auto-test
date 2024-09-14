# Quality Assurance of Data Imported from the Genomic Data Commons (GDC) into the Xena Browser

The Xena Browser is a web-based tool that enables researchers to visualize
and analyze cancer functional genomics data by integrating genomic,
clinical, and survival data from different sources, including the Genomic
Data Commons (GDC). To ensure the integrity of the data imported
from the GDC via the Xena GDC Extract Transform Load (ETL)
pipeline, automated quality assurance testing code has been developed
here to validate this data. This testing code utilizes orthogonal data
retrieval and comparison methods, independently verifying genomic, clinical, 
and survival data against data retrieved from the GDC. Multiple
methods are implemented to accomplish this, such as rounding optimizations
for gene expression data, handling multiple data values for
single samples, addressing inconsistencies in the clinical and mutation
files, and creating a logger file to document all inconsistencies. The
testing code has been rigorously tested using flawed data and files to
ensure its effectiveness in identifying errors. This automated approach
not only saves time in QAing this data but also ensures the quality
of important cancer functional genomics data for research purposes
