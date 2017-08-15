#!/bin/bash

# Download all required data to initialize the pipeline
download_folder='data/'

wget -i data/data_urls.txt '--directory-prefix='$download_folder
tar -xf data/gbm_tcga.tar.gz -d $download_folder
