#!/bin/bash

INDEX_TYPES="block_interpolative block_simple16 block_anspacked block_ansmsb block_ansmsbminmax block_ansmsbmed90p block_ansmsbmed90pmerged"

# GOV2
GOV2_PATH="/scratch/mpetri/gov2-d2si-new/gov2"
GOV2_OUTPATH="/scratch/mpetri/experiments_wsdm/gov2-"
for INDEX in $INDEX_TYPES
do
    echo GOV2 $INDEX
    echo ./build/create_freq_index $INDEX $GOV2_PATH $GOV2_OUTPATH$INDEX  --check &> output_GOV2.txt
done

# NEWS
NEWS_PATH="/scratch/mpetri/CC-NEWS-NEW/CC-NEWS-20160901-20170228-en"
NEWS_OUTPATH="/scratch/mpetri/experiments_wsdm/news-"
for INDEX in $INDEX_TYPES
do
    echo NEWS $INDEX
    echo ./build/create_freq_index $INDEX $NEWS_PATH $NEWS_OUTPATH$INDEX --check &> output_NEWS.txt
done
