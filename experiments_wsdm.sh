#!/bin/bash

INDEX_TYPES="block_ansmsbmedmaxmerged"
#INDEX_TYPES="block_interpolative block_optpfor block_anspacked block_ansmsbmedmaxmerged opt"


# GOV2
GOV2_PATH="/ssd/gov2-d2si-new/gov2"
GOV2_OUTPATH="/ssd/experiments_wsdm/gov2-"
for INDEX in $INDEX_TYPES
do
    echo GOV2 $INDEX
    ./build/create_freq_index $INDEX $GOV2_PATH $GOV2_OUTPATH$INDEX  --check &>> output_GOV2.txt
done

# NEWS
NEWS_PATH="/ssd/CC-NEWS-NEW/CC-NEWS-20160901-20170228-en"
NEWS_OUTPATH="/ssd/experiments_wsdm/news-"
for INDEX in $INDEX_TYPES
do
    echo NEWS $INDEX
    ./build/create_freq_index $INDEX $NEWS_PATH $NEWS_OUTPATH$INDEX --check &>> output_NEWS.txt
done

