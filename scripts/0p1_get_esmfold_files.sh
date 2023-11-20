#!/bin/bash
cd .. # file assumes we are in the project-benchmark-seq2struct directory
mkdir ./pepstructure/esmfold_v1/

huggingface_url="https://huggingface.co/facebook/esmfold_v1/resolve/main/"

for file in "README.md" "config.json" "special_tokens_map.json" "tokenizer_config.json" "vocab.txt"; do
	wget -P ./pepstructure/esmfold_v1/ $huggingface_url$file
done
