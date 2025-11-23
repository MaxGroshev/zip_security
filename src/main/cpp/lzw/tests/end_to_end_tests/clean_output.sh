#!/bin/bash

rm -r $(dirname $0)/my_decompress_dat/*.decode
rm -r $(dirname $0)/my_compress_dat/*.lzw
rm -r $(dirname $0)/my_compress_dat/*.zstd
rm -r $(dirname $0)/my_compress_dat/*.zst

