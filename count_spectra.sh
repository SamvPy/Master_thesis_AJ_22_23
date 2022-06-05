#! /bin/sh
gzip -dk -S .gzip -r *
grep -r --include \*.mgf -c "^BEGIN IONS"  >>output.txt
find . -name "*.mgf" -type f -delete
find . -name "*.gzip" -type f -delete
find . -name "*.gz" -type f -delete