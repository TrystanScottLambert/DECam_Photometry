# Script to reduce all .sed files to below 25 000 angstrom. Anything else would be unphysical

rm *trimmed.sed
for file in *.sed; do
    awk '$1 < 25000' "$file" > "${file%.sed}_trimmed.sed"
done
