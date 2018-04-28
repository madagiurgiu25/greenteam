
file=$1
out="filtered.txt"

#cat "$file" | awk -F"\t" '$7 !~ /Data is not available/ {print}' | grep -iP "(plaque*)|(athero*)|(endothelial*)|(lipo*)|(macrophag*)|(vascular*)|(aorta*)|(immune*)|(muscle*)" > tmp.txt
cat "$file" | grep -iP "(plaque*)|(athero*)|(endothelial*)|(lipo*)|(macrophag*)|(vascular*)|(aorta*)|(immune*)|(muscle*)" > tmp.txt
