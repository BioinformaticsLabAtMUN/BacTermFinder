# All the bed files for one bacteria have the same reference genome so the name should've been the same. 
# I had underscore index for each line which had to get removed

# Step 1
for dir in $(ls -d */)
do
    # go through each dir
    cd $dir
    # go through each .bed file
    for file in *.bed
    do
        # check that file name has ".bed"
        if [[ $file == *.bed ]]
        then
            echo $file
            # delete the underscore and what's after that before tab in the bed file column 1 for each line
            sed -i 's|\(\.[0-9]\)_[0-9]\+	|\1	|g' $file
        fi
    done
    # go back to the parent directory
    cd ..
done