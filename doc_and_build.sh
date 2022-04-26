#!/bin/bash

# This script will generate a new 
# html documentation as well as 
# build a new package binary

# ----------------------------------------
#  Make an HTML documentation using pydoc
# ----------------------------------------

# make a docs folder if not yet present
if [[ ! -d ./docs ]]; then
    mkdir -p ./docs
fi

# run pdoc 
pdoc ./EEGToolkit/ --output-dir ./docs --html --force

echo "Documentation generated successfully!"

# ----------------------------------------
#             Build a new binary
# ----------------------------------------


# first make new requirements
pipreqs ./EEGToolkit --force --savepath ./requirements.txt
echo "Requirements generated successfully!"

# now build
python -m build .
echo "Binary generated successfully!"
