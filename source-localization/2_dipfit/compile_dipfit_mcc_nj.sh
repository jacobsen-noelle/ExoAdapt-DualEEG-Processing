#!/bin/sh
# Submodules that you want to add to compile your function
SUBMODS_IN=("/blue/dferris/jacobsen.noelle/fieldtrip-20210614"
            "/blue/dferris/jacobsen.noelle/eeglab2022.0/")
# The location of the function
FILE_IN=/blue/dferris/jacobsen.noelle/CustomHeadModel/2_dipfit
# The location for the output from the compiler
FILE_OUT=/blue/dferris/jacobsen.noelle/CustomHeadModel/2_dipfit/compiled
# Check to make sure directories exist
if test -f $FILE_OUT; then
    echo "$FILE_OUT exists."
else
    mkdir $FILE_OUT
fi

# purge modules currently loaded and load matlab compiler resources
module purge
module load matlab/2020a

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"

cd $FILE_IN

mcc -m dipfit_customHeadModel_hpg.m -d $FILE_OUT -a $SUBMODS_IN
