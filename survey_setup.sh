#!/bin/bash
# Description: This script creates the startup folders to
# Author: Angel-Emilio Villegas Sanchez

# Source directory for files you want to place in your slurm job directories
src="/mnt/home/k0052095/organic_survey/src_data"
extradir="/mnt/home/k0052095/organic_survey/extraction2"
success="/mnt/home/k0052095/organic_survey/test_success"
# Source of SMILES structures
molecules=$src/7.cls.smi
# Get the date-time in order to show log
current_datetime=$(date +"%m-%d-%y %H:%M")
printf "$current_datetime:\n" >> SURVEY_LOG
# Indices to loop over from SMILES structures.
indeces=$(seq 102 3 132)
# Declare groups and processors per group for your FLOSIC calculation
groups=15
procs=25
# SLURM SETTINGS
time="6:30:00"

InitialDFT () {
        # Main Loop
        # Finding all directories in file. ls returns them with a /
        for index in $(ls -d */ | sed 's/\///g')
        do

                mol=$(awk "NR==$index" $molecules)
                # Check if folder exists, if not create it
                if [ ! -f $index/RUNS ]; then
                        cd $index

                        # Copy/create startup files
                        cp $src/survey_batch.sh .
                        cp $src/NRLMOL_INPUT.DAT .

                        # Generate CLUSTER and FRMORB
                        ~/.local/bin/FODLego/findfods.py $mol
                        mv FRMORB tmpfo

                        # Execute FLOSIC DFT RUN
                        sbatch -n 28 -J "${index}_dft" -t 0:15:00 survey_batch.sh
                        printf "Submitting in $index\n" >> ../SURVEY_LOG

                        # Leave directory and go to next one
                        cd ..
                else
                        #RedoDFT?
                        echo "There is a RUNS file already here."
                fi
        done
}

# Description: Find folders that have an unfinished RUNS file. Redo the submission for DFT.
RedoDFT() {
        # Main Loop
        for index in $(ls -d */)
        do
                cd $index
                        runs=$(awk 'NR==2 {print $1}' RUNS)
                        # Check RUNS output
                        case $runs in
                                3)
                                        msg="Attempting a DFT rerun at $index"
                                        echo "$msg"
                                        printf "$msg\n" >> ../SURVEY_LOG
                                        sbatch -n 28 -J "${index}" -t 0:30:00 survey_batch.sh
                                        ;;
                        esac
                cd ..
        done
}

# Description: Check if fande.dat exists, if so, check for convergence
fande_converge() {
        if [ -f "fande.dat" ]; then
                last_bytes=$(tail -c 2 fande.dat)
                if (( $last_bytes == 4 )); then
                        return 0 # Condition met, return true
                else
                        return 1 # Condition not met, return false
                fi
        else
                return 2 # File does not exist, return false
        fi
}

# Description: Find out if there is a pending or ongoing slurm batch in this folder
ongoing_batch() {
        local folder_name="$1"
        local slurm_list=$(squeue -u $USER | awk 'NR >= 2 {print $3}')

        for dir in $slurm_list; do
                if [ "$dir" == "$folder_name" ]; then
                        return 0 # Folder found in the list
                fi
        done

        return 1 # Folder not found in the list
}

submit_flo() {
        # Change the igroup file
        tasks=$(($groups*$procs))
        echo -e "$groups\n$procs" > igroup

        # Execute FLOSIC if there is no E-04 convergence
        if ongoing_batch $index
        then
                echo "$index: Folder found in SLURM RUNNINg/PENDING list."
        else
                log="Running a new FLOSIC job at $index"
                printf "$log\n" >> ../SURVEY_LOG
                sbatch -n $tasks -J "${index}" -t $time survey_batch.sh
        fi
}

# Description: Submit jobs with FLOSIC calculations
FLOSIC() {
        for index in $(ls -d */)
            do
                    # Check if folder exists, if not create it
                    if [ -d $index ]; then

                                        # Step into directory
                                        cd $index
                                        runs=$(awk 'NR==2 {print $1}' RUNS)
                                        autoconv=" ************************************************"
                                        checklog=$(sed -n '$ p' log)

                                        # Check the convergence
                                        if [[ 4 == $runs && $checklog != $autoconv ]]; then
                                                            # Change FRMORB name to execute FLOSIC
                                                            if [ -f tmpfo ]; then
                                                                                mv tmpfo FRMORB
                                                            fi

                                                            submit_flo

                                        else
                                                echo "Folder $index does not have wavefunctions. DFT must be finished."
                                                echo "Or, the last calculation was succesful and fully converged"
                                        fi
                                        # Leave directory and go to next one
                                        cd ..
                    else
                            echo "Folder $index does not exist."
                    fi
        done
}


CHECKFLO() {
    tally=0
    for index in $(ls -d */)
        do
                if [ -d $index ]; then
                # Step into directory
                cd $index
                # Get the current date and time
                if [ -f "fande.out" ]; then
                    # Count the number of entries matching the pattern in fande.out
                    num_entries=$(grep -Ec "(-04)" fande.out)

                    # Log the current date and number of successful convergences
                    log_entry="$num_entries successful convergences in  $index\n"

                    # Append log entry to SURVEY_LOG file
                    printf "$log_entry" >> ../SURVEY_LOG

                    # Output the number of successful convergences
                    echo "Number of successful convergences: $num_entries in $index"

                    # Create more runs if you only have one
                    if [[ $num_entries -le 3 ]]; then
                        echo "Submitting batch in $index"
                        submit_flo
                    fi
                else
                        echo "fande.out does not exist in $index."
                                            log="REDO calculation in $index\n"
                                            printf "$log" >> ../SURVEY_LOG
                                    fi
                                    cd ..
                    else
                            echo "Folder $index does not exist."
                    fi
        done
}

EXTRACT(){
	mkdir $extradir
        for dir in $(ls); do
                fande="$dir/fande.out"
                if [ -f $fande ]; then
                            # Get counts of
                            num_entries=$(grep -Ec "(-04)" $fande)
                            if [ $num_entries -gt 0 ]; then
                                        echo $dir
                                        # Copy relevant files!
                                        cp -n $dir/XMOL.xyz $extradir/${dir}_XMOL.xyz
                                        cp -n $dir/FRMORB $extradir/${dir}_FRMORB
                                        SMI=$(sed -n '2p' ${dir}/out.xyz)
                                        printf "${SMI}\textraction/${dir}_FRMORB\textraction/${dir}_XMOL.xyz\n" >> files_2nd_extraction
                            fi
                else
                        echo "No FLOSIC calculation found at $dir"
                fi
        done
}

MOVESUCCESS(){
        for dir in $(ls); do
                fande="$dir/fande.out"
                if [ -f $fande ]; then
                            # Get counts of
                            num_entries=$(grep -Ec "(-04)" $fande)
                            if [ $num_entries -gt 0 ]; then
                                        echo "Moving $dir to test_success"
                                        # Copy relevant files!
                                        mv $dir $success
                            fi
                else
                        echo "No FLOSIC calculation found at $dir"
                fi
        done
}

RMINVAL(){
        FODLego="/mnt/home/k0052095/.local/bin/FODLego/findfods.py"
        mkdir trashruns
        for dir in $(ls); do
                mol=$(awk "NR==$dir" ../src_data/7.cls.smi)
                check=$($FODLego checksmiles $mol)
                if [[ $check == 'nonvalid' ]]; then
                            # Mv away from directory!
                            mv $dir trashruns/$dir
                fi
        done
}

CHECKVALIDRUNS(){
        FODLego="/mnt/home/k0052095/.local/bin/FODLego/findfods.py"
        for mol in $(cat $molecules); do
                check=$($FODLego checksmiles $mol)
                if [[ $check == 'nonvalid' ]]; then
                            # Mv away from directory!
                            echo "$mol" >> goodstructs
                fi
        done
}

CHECKDFT() {
        for index in $(ls -d */)
        do
                # Check if folder exists, if not create it
                if [ -d $index ]; then
                        # Step into directory
                        cd $index
                        runs=$(awk 'NR==2 {print $1}' RUNS)
                        # Check RUNS output
                        case $runs in
                                3)
                                        gp=$(grep "REORMSH" log)

                                        # Check if $gp is not empty
                                        if [ -n "$gp" ]; then
                                                # Print message to SURVEY_LOG
                                                message="$index: $gp\n"
                                                echo -e $message
                                                printf  "$message" >> ../SURVEY_LOG
                                        fi
                                        log="Wavefunctions not finished in $index"
                                        echo -E "$log"
                                        printf "$log\n" >> ../SURVEY_LOG
                                        ;;
                                4)
                                        echo "Wavefunctions finished at $index"
                                        ;;
                        esac
                        cd ..
                fi
        done
}
# Function to log line 7 of survey_setup.sh file
# Written by ChatGPT. Prompted by AEVS
logsurvey() {
        # Append log entry to SURVEY_LOG file
        log_entry=$(basename $PWD)
        printf "$log_entry\n" >> SURVEY_LOG
}

GETRECORDS(){
        # From success, get the records with their prefix index
        if [ ! -d getrecords ]; then
                    mkdir getrecords/
        fi
        for struct in $(ls -d $success/*/); do
                    base=$(basename $struct)
                    cp $struct/records getrecords/${base}_record
        done
}

# This is the meat of this script. Based of user input, do the DFT or FLOSIC loop
case $# in
        0)
                echo "No arguments passed" ;;
        1)
                case $1 in
                        DFT)
                                echo Creating initial DFT Loop
                                InitialDFT
                                logsurvey
                                ;;
                        RedoDFT)
                                echo "Redoing DFT in failed directories"
                                RedoDFT
                                ;;
                        FLO)
                                echo Running FLOSIC calculations
                                FLOSIC
                                logsurvey
                                ;;
                        CHECKFLO)
                                CHECKFLO
                                ;;
                        CHECKDFT)
                                CHECKDFT
                                ;;
                        MOVESUCCESS)
                                MOVESUCCESS
                                ;;
                        RMINVAL)
                                RMINVAL
                                ;;
                        EXTRACT)
                                EXTRACT
                                ;;
                        CHECKVALIDRUNS)
                                CHECKVALIDRUNS
                                ;;
                        GETRECORDS)
                                GETRECORDS
                                ;;
                esac
esac
