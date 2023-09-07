#!/bin/bash

# Function to checkpoint the pipeline
checkpoint() {
    echo "Checkpoint: $1"
    echo "$1" > checkpoint.txt
}

# Function to handle errors
handle_error() {
    echo "Error occurred. Checkpoint: $1"
    exit 1
}

# Check if a checkpoint file exists, indicating a previous run
if [ -e checkpoint.txt ]; then
    last_checkpoint=$(cat checkpoint.txt)
    echo "Resuming from checkpoint: $last_checkpoint"
else
    last_checkpoint="start"
fi

# Main pipeline
case "$last_checkpoint" in
    "start")
        echo "Starting the pipeline..."
        # Your code for the first step here
        if [ $? -ne 0 ]; then
            handle_error "step1"
        fi
        checkpoint "step1"

    # Add more steps here as needed, using the same pattern
    "step1")
        echo "Continuing from step1..."
        # Your code for the second step here
        if [ $? -ne 0 ]; then
            handle_error "step2"
        fi
        checkpoint "step2"

    "step2")
        echo "Continuing from step2..."
        # Your code for the third step here
        if [ $? -ne 0 ]; then
            handle_error "step3"
        fi
        checkpoint "step3"

    "step3")
        echo "Pipeline completed successfully!"
        # Clean up the checkpoint file
        rm checkpoint.txt
        ;;

    *)
        echo "Unknown checkpoint: $last_checkpoint"
        handle_error "unknown"
        ;;
esac