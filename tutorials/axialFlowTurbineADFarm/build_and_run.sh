#!/bin/bash

# Save the original directory
ORIGINAL_DIR=$(pwd)

# Go up two levels
cd ../.. || {
    echo "Could not go up two levels. Aborting."
    exit 1
}

echo "Running Allwmake in $(pwd) ..."
./Allwmake
MAKE_STATUS=$?

# Return to the original directory regardless of the result
cd "$ORIGINAL_DIR" || {
    echo "Could not return to the original directory. Aborting."
    exit 1
}

if [ $MAKE_STATUS -eq 0 ]; then
    echo "Allwmake succeeded. Running Allclean and Allrun in $ORIGINAL_DIR ..."
    ./Allclean
    time ./Allrun "$@"
else
    echo "Allwmake failed (exit code $MAKE_STATUS). Skipping Allclean/Allrun."
    exit $MAKE_STATUS
fi
