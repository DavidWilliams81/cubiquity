#!/usr/bin/env bash
set -euo pipefail # exit on errors, unset vars, or failed pipes

mkdir temp
trap "rm -rf temp" EXIT # Delete temp on exit

init_test() {
    mkdir temp/src # For source data
    mkdir temp/dst # For exported data

    # Generate the source dag
    ./cubiquity generate fractal_noise \
        --size=100 \
        --output=temp/src/fractal_noise.dag
}

# --- Define tests ------------------------------------------------------

# Create a tar file containing a compressed stream and uncompressed metadata.
# The user can extract this without taking up excessive space (because the
# stream is compressed), and then pipe the stream through a decompressor and
# into their application.
test_tar_with_compressed_bin() {
    init_test

    # Write the bin file to stdout (and compress with
    # gzip) while writing the metadata as a normal file.
    # FIXME - This should error if temp doesn't exist but it silently fails
    ./cubiquity export bin temp/src/fractal_noise.dag \
        --output=- \
        --output-metadata=temp/dst/fractal_noise.txt \
        --verbose \
            | gzip > temp/dst/fractal_noise.bin.gz \

    # Run tar from inside the temp/dst folder (-C).
    tar -cvf fractal_noise.tar -C temp/dst fractal_noise.bin.gz fractal_noise.txt
}

test_zip() {
    init_test

    # We populate out zip file from named pipes (FIFOs), as trying to pass the
    # data via stdin is a bit of a minefield. The zip tool expects that stdin
    # input is a list of files to compress (rather than file contents), and if
    # you get it to accept it as file contents it calls the resulting file '-'.
    # This can be renamed via 'zipnote', but at the expense of recompressing the
    # archive (I believe). Also you then have to add the text file, as only one
    # file can go in via stdin. Named pipes are a much cleaner solution.
    #
    # Other tools such as bsdtar and 7zip have similar issues - stdin is just
    # not well suited for populating multifile archives. And these don't seem to
    # support even named pipes for .zip files
    mkfifo temp/dst/fractal_noise.bin temp/dst/fractal_noise.txt

    # Start exporting to named pipes (will block until reader is active)
    ./cubiquity export bin temp/src/fractal_noise.dag \
        --output=temp/dst/fractal_noise.bin \
        --output-metadata=temp/dst/fractal_noise.txt \
        --verbose &

    # Sleep for testing purposes (we see Cubiquity blocks waiting for reader).
    sleep 2

    # Read fromm the named pipes (fifos) into zip file while discarding paths
    zip --junk-paths --fifo fractal_noise.zip \
        temp/dst/fractal_noise.bin temp/dst/fractal_noise.txt
}

# --- Menu --------------------------------------------------------------

echo "Select a test to run:"
echo "  1) Export .tar file (compressed bin with uncompressed metadata)"
echo "  2) Export .zip file (both bin and metadata are compressed)"
echo "  q) Quit"
echo

read -rp "Enter choice: " choice
echo

case "$choice" in
	1) test_tar_with_compressed_bin ;;
	2) test_zip ;;
	q|Q) echo "Exiting."; exit 0 ;;
	*) echo "Invalid choice." >&2; exit 1 ;;
esac
