# HyperMut

This repository contains a program to identify sequences
that have multiple T-to-C or G-to-A mutations in a fragment.

Primary inputs are a sam file and a reference sequence.
This program have several configurable parameters including
the quality threshold, maximum number of undetermined nucleotides,
threshold to count as "many" mutations, and specify the range
of sequence that should be considered and ignored due to PCR primers.

An example would look like:

  ruby sammmcontext-f4.rb --ref=ref.fa --sam=sample.sam \
    --min-qv=60 --max-N=3 -m 600 -t 2 \
    -p 1 -e 4264 -c 27 -d 4246 -g 3699 -h 3911 \
    > sample.refh.ighsubstcuorga.out60c

See the program for detail of the parameters.

