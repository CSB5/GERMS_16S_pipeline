#!/bin/bash

blastresult="$1"
test -z "$blastresult" && { echo "ERROR: Need BLAST result file as argument (format -m8 or 9)" 1>&2; exit 1; }
test -f "$blastresult" || { echo "ERROR: File $blastresult does not exist"; exit 1; }

# idea is as follows
# remove comments lines
# print only query, subject, identity, alnlength, evalue and bitscore
# use negative alnlength and bitscore since evalue is most important and has a reversed order
# sort according to evalue, bit-score, alnlength in that order (using scientific ordering -g)
# then print best hit or hits with the one best score
#



grep -v '^#' $blastresult | \
    awk '{print $1, $2, -$3, -$4, $11, -$12}' | \
    sort -g -k 1 -k 5 -k 6 -k 3 -k 4 | \
    awk 'BEGIN {print "# query best-hit"; scorestr_last_seen="EMPTY"; query_last_seen="EMPTY"}
         {query=$1; scorestr=sprintf("%s%s%s%s%s", $1, $3, $4, $5, $6);
          if (query != query_last_seen) {print $1, $2, $3, $4, $5, $6; query_last_seen=query; scorestr_last_seen=scorestr}}'

# if you want to keep all equally good hits use this instead
#          if (query != query_last_seen || scorestr == scorestr_last_seen) {print $1, $2; query_last_seen=query; scorestr_last_seen=scorestr}}'
 
