f=$1
cat $f > `basename $f .fa`.pp.fa	

BASE=`basename $f .fa`
PP=`basename $f .fa`.pp.fa
echo "Running bwtdisk........"
	
# Prepare the file by concatenating the contigs/sequences into one long string
# separated by $
dbgfm/bwtdisk-prepare $PP > $PP.joined

# Reverse the text, because bwtdisk constructs the BWT of the reverse text by default.
bwtdisk/text_rev -vvv $PP.joined $PP.joined.rev

# Generate the BWT.
bwtdisk/bwte -vvv $PP.joined.rev

# Rename the bwt
mv $PP.joined.rev.bwt $BASE.pp.bwtdisk
echo "Done converting."

# To query dbgfm
#dbgfm/dbgfm $BASE.pp 
