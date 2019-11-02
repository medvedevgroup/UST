DBGFM_DIRECTORY=/home/user/dbgfm
BWTDISK_DIRECTORY=/home/user/bwtdisk

f=$1

cat $f > `basename $f .fa`.pp.fa

BASE=`basename $f .fa`
PP=`basename $f .fa`.pp.fa
echo "Running bwtdisk........"

# Prepare the file by concatenating the contigs/sequences into one long string
# separated by $
$DBGFM_DIRECTORY/bwtdisk-prepare $PP > $PP.joined

# Reverse the text, because bwtdisk constructs the BWT of the reverse text by default.
$BWTDISK_DIRECTORY/text_rev -vvv $PP.joined $PP.joined.rev

# Generate the BWT.
$BWTDISK_DIRECTORY/bwte -vvv $PP.joined.rev

# Rename the bwt
mv $PP.joined.rev.bwt $BASE.pp.bwtdisk
echo "Done converting."

# To query dbgfm
#DBGFM_DIRECTORY/dbgfm $BASE.pp
