# Get build info
PWD=`pwd`
TEMP=$PWD/.tmp.buildInfo.hpp

MFEM=false
HYPRE=false

if [ $# -le 0 ]; then
  echo "Illegal number of inputs"
  echo "Provide at least the RBVMS source directory"
  exit -1
fi
RBVMS=$1
shift
for var in $@; do
  if [ ${var^^} == HYPRE ]; then
     HYPRE=true
  fi
  if [ ${var^^} == MFEM ]; then
     MFEM=true
  fi
done

# Create new hpp file from git info
echo "#include <sstream>"                  > $TEMP

echo 'std::istringstream buildInfo(R"~('   >>$TEMP
echo "------------------------------------"| tee -a $TEMP
echo " - Build info:"                      | tee -a $TEMP
hostname                                   | tee -a $TEMP
date +%F                                   | tee -a $TEMP

echo " - RBVMS git info:"                  | tee -a $TEMP
cd $RBVMS
pwd | tee -a $TEMP
git remote get-url origin | tee -a $TEMP
git rev-parse --abbrev-ref HEAD | tee -a $TEMP
git rev-parse HEAD | tee -a $TEMP
git status -sb | tee -a $TEMP

if $HYPRE; then
  echo " - HYPRE git info:"| tee -a $TEMP
  cd $RBVMS/../hypre
  pwd | tee -a $TEMP
  git remote get-url origin | tee -a $TEMP
  git rev-parse --abbrev-ref HEAD | tee -a $TEMP
  git rev-parse HEAD | tee -a $TEMP
  git status -sb | tee -a $TEMP

fi

if $MFEM; then
  echo " - MFEM git info:"| tee -a $TEMP
  cd $RBVMS/../mfem
  pwd | tee -a $TEMP
  git remote get-url origin | tee -a $TEMP
  git rev-parse --abbrev-ref HEAD | tee -a $TEMP
  git rev-parse HEAD | tee -a $TEMP
  git status -sb | tee -a $TEMP
fi

echo "------------------------------------"| tee -a $TEMP
echo ')~");'>>$TEMP

# Replace hpp if necessary
cd $RBVMS
if test -f buildInfo.hpp; then
  cmp -s $TEMP buildInfo.hpp || mv $TEMP buildInfo.hpp
else
  mv $TEMP buildInfo.hpp
fi

# Cleanup tmp file
rm -f $TEMP

ls -ltrh $PWD/buildInfo.hpp
