# Get build info
HOST=`hostname`
DATE=`date +%F`
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
  echo ${var^^}
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
echo "Build info:"                         | tee -a $TEMP
echo "  "$HOST                             | tee -a $TEMP
echo "  "$DATE                             | tee -a $TEMP
echo "RBVMS git info:"                     | tee -a $TEMP

cd $RBVMS
pwd
RBVMS_GIT_REPO=`git remote get-url origin`
RBVMS_GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
RBVMS_GIT_COMMIT=`git rev-parse HEAD`
RBVMS_GIT_STATUS=`git status -s`

echo "  "$RBVMS_GIT_REPO                   | tee -a $TEMP
echo "  "$RBVMS_GIT_BRANCH                 | tee -a $TEMP
echo "  "$RBVMS_GIT_COMMIT                 | tee -a $TEMP
echo "  "$RBVMS_GIT_STATUS                 | tee -a $TEMP

if $HYPRE; then
  echo "HYPRE git info:"                     | tee -a $TEMP
  cd $RBVMS/../hypre
  pwd
  HYPRE_GIT_REPO=`git remote get-url origin`
  HYPRE_GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
  HYPRE_GIT_COMMIT=`git rev-parse HEAD`
  HYPRE_GIT_STATUS=`git status -s`

  echo "  "$HYPRE_GIT_REPO                   | tee -a $TEMP
  echo "  "$HYPRE_GIT_BRANCH                 | tee -a $TEMP
  echo "  "$HYPRE_GIT_COMMIT                 | tee -a $TEMP
  echo "  "$HYPRE_GIT_STATUS                 | tee -a $TEMP

fi

if $MFEM; then
  echo "MFEM git info:"                      | tee -a $TEMP
  cd $RBVMS/../mfem
  pwd
  MFEM_GIT_REPO=`git remote get-url origin`
  MFEM_GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
  MFEM_GIT_COMMIT=`git rev-parse HEAD`
  MFEM_GIT_STATUS=`git status -s`

  echo "  "$MFEM_GIT_REPO                    | tee -a $TEMP
  echo "  "$MFEM_GIT_BRANCH                  | tee -a $TEMP
  echo "  "$MFEM_GIT_COMMIT                  | tee -a $TEMP
  echo "  "$MFEM_GIT_STATUS                  | tee -a $TEMP
fi

echo "------------------------------------"| tee -a $TEMP
echo ')~");'                               >>$TEMP

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
