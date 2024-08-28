# Get build info
HOST=`hostname`
DATE=`date +%F`



pushd $1
RBVMS_GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
RBVMS_GIT_COMMIT=`git rev-parse HEAD`
RBVMS_GIT_STATUS=`git status -s`
popd

pushd $2
MFEM_GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
MFEM_GIT_COMMIT=`git rev-parse HEAD`
MFEM_GIT_STATUS=`git status -s`
popd


# Print build info
echo "------------------------------------"
echo "Build info:"
echo  "  "$HOST
echo  "  "$DATE
echo "RBVMS git info:"
echo  "  "$RBVMS_GIT_BRANCH
echo  "  "$RBVMS_GIT_COMMIT
echo  "  "$RBVMS_GIT_STATUS
echo "MFEM git info:"
echo  "  "$MFEM_GIT_BRANCH
echo  "  "$MFEM_GIT_COMMIT
echo  "  "$MFEM_GIT_STATUS
echo "------------------------------------"

# Create new hpp file from git info
echo "#include <string>"                   > .tmp.buildInfo.hpp
echo "#include <fstream>"                  >>.tmp.buildInfo.hpp
echo "#include <sstream>"                  >>.tmp.buildInfo.hpp

echo 'std::istringstream buildInfo(R"~('   >>.tmp.buildInfo.hpp
echo "------------------------------------">>.tmp.buildInfo.hpp
echo "Build info:"                         >>.tmp.buildInfo.hpp
echo "  "$HOST                             >>.tmp.buildInfo.hpp
echo "  "$DATE                             >>.tmp.buildInfo.hpp
echo "RBVMS git info:"                     >>.tmp.buildInfo.hpp
echo "  "$RBVMS_GIT_BRANCH                 >>.tmp.buildInfo.hpp
echo "  "$RBVMS_GIT_COMMIT                 >>.tmp.buildInfo.hpp
echo "  "$RBVMS_GIT_STATUS                 >>.tmp.buildInfo.hpp
echo "MFEM git info:"                      >>.tmp.buildInfo.hpp
echo "  "$MFEM_GIT_BRANCH                  >>.tmp.buildInfo.hpp
echo "  "$MFEM_GIT_COMMIT                  >>.tmp.buildInfo.hpp
echo "  "$MFEM_GIT_STATUS                  >>.tmp.buildInfo.hpp
echo "------------------------------------">>.tmp.buildInfo.hpp
echo ')~");'                               >>.tmp.buildInfo.hpp

# Replace hpp if necessary
pwd
if test -f buildInfo.hpp; then
  cmp -s  .tmp.buildInfo.hpp buildInfo.hpp || mv .tmp.buildInfo.hpp buildInfo.hpp
else
  mv .tmp.buildInfo.hpp buildInfo.hpp
fi

# Cleanup tmp file
rm -f .tmp.buildInfo.hpp 



