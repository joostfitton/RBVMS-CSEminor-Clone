#include "coefficients.hpp"
#include <dlfcn.h>

using namespace std;
using namespace mfem;

/// Constructor
LibFunctionCoefficient::LibFunctionCoefficient(string libName,
                                               string funName,
                                               bool required)
{
   // Open library
   libHandle = dlopen (libName.c_str(), RTLD_LAZY);
   if (!libHandle)
   {
      cout<<"Library "<<libName<<" not found."<<endl;
      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout<<"Function "<<funName<<" will be set to homogenous."<<endl;
      return;
   }

   // Get Function
   TDFunction = (TDFunPtr)dlsym(libHandle, funName.c_str());
   if (!TDFunction)
   {
      cout<<"Function "<<funName<<" not found in "<<libName<<endl;
      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout<<"Function will be set to homogenous."<<endl;
   }
}

/// class for C-function coefficient in seperate library
LibFunctionCoefficient::LibFunctionCoefficient(string libName,
                                               vector<string> funNames,
                                               bool required)
{
   // Open library
   libHandle = dlopen (libName.c_str(), RTLD_LAZY);
   if (!libHandle)
   {
      cout<<"Library "<<libName<<" not found."<<endl;
      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout <<"Functions:";
      for (int i = 0; i < funNames.size();i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" will be set to homogenous."<<endl;
      return;
   }

   // Search Function
   TDFunction = nullptr;
   for (int i = 0; i < funNames.size();i++)
   {
      TDFunction = (TDFunPtr)dlsym(libHandle, funNames[i].c_str());
      if (TDFunction) break;
   }

   // Check if function is found
   if (!TDFunction)
   {
      cout <<"Functions:";
      for (int i = 0; i < funNames.size();i++)
      {
         cout<<" "<<funNames[i];
      }
      cout<<" can not be found in "<<libName<<endl;

      if (required)
      {
         mfem_error("Can not obtain required function.\n");
      }
      cout<<"Function will be set to homogenous."<<endl;
   }
}

/// Evaluate coefficient
double LibFunctionCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
   if (!TDFunction) return 0.0;

   T.Transform(ip, x);
   return ((*TDFunction)(x.GetData(),x.Size(),GetTime()));
}

/// Destructor
LibFunctionCoefficient::~LibFunctionCoefficient()
{
   if (libHandle) dlclose(libHandle);
}
