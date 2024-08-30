#include "coefficients.hpp"
#include <dlfcn.h>

using namespace std;
using namespace mfem;

/// Constructor
LibCoefficient::LibCoefficient(string libName,
                               string funName,
                               bool required,
                               real_t val) : val(val)
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
      cout<<"Function will be set to default"<<endl;
      cout<<funName<<" = "<<val<< endl;
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
      cout<<"Function will be set to default"<<endl;
      cout<<funName<<" = "<<val<< endl;
   }
}

/// class for C-function coefficient in seperate library
LibCoefficient::LibCoefficient(string libName,
                               vector<string> funNames,
                               bool required,
                               real_t val) : val(val)
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
      cout<<" will be set to val = "<<val<< endl;
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
      cout<<"Function will be set to val = "<<val<< endl;
   }
}

/// Evaluate coefficient
real_t LibCoefficient::Eval(ElementTransformation &T,
                            const IntegrationPoint &ip)
{
   // Homegenous if not defined
   if (!TDFunction) return val;

   // Evaluate library function
   T.Transform(ip, x);
   return ((*TDFunction)(x.GetData(),x.Size(),GetTime()));
}

/// Destructor
LibCoefficient::~LibCoefficient()
{
   if (libHandle) dlclose(libHandle);
}


/// Constructor
LibVectorCoefficient::LibVectorCoefficient(int vdim,
                                           string libName,
                                           string funName,
                                           bool required) 
  : VectorCoefficient(vdim)
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
      cout<<"Function "<<funName<<
            " will be set to be homogenous."<<endl;
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
      cout<<"Function set to be homogenous."<<endl;
   }
}

/// class for C-function coefficient in seperate library
LibVectorCoefficient::LibVectorCoefficient(int vdim,
                                           string libName,
                                           vector<string> funNames,
                                           bool required)
  : VectorCoefficient(vdim)
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
void LibVectorCoefficient::Eval(Vector &V,
                                ElementTransformation &T,
                                const IntegrationPoint &ip)
{
   // Homegenous if not defined
   if (!TDFunction)
   {
      V = 0.0;
      return;
   }

   // Evaluate library function
   T.Transform(ip, x);
   (*TDFunction)(x.GetData(), x.Size(), GetTime(), V.GetData(), V.Size());
}

/// Destructor
LibVectorCoefficient::~LibVectorCoefficient()
{
   if (libHandle) dlclose(libHandle);
}
