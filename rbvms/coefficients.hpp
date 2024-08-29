// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#ifndef RBVMS_COEFF_HPP
#define RBVMS_COEFF_HPP

#include "mfem.hpp"

using namespace mfem;
using namespace std;

/// class for C-function coefficient in seperate library
class LibCoefficient : public Coefficient
{
protected:
   typedef real_t (*TDFunPtr)(real_t *, int, real_t);
   TDFunPtr TDFunction;
   void *libHandle;
   Vector x;
   real_t val;

public:
   /// Define a time-independent coefficient from a C-library
   LibCoefficient(string libName, string funName,
                  bool required = true, real_t val = 0.0);
   LibCoefficient(string libName, vector<string> funNames,
                  bool required = true, real_t val = 0.0);

   /// Evaluate coefficient
   virtual real_t Eval(ElementTransformation &T,
                       const IntegrationPoint &ip) override;
   /// Destructor
   ~LibCoefficient();
};

/// class for C-function coefficient in seperate library
class LibVectorCoefficient : public VectorCoefficient
{
protected:
   typedef void (*TDFunPtr)(real_t *, int, real_t, real_t *, int);
   TDFunPtr TDFunction;
   void *libHandle;
   Vector x;

public:
   /// Define a time-independent coefficient from a C-library
   LibVectorCoefficient(int vdim, string libName, string funName, bool required = true);
   LibVectorCoefficient(int vdim, string libName, vector<string> funNames, bool required = true);

   /// Evaluate coefficient
   virtual void Eval(Vector &V,
                     ElementTransformation &T,
                     const IntegrationPoint &ip) override;
   /// Destructor
   ~LibVectorCoefficient();
};



#endif
