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
class LibFunctionCoefficient : public Coefficient
{
protected:
   typedef double (*TDFunPtr)(double *, int, double);
   TDFunPtr TDFunction;
   void *libHandle;
   Vector x;



public:
   /// Define a time-independent coefficient from a C-library
   LibFunctionCoefficient(string libName, string funName, bool required = true);
   LibFunctionCoefficient(string libName, vector<string> funNames, bool required = true);

   /// Evaluate coefficient
   virtual double Eval(ElementTransformation &T,
                       const IntegrationPoint &ip) override;
   /// Destructor
   ~LibFunctionCoefficient();
};

#endif
