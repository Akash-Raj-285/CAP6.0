// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__Jets
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "JetHistos.hpp"
#include "JetSingleHistos.hpp"
#include "JetPairHistos.hpp"
#include "JetAnalyzer.hpp"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_CAPcLcLJetHistos(void *p = nullptr);
   static void *newArray_CAPcLcLJetHistos(Long_t size, void *p);
   static void delete_CAPcLcLJetHistos(void *p);
   static void deleteArray_CAPcLcLJetHistos(void *p);
   static void destruct_CAPcLcLJetHistos(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CAP::JetHistos*)
   {
      ::CAP::JetHistos *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CAP::JetHistos >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("CAP::JetHistos", ::CAP::JetHistos::Class_Version(), "JetHistos.hpp", 16,
                  typeid(::CAP::JetHistos), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CAP::JetHistos::Dictionary, isa_proxy, 4,
                  sizeof(::CAP::JetHistos) );
      instance.SetNew(&new_CAPcLcLJetHistos);
      instance.SetNewArray(&newArray_CAPcLcLJetHistos);
      instance.SetDelete(&delete_CAPcLcLJetHistos);
      instance.SetDeleteArray(&deleteArray_CAPcLcLJetHistos);
      instance.SetDestructor(&destruct_CAPcLcLJetHistos);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CAP::JetHistos*)
   {
      return GenerateInitInstanceLocal(static_cast<::CAP::JetHistos*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CAP::JetHistos*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CAPcLcLJetSingleHistos(void *p = nullptr);
   static void *newArray_CAPcLcLJetSingleHistos(Long_t size, void *p);
   static void delete_CAPcLcLJetSingleHistos(void *p);
   static void deleteArray_CAPcLcLJetSingleHistos(void *p);
   static void destruct_CAPcLcLJetSingleHistos(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CAP::JetSingleHistos*)
   {
      ::CAP::JetSingleHistos *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CAP::JetSingleHistos >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("CAP::JetSingleHistos", ::CAP::JetSingleHistos::Class_Version(), "JetSingleHistos.hpp", 15,
                  typeid(::CAP::JetSingleHistos), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CAP::JetSingleHistos::Dictionary, isa_proxy, 4,
                  sizeof(::CAP::JetSingleHistos) );
      instance.SetNew(&new_CAPcLcLJetSingleHistos);
      instance.SetNewArray(&newArray_CAPcLcLJetSingleHistos);
      instance.SetDelete(&delete_CAPcLcLJetSingleHistos);
      instance.SetDeleteArray(&deleteArray_CAPcLcLJetSingleHistos);
      instance.SetDestructor(&destruct_CAPcLcLJetSingleHistos);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CAP::JetSingleHistos*)
   {
      return GenerateInitInstanceLocal(static_cast<::CAP::JetSingleHistos*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CAP::JetSingleHistos*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CAPcLcLJetPairHistos(void *p = nullptr);
   static void *newArray_CAPcLcLJetPairHistos(Long_t size, void *p);
   static void delete_CAPcLcLJetPairHistos(void *p);
   static void deleteArray_CAPcLcLJetPairHistos(void *p);
   static void destruct_CAPcLcLJetPairHistos(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CAP::JetPairHistos*)
   {
      ::CAP::JetPairHistos *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CAP::JetPairHistos >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("CAP::JetPairHistos", ::CAP::JetPairHistos::Class_Version(), "JetPairHistos.hpp", 13,
                  typeid(::CAP::JetPairHistos), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CAP::JetPairHistos::Dictionary, isa_proxy, 4,
                  sizeof(::CAP::JetPairHistos) );
      instance.SetNew(&new_CAPcLcLJetPairHistos);
      instance.SetNewArray(&newArray_CAPcLcLJetPairHistos);
      instance.SetDelete(&delete_CAPcLcLJetPairHistos);
      instance.SetDeleteArray(&deleteArray_CAPcLcLJetPairHistos);
      instance.SetDestructor(&destruct_CAPcLcLJetPairHistos);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CAP::JetPairHistos*)
   {
      return GenerateInitInstanceLocal(static_cast<::CAP::JetPairHistos*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CAP::JetPairHistos*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CAPcLcLJetAnalyzer(void *p = nullptr);
   static void *newArray_CAPcLcLJetAnalyzer(Long_t size, void *p);
   static void delete_CAPcLcLJetAnalyzer(void *p);
   static void deleteArray_CAPcLcLJetAnalyzer(void *p);
   static void destruct_CAPcLcLJetAnalyzer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CAP::JetAnalyzer*)
   {
      ::CAP::JetAnalyzer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CAP::JetAnalyzer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("CAP::JetAnalyzer", ::CAP::JetAnalyzer::Class_Version(), "JetAnalyzer.hpp", 19,
                  typeid(::CAP::JetAnalyzer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CAP::JetAnalyzer::Dictionary, isa_proxy, 4,
                  sizeof(::CAP::JetAnalyzer) );
      instance.SetNew(&new_CAPcLcLJetAnalyzer);
      instance.SetNewArray(&newArray_CAPcLcLJetAnalyzer);
      instance.SetDelete(&delete_CAPcLcLJetAnalyzer);
      instance.SetDeleteArray(&deleteArray_CAPcLcLJetAnalyzer);
      instance.SetDestructor(&destruct_CAPcLcLJetAnalyzer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CAP::JetAnalyzer*)
   {
      return GenerateInitInstanceLocal(static_cast<::CAP::JetAnalyzer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CAP::JetAnalyzer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace CAP {
//______________________________________________________________________________
atomic_TClass_ptr JetHistos::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JetHistos::Class_Name()
{
   return "CAP::JetHistos";
}

//______________________________________________________________________________
const char *JetHistos::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetHistos*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JetHistos::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetHistos*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetHistos::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetHistos*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetHistos::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetHistos*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace CAP
namespace CAP {
//______________________________________________________________________________
atomic_TClass_ptr JetSingleHistos::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JetSingleHistos::Class_Name()
{
   return "CAP::JetSingleHistos";
}

//______________________________________________________________________________
const char *JetSingleHistos::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetSingleHistos*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JetSingleHistos::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetSingleHistos*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetSingleHistos::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetSingleHistos*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetSingleHistos::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetSingleHistos*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace CAP
namespace CAP {
//______________________________________________________________________________
atomic_TClass_ptr JetPairHistos::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JetPairHistos::Class_Name()
{
   return "CAP::JetPairHistos";
}

//______________________________________________________________________________
const char *JetPairHistos::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetPairHistos*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JetPairHistos::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetPairHistos*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetPairHistos::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetPairHistos*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetPairHistos::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetPairHistos*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace CAP
namespace CAP {
//______________________________________________________________________________
atomic_TClass_ptr JetAnalyzer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JetAnalyzer::Class_Name()
{
   return "CAP::JetAnalyzer";
}

//______________________________________________________________________________
const char *JetAnalyzer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetAnalyzer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JetAnalyzer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetAnalyzer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetAnalyzer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetAnalyzer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetAnalyzer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CAP::JetAnalyzer*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace CAP
namespace CAP {
//______________________________________________________________________________
void JetHistos::Streamer(TBuffer &R__b)
{
   // Stream an object of class CAP::JetHistos.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(CAP::JetHistos::Class(),this);
   } else {
      R__b.WriteClassBuffer(CAP::JetHistos::Class(),this);
   }
}

} // namespace CAP
namespace ROOT {
   // Wrappers around operator new
   static void *new_CAPcLcLJetHistos(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetHistos : new ::CAP::JetHistos;
   }
   static void *newArray_CAPcLcLJetHistos(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetHistos[nElements] : new ::CAP::JetHistos[nElements];
   }
   // Wrapper around operator delete
   static void delete_CAPcLcLJetHistos(void *p) {
      delete (static_cast<::CAP::JetHistos*>(p));
   }
   static void deleteArray_CAPcLcLJetHistos(void *p) {
      delete [] (static_cast<::CAP::JetHistos*>(p));
   }
   static void destruct_CAPcLcLJetHistos(void *p) {
      typedef ::CAP::JetHistos current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CAP::JetHistos

namespace CAP {
//______________________________________________________________________________
void JetSingleHistos::Streamer(TBuffer &R__b)
{
   // Stream an object of class CAP::JetSingleHistos.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(CAP::JetSingleHistos::Class(),this);
   } else {
      R__b.WriteClassBuffer(CAP::JetSingleHistos::Class(),this);
   }
}

} // namespace CAP
namespace ROOT {
   // Wrappers around operator new
   static void *new_CAPcLcLJetSingleHistos(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetSingleHistos : new ::CAP::JetSingleHistos;
   }
   static void *newArray_CAPcLcLJetSingleHistos(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetSingleHistos[nElements] : new ::CAP::JetSingleHistos[nElements];
   }
   // Wrapper around operator delete
   static void delete_CAPcLcLJetSingleHistos(void *p) {
      delete (static_cast<::CAP::JetSingleHistos*>(p));
   }
   static void deleteArray_CAPcLcLJetSingleHistos(void *p) {
      delete [] (static_cast<::CAP::JetSingleHistos*>(p));
   }
   static void destruct_CAPcLcLJetSingleHistos(void *p) {
      typedef ::CAP::JetSingleHistos current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CAP::JetSingleHistos

namespace CAP {
//______________________________________________________________________________
void JetPairHistos::Streamer(TBuffer &R__b)
{
   // Stream an object of class CAP::JetPairHistos.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(CAP::JetPairHistos::Class(),this);
   } else {
      R__b.WriteClassBuffer(CAP::JetPairHistos::Class(),this);
   }
}

} // namespace CAP
namespace ROOT {
   // Wrappers around operator new
   static void *new_CAPcLcLJetPairHistos(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetPairHistos : new ::CAP::JetPairHistos;
   }
   static void *newArray_CAPcLcLJetPairHistos(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetPairHistos[nElements] : new ::CAP::JetPairHistos[nElements];
   }
   // Wrapper around operator delete
   static void delete_CAPcLcLJetPairHistos(void *p) {
      delete (static_cast<::CAP::JetPairHistos*>(p));
   }
   static void deleteArray_CAPcLcLJetPairHistos(void *p) {
      delete [] (static_cast<::CAP::JetPairHistos*>(p));
   }
   static void destruct_CAPcLcLJetPairHistos(void *p) {
      typedef ::CAP::JetPairHistos current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CAP::JetPairHistos

namespace CAP {
//______________________________________________________________________________
void JetAnalyzer::Streamer(TBuffer &R__b)
{
   // Stream an object of class CAP::JetAnalyzer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(CAP::JetAnalyzer::Class(),this);
   } else {
      R__b.WriteClassBuffer(CAP::JetAnalyzer::Class(),this);
   }
}

} // namespace CAP
namespace ROOT {
   // Wrappers around operator new
   static void *new_CAPcLcLJetAnalyzer(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetAnalyzer : new ::CAP::JetAnalyzer;
   }
   static void *newArray_CAPcLcLJetAnalyzer(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::CAP::JetAnalyzer[nElements] : new ::CAP::JetAnalyzer[nElements];
   }
   // Wrapper around operator delete
   static void delete_CAPcLcLJetAnalyzer(void *p) {
      delete (static_cast<::CAP::JetAnalyzer*>(p));
   }
   static void deleteArray_CAPcLcLJetAnalyzer(void *p) {
      delete [] (static_cast<::CAP::JetAnalyzer*>(p));
   }
   static void destruct_CAPcLcLJetAnalyzer(void *p) {
      typedef ::CAP::JetAnalyzer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CAP::JetAnalyzer

namespace {
  void TriggerDictionaryInitialization_libJets_Impl() {
    static const char* headers[] = {
"JetHistos.hpp",
"JetSingleHistos.hpp",
"JetPairHistos.hpp",
"JetAnalyzer.hpp",
nullptr
    };
    static const char* includePaths[] = {
"/usr/local/Cellar/root/6.32.06/include/root",
"/Users/aa7526/Documents/GitHub/CAP6.1/src",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Base",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Math",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Xml",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/ParticleDb",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Particles",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/SubSample",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/CAPPythia",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Global",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Jets",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Spherocity",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/ParticleSingle",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/ParticlePair",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/ParticlePair3D",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/NuDyn",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/PtPt",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Plotting",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Therminator",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Performance",
"/Users/aa7526/Documents/GitHub/CAP6.1/src/Exec",
"/include",
"/Users/aa7526/opt/fastjet-3.4.3/include",
"/Users/aa7526/opt/Pythia/pythia8307/include",
"/Users/aa7526/opt/Pythia/pythia8307/include/Pythia8",
"/usr/local/Cellar/root/6.32.06/include/root",
"/Users/aa7526/Documents/GitHub/CAP6.1/homeBuild/Jets/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libJets dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace CAP{class __attribute__((annotate("$clingAutoload$JetHistos.hpp")))  JetHistos;}
namespace CAP{class __attribute__((annotate("$clingAutoload$JetSingleHistos.hpp")))  JetSingleHistos;}
namespace CAP{class __attribute__((annotate("$clingAutoload$JetPairHistos.hpp")))  JetPairHistos;}
namespace CAP{class __attribute__((annotate("$clingAutoload$JetAnalyzer.hpp")))  JetAnalyzer;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libJets dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "JetHistos.hpp"
#include "JetSingleHistos.hpp"
#include "JetPairHistos.hpp"
#include "JetAnalyzer.hpp"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CAP::JetAnalyzer", payloadCode, "@",
"CAP::JetHistos", payloadCode, "@",
"CAP::JetPairHistos", payloadCode, "@",
"CAP::JetSingleHistos", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libJets",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libJets_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libJets_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libJets() {
  TriggerDictionaryInitialization_libJets_Impl();
}
