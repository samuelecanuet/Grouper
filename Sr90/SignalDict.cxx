// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME SignalDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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
#include "Signal.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_Signal(void *p = nullptr);
   static void *newArray_Signal(Long_t size, void *p);
   static void delete_Signal(void *p);
   static void deleteArray_Signal(void *p);
   static void destruct_Signal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Signal*)
   {
      ::Signal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Signal >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Signal", ::Signal::Class_Version(), "Signal.h", 7,
                  typeid(::Signal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Signal::Dictionary, isa_proxy, 4,
                  sizeof(::Signal) );
      instance.SetNew(&new_Signal);
      instance.SetNewArray(&newArray_Signal);
      instance.SetDelete(&delete_Signal);
      instance.SetDeleteArray(&deleteArray_Signal);
      instance.SetDestructor(&destruct_Signal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Signal*)
   {
      return GenerateInitInstanceLocal((::Signal*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Signal*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Signal::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Signal::Class_Name()
{
   return "Signal";
}

//______________________________________________________________________________
const char *Signal::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Signal::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Signal::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Signal::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Signal::Streamer(TBuffer &R__b)
{
   // Stream an object of class Signal.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Signal::Class(),this);
   } else {
      R__b.WriteClassBuffer(Signal::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Signal(void *p) {
      return  p ? new(p) ::Signal : new ::Signal;
   }
   static void *newArray_Signal(Long_t nElements, void *p) {
      return p ? new(p) ::Signal[nElements] : new ::Signal[nElements];
   }
   // Wrapper around operator delete
   static void delete_Signal(void *p) {
      delete ((::Signal*)p);
   }
   static void deleteArray_Signal(void *p) {
      delete [] ((::Signal*)p);
   }
   static void destruct_Signal(void *p) {
      typedef ::Signal current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Signal

namespace ROOT {
   static TClass *vectorlESignalgR_Dictionary();
   static void vectorlESignalgR_TClassManip(TClass*);
   static void *new_vectorlESignalgR(void *p = nullptr);
   static void *newArray_vectorlESignalgR(Long_t size, void *p);
   static void delete_vectorlESignalgR(void *p);
   static void deleteArray_vectorlESignalgR(void *p);
   static void destruct_vectorlESignalgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Signal>*)
   {
      vector<Signal> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Signal>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Signal>", -2, "vector", 339,
                  typeid(vector<Signal>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlESignalgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Signal>) );
      instance.SetNew(&new_vectorlESignalgR);
      instance.SetNewArray(&newArray_vectorlESignalgR);
      instance.SetDelete(&delete_vectorlESignalgR);
      instance.SetDeleteArray(&deleteArray_vectorlESignalgR);
      instance.SetDestructor(&destruct_vectorlESignalgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Signal> >()));

      ::ROOT::AddClassAlternate("vector<Signal>","std::vector<Signal, std::allocator<Signal> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Signal>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlESignalgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Signal>*)nullptr)->GetClass();
      vectorlESignalgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlESignalgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlESignalgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Signal> : new vector<Signal>;
   }
   static void *newArray_vectorlESignalgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Signal>[nElements] : new vector<Signal>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlESignalgR(void *p) {
      delete ((vector<Signal>*)p);
   }
   static void deleteArray_vectorlESignalgR(void *p) {
      delete [] ((vector<Signal>*)p);
   }
   static void destruct_vectorlESignalgR(void *p) {
      typedef vector<Signal> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Signal>

namespace {
  void TriggerDictionaryInitialization_SignalDict_Impl() {
    static const char* headers[] = {
"Signal.h",
nullptr
    };
    static const char* includePaths[] = {
"/softs/root/6.26.10/include/",
"/home/local1/Documents/2021_Analysis/Grouper/Sr90/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "SignalDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate(R"ATTRDUMP(Signal)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$Signal.h")))  Signal;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "SignalDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Signal.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Signal", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("SignalDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_SignalDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_SignalDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_SignalDict() {
  TriggerDictionaryInitialization_SignalDict_Impl();
}
