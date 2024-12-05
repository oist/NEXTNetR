#pragma once

#include <cstddef>      // for nullptr_t, NULL
#include <memory>       // for bad_weak_ptr
#include <type_traits>  // for add_lvalue_reference
#include <string>

#include <boost/core/demangle.hpp>

#include "cpp11/R.hpp"         // for SEXP, SEXPREC, TYPEOF, R_NilValue, R_C...
#include "cpp11/protect.hpp"   // for protect, safe, protect::function
#include "cpp11/r_bool.hpp"    // for r_bool
#include "cpp11/r_vector.hpp"  // for type_error
#include "cpp11/sexp.hpp"      // for sexp

namespace nextnetR {


template <typename T>
void default_deleter(T* obj) {
  delete obj;
}

class tag_error : public std::exception {
 public:
  tag_error(SEXP expected, SEXP actual) : expected_(expected), actual_(actual) {}
  virtual const char* what() const noexcept override {
    snprintf(str_, 128, "object tag '%s' invalid, expected '%s'",
             CHAR(PRINTNAME(actual_)), CHAR(PRINTNAME(expected_)));
    return str_;
  }

 private:
  SEXP expected_;
  SEXP actual_;
  mutable char str_[128];
};

template <typename T, void Deleter(T*) = default_deleter<T>>
class safe_external_pointer {
 private:
  static SEXP tag_; // Stores a symbol (SYMSXP) so don't worry about lifetimes

  cpp11::sexp data_ = nullptr;

  static SEXP tag() {
    if (tag_ == nullptr) {
      const std::string name = boost::core::demangle(typeid(T).name());
      tag_ = cpp11::safe[Rf_install](name.c_str());
    }
    return tag_;
  }

  static SEXP valid_data(SEXP data) {
    if (data == nullptr) {
      throw cpp11::type_error(EXTPTRSXP, NILSXP);
    }
    if (TYPEOF(data) != EXTPTRSXP) {
      throw cpp11::type_error(EXTPTRSXP, TYPEOF(data));
    }
    const SEXP data_tag = R_ExternalPtrTag(data);
    if (TYPEOF(data_tag) != SYMSXP) {
      throw cpp11::type_error(SYMSXP, TYPEOF(data_tag));
    }
    if (data_tag != tag())
      throw tag_error(tag(), data_tag);

    return data;
  }

  static void r_deleter(SEXP data) {
    if (data == nullptr)  return;
    if (TYPEOF(data) != EXTPTRSXP) return;

    const SEXP data_tag = R_ExternalPtrTag(data);
    if (TYPEOF(data_tag) != SYMSXP) return;
    if (data_tag != tag()) return;

    T* ptr = static_cast<T*>(R_ExternalPtrAddr(data));
    if (ptr == NULL) return;

    R_ClearExternalPtr(data);

    Deleter(ptr);
  }

 public:
  using pointer = T*;

  safe_external_pointer() noexcept {}
  safe_external_pointer(std::nullptr_t) noexcept {}

  safe_external_pointer(SEXP data) : data_(valid_data(data)) {}

  safe_external_pointer(pointer p, bool use_deleter = true, bool finalize_on_exit = true)
      : data_(cpp11::safe[R_MakeExternalPtr]((void*)p, tag(), R_NilValue))
  {
    if (use_deleter) {
      R_RegisterCFinalizerEx(data_, r_deleter, static_cast<cpp11::r_bool>(finalize_on_exit));
    }
  }

  safe_external_pointer(pointer p, cpp11::sexp prot_data, bool use_deleter = true, bool finalize_on_exit = true)
      : data_(cpp11::safe[R_MakeExternalPtr]((void*)p, tag(), prot_data))
  {
    if (use_deleter) {
      R_RegisterCFinalizerEx(data_, r_deleter, static_cast<cpp11::r_bool>(finalize_on_exit));
    }
  }

  safe_external_pointer(const safe_external_pointer& rhs)
    :data_(rhs.data_)
  {}

  safe_external_pointer(safe_external_pointer&& rhs)
    :data_(std::move(rhs.data_))
  {}

  safe_external_pointer& operator=(safe_external_pointer& rhs) noexcept {
    data_ = rhs.data_;
  }

  safe_external_pointer& operator=(safe_external_pointer&& rhs) noexcept {
    data_ = std::move(rhs.data_);
  }

  safe_external_pointer& operator=(std::nullptr_t) noexcept {
    data_ = nullptr;
  };

  operator SEXP() const noexcept { return data_; }

  pointer get() const noexcept {
    return static_cast<T*>(R_ExternalPtrAddr(data_));
  }

  SEXP protected_data() const noexcept {
    return R_ExternalPtrProtected(data_);
  }

  typename std::add_lvalue_reference<T>::type operator*() {
    pointer addr = get();
    if (addr == nullptr) {
      throw std::bad_weak_ptr();
    }
    return *get();
  }

  pointer operator->() const {
    pointer addr = get();
    if (addr == nullptr) {
      throw std::bad_weak_ptr();
    }
    return get();
  }

  void reset(pointer ptr = pointer()) {
    *this = safe_external_pointer(ptr);
  }

  void reset(pointer ptr = pointer(), cpp11::sexp prot_data = cpp11::sexp()) {
    *this = safe_external_pointer(ptr, prot_data);
  }

  void swap(safe_external_pointer& other) noexcept {
    SEXP tmp = other.data_;
    other.data_ = data_;
    data_ = tmp;
  }

  operator bool() noexcept { return data_ != nullptr; }
};

template <typename T, void Deleter(T*)>
SEXP safe_external_pointer<T, Deleter>::tag_ = nullptr;

template <class T, void Deleter(T*)>
void swap(safe_external_pointer<T, Deleter>& lhs, safe_external_pointer<T, Deleter>& rhs) noexcept {
  lhs.swap(rhs);
}

template <class T, void Deleter(T*)>
bool operator==(const safe_external_pointer<T, Deleter>& x,
                const safe_external_pointer<T, Deleter>& y) {
  return x.data_ == y.data_;
}

template <class T, void Deleter(T*)>
bool operator!=(const safe_external_pointer<T, Deleter>& x,
                const safe_external_pointer<T, Deleter>& y) {
  return x.data_ != y.data_;
}

template <class T, void Deleter(T*)>
bool operator<(const safe_external_pointer<T, Deleter>& x,
               const safe_external_pointer<T, Deleter>& y) {
  return x.data_ < y.data_;
}

template <class T, void Deleter(T*)>
bool operator<=(const safe_external_pointer<T, Deleter>& x,
                const safe_external_pointer<T, Deleter>& y) {
  return x.data_ <= y.data_;
}

template <class T, void Deleter(T*)>
bool operator>(const safe_external_pointer<T, Deleter>& x,
               const safe_external_pointer<T, Deleter>& y) {
  return x.data_ > y.data_;
}

template <class T, void Deleter(T*)>
bool operator>=(const safe_external_pointer<T, Deleter>& x,
                const safe_external_pointer<T, Deleter>& y) {
  return x.data_ >= y.data_;
}

}  // namespace cpp11
