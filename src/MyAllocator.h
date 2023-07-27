#ifndef MY_ALLOCATOR_H_
#define MY_ALLOCATOR_H_

#include "my_cuda_utils.h"
#include <limits>
#include <iostream>

   template <class T>
   class MyAlloc {
     public:
       // type definitions
       typedef T        value_type;
       typedef T*       pointer;
       typedef const T* const_pointer;
       typedef T&       reference;
       typedef const T& const_reference;
       typedef std::size_t    size_type;
       typedef std::ptrdiff_t difference_type;

       // rebind allocator to type U
       template <class U>
       struct rebind {
           typedef MyAlloc<U> other;
       };

       // return address of values
       pointer address (reference value) const {
           return &value;
       }
       const_pointer address (const_reference value) const {
           return &value;
       }

       /* constructors and destructor
        * - nothing to do because the allocator has no state
        */
       MyAlloc() throw() {
       }
       MyAlloc(const MyAlloc&) throw() {
       }
       template <class U>
         MyAlloc (const MyAlloc<U>&) throw() {
       }
       ~MyAlloc() throw() {
       }

       // return maximum number of elements that can be allocated
       size_type max_size () const throw() {
           return std::numeric_limits<std::size_t>::max() / sizeof(T); //TODO: Ask Khaled if I can assume it
       }

       // allocate but don't initialize num elements of type T
       pointer allocate (size_type num, const void* = 0) {
           // print message and allocate memory with global new
           //std::cerr << "allocate " << num << " element(s)"
                //     << " of size " << sizeof(T) << std::endl;
           
	   if (0 == s)
		   return NULL;

	   pointer temp;
	   cudaMallocHost(reinterpret_cast<void**>&temp, s * sizeof(T));
	   cuda_check_last(__FILE__,__LINE__);

	   if (temp == NULL)
  		throw std::bad_alloc();

	   return temp;		
	       
           
	   
	   
	   //std::cerr << " allocated at: " << (void*)ret << std::endl;
           return ret;
       }

       // initialize elements of allocated storage p with value value
       void construct (pointer p, const T& value) {
           // initialize memory with placement new
           new((void*)p)T(value);
       }

       // destroy elements of initialized storage p
       void destroy (pointer p) {
           // destroy objects by calling their destructor
           p->~T();
       }

       // deallocate storage p of deleted elements
       void deallocate (pointer p, size_type num) {
       		cudaFreeHost(p);
       }
   };


#endif
