//
//  DerivedFinalClass.hpp
//  testRTTypes
//
//  Created by Omar Duran on 12/16/24.
//

#ifndef DerivedFinalClass_hpp
#define DerivedFinalClass_hpp

#include <stdio.h>
#include "BaseClass.hpp"

class DerivedFinalClass final : public BaseClass
{
public:
    explicit DerivedFinalClass();
    
    virtual ~DerivedFinalClass() noexcept override {
        
    };
};

#endif /* DerivedFinalClass_hpp */
